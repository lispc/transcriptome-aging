#!/usr/bin/env python3
"""
Process GSE288795 (Rapamycin + Trametinib combo in mouse) through tAge pipeline.
This is a high-value validation dataset while we wait for ARCHS4.
"""
import os, re
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
import warnings
from sklearn.exceptions import InconsistentVersionWarning
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings("ignore", category=UserWarning)

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
MODEL_PATH = PROJECT_ROOT / "models" / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
DATA_FILE = PROJECT_ROOT / "data" / "geo_expansion" / "counts" / "GSE288795_count_data.xlsx"
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
FIGURES_DIR = PROJECT_ROOT / "figures" / "direction3"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 150
sns.set_style("whitegrid")

def load_model_info():
    model = joblib.load(MODEL_PATH)
    imputer = model.named_steps["imputation"]
    scaler = model.named_steps["scaler"]
    estimator = model.named_steps["estimator"]
    feature_names = list(model.feature_names_in_)
    return {
        "model": model,
        "feature_names": feature_names,
        "coef": estimator.coef_,
        "intercept": float(estimator.intercept_),
        "impute_stats": imputer.statistics_,
        "center_mean": scaler.mean_,
    }

def parse_meta(desc):
    """Parse 'Male Muscle Rapamycin Biol Rep 1' into dict."""
    m = re.match(r'(Male|Female)\s+(\w+)\s+(Rapamycin/Trametinib|Rapamycin|Trametinib|Control)\s+Biol Rep\s+(\d+)', desc)
    if not m:
        return {"sex": "", "tissue": "", "treatment": "", "rep": ""}
    return {
        "sex": m.group(1),
        "tissue": m.group(2),
        "treatment": m.group(3),
        "rep": int(m.group(4)),
    }

def map_ensembl_to_entrez(ensembl_ids, species='mouse', batch_size=1000):
    try:
        import mygene
    except ImportError:
        print("ERROR: mygene not installed")
        return {}
    mg = mygene.MyGeneInfo()
    clean = [g.split('.')[0] for g in ensembl_ids]
    mapping = {}
    for i in range(0, len(clean), batch_size):
        batch = clean[i:i+batch_size]
        result = mg.querymany(batch, scopes='ensembl.gene', fields='entrezgene',
                              species=species, verbose=False, as_dataframe=True)
        for idx, row in result.iterrows():
            if 'entrezgene' in row and pd.notna(row['entrezgene']):
                mapping[idx] = int(row['entrezgene'])
    return mapping

def main():
    print("="*60)
    print("Processing GSE288795")
    print("="*60)

    # Load data
    print("\nLoading Excel...")
    df_raw = pd.read_excel(DATA_FILE, sheet_name=0, header=None)
    print(f"  Raw shape: {df_raw.shape}")

    # Row 0 = descriptions, Row 1 = library names, Row 2+ = data
    descriptions = df_raw.iloc[0, 1:].tolist()
    sample_ids = [f"S{i+1}" for i in range(len(descriptions))]
    data = df_raw.iloc[2:, 1:].values.astype(float)
    genes = df_raw.iloc[2:, 0].astype(str).tolist()
    print(f"  Samples: {len(descriptions)}")
    print(f"  Genes: {len(genes)}")

    # Build metadata
    meta = pd.DataFrame([parse_meta(d) for d in descriptions])
    meta["sample_id"] = sample_ids
    meta["description"] = descriptions
    print("\nSample breakdown:")
    print(meta.groupby(["tissue", "treatment", "sex"]).size().to_string())

    # Map genes
    print("\nMapping Ensembl -> Entrez...")
    cache = RESULTS_DIR / "gse288795_ensembl_to_entrez.pkl"
    if cache.exists():
        import pickle
        mapping = pickle.load(open(cache, "rb"))
        print(f"  Loaded cache: {len(mapping)} mappings")
    else:
        mapping = map_ensembl_to_entrez(genes)
        import pickle
        pickle.dump(mapping, open(cache, "wb"))
        print(f"  Mapped: {len(mapping)}/{len(genes)}")

    # Build expression matrix with Entrez IDs
    entrez_list = []
    expr_list = []
    for i, g in enumerate(genes):
        cg = g.split('.')[0]
        if cg in mapping:
            entrez_list.append(str(mapping[cg]))
            expr_list.append(data[i])
    expr = pd.DataFrame(np.array(expr_list), columns=sample_ids, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    # Log2 transform (normalized counts, so use log2(x+1))
    log_expr = np.log2(expr + 1)

    # tAge prediction
    print("\nRunning tAge prediction...")
    model_info = load_model_info()
    feature_names = model_info["feature_names"]
    n_samples = log_expr.shape[1]
    n_features = len(feature_names)
    X = np.zeros((n_samples, n_features), dtype=float)
    feat_to_idx = {name: i for i, name in enumerate(feature_names)}

    common = 0
    for gene in log_expr.index:
        if gene in feat_to_idx:
            X[:, feat_to_idx[gene]] = log_expr.loc[gene].values
            common += 1
    print(f"  Common genes: {common}/{n_features}")

    for i in range(n_features):
        if np.all(X[:, i] == 0):
            X[:, i] = model_info["impute_stats"][i]

    X_centered = X - model_info["center_mean"][np.newaxis, :]
    preds = model_info["intercept"] + np.dot(X_centered, model_info["coef"])
    preds_adj = preds * 48.0  # mouse adjustment

    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj

    # Save predictions
    meta.to_csv(RESULTS_DIR / "gse288795_predictions.csv", index=False)
    print(f"\nPredictions saved.")

    # Summary by group
    print("\nMean tAge by group:")
    summary = meta.groupby(["tissue", "treatment"])["tAge_adj"].agg(["mean", "std", "count"])
    print(summary.to_string())

    # Drug - Control differences (within tissue/sex)
    print("\nDrug - Control differences (months, negative=rejuvenation):")
    diff_rows = []
    for tissue in meta["tissue"].unique():
        for sex in meta["sex"].unique():
            sub = meta[(meta["tissue"] == tissue) & (meta["sex"] == sex)]
            if len(sub) == 0:
                continue
            ctrl = sub[sub["treatment"] == "Control"]["tAge_adj"]
            for drug in ["Rapamycin", "Trametinib", "Rapamycin/Trametinib"]:
                drug_vals = sub[sub["treatment"] == drug]["tAge_adj"]
                if len(drug_vals) == 0 or len(ctrl) == 0:
                    continue
                diff = drug_vals.mean() - ctrl.mean()
                diff_rows.append({
                    "tissue": tissue,
                    "sex": sex,
                    "drug": drug,
                    "n_drug": len(drug_vals),
                    "n_ctrl": len(ctrl),
                    "drug_mean": drug_vals.mean(),
                    "ctrl_mean": ctrl.mean(),
                    "difference": diff,
                })
    diff_df = pd.DataFrame(diff_rows)
    print(diff_df.to_string(index=False))
    diff_df.to_csv(RESULTS_DIR / "gse288795_drug_control_differences.csv", index=False)

    # Plotting
    print("\nGenerating figures...")
    for tissue in meta["tissue"].unique():
        sub = meta[meta["tissue"] == tissue]
        plt.figure(figsize=(8, 5))
        sns.boxplot(data=sub, x="treatment", y="tAge_adj", hue="sex")
        plt.title(f"GSE288795 {tissue} — tAge by Treatment")
        plt.ylabel("Predicted tAge (months)")
        plt.axhline(sub[sub["treatment"]=="Control"]["tAge_adj"].mean(), color="gray", linestyle="--")
        plt.tight_layout()
        out = FIGURES_DIR / f"gse288795_{tissue.lower()}_tage.png"
        plt.savefig(out)
        plt.close()
        print(f"  Saved {out}")

    # Overall drug effect barplot
    plt.figure(figsize=(8, 5))
    pivot = diff_df.groupby("drug")["difference"].mean().reindex(["Rapamycin", "Trametinib", "Rapamycin/Trametinib"])
    colors = ["green" if v < 0 else "red" for v in pivot.values]
    pivot.plot(kind="bar", color=colors)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.ylabel("Mean Drug - Control tAge (months)")
    plt.title("GSE288795 Drug Effects (pooled across tissues/sex)")
    plt.xticks(rotation=0)
    plt.tight_layout()
    out = FIGURES_DIR / "gse288795_drug_effects_pooled.png"
    plt.savefig(out)
    plt.close()
    print(f"  Saved {out}")

    print("\nDone!")

if __name__ == "__main__":
    main()
