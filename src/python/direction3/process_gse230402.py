#!/usr/bin/env python3
"""
Process GSE230402 (Sex differences in CR on hepatic gene expression) through tAge pipeline.
22 samples: WT mice, ad libitum vs 30% caloric restriction.
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
DATA_FILE = PROJECT_ROOT / "data" / "geo_expansion" / "counts" / "GSE230402_Counts.normalized_WT.txt.gz"
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
        "feature_names": feature_names,
        "coef": estimator.coef_,
        "intercept": float(estimator.intercept_),
        "impute_stats": imputer.statistics_,
        "center_mean": scaler.mean_,
    }

def map_symbols_to_entrez(symbols, species='mouse', batch_size=1000):
    try:
        import mygene
    except ImportError:
        print("ERROR: mygene not installed")
        return {}
    mg = mygene.MyGeneInfo()
    mapping = {}
    for i in range(0, len(symbols), batch_size):
        batch = symbols[i:i+batch_size]
        result = mg.querymany(batch, scopes='symbol,alias', fields='entrezgene',
                              species=species, verbose=False, as_dataframe=True)
        for idx, row in result.iterrows():
            if 'entrezgene' in row and pd.notna(row['entrezgene']):
                mapping[idx] = int(row['entrezgene'])
    return mapping

def main():
    print("="*60)
    print("Processing GSE230402")
    print("="*60)

    # Load normalized counts
    print("\nLoading counts...")
    df = pd.read_csv(DATA_FILE, sep='\t', compression='gzip')
    print(f"  Shape: {df.shape}")
    genes = df.iloc[:, 0].astype(str).tolist()
    counts = df.iloc[:, 1:].values.astype(float)
    samples = df.columns[1:].tolist()
    print(f"  Samples: {samples}")

    # Real metadata from GEO series matrix (interleaved order)
    meta = pd.DataFrame({
        "sample_id": samples,
        "sex": ["Female","Male","Male","Female","Female","Female","Male","Female","Male","Female","Male","Male","Male","Female","Male","Male","Female","Female","Female","Female","Male","Male"],
        "diet":  ["CR","AL","CR","AL","CR","AL","AL","CR","CR","CR","AL","AL","CR","AL","CR","AL","CR","AL","AL","CR","AL","CR"],
    })
    print("\nSample design (inferred):")
    print(meta.groupby(["sex", "diet"]).size().to_string())

    # Map gene symbols to Entrez
    print("\nMapping Gene Symbols -> Entrez...")
    cache = RESULTS_DIR / "gse230402_symbol_to_entrez.pkl"
    if cache.exists():
        import pickle
        mapping = pickle.load(open(cache, "rb"))
        print(f"  Loaded cache: {len(mapping)} mappings")
    else:
        mapping = map_symbols_to_entrez(genes)
        import pickle
        pickle.dump(mapping, open(cache, "wb"))
        print(f"  Mapped: {len(mapping)}/{len(genes)}")

    entrez_list = []
    expr_list = []
    for i, g in enumerate(genes):
        if g in mapping:
            entrez_list.append(str(mapping[g]))
            expr_list.append(counts[i])
    expr = pd.DataFrame(np.array(expr_list), columns=samples, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    # Log2 transform (normalized counts)
    log_expr = np.log2(expr + 1)

    # Predict tAge
    print("\nRunning tAge prediction...")
    mi = load_model_info()
    n_samples = log_expr.shape[1]
    n_features = len(mi["feature_names"])
    X = np.zeros((n_samples, n_features), dtype=float)
    f2i = {name: i for i, name in enumerate(mi["feature_names"])}
    common = 0
    for gene in log_expr.index:
        if gene in f2i:
            X[:, f2i[gene]] = log_expr.loc[gene].values
            common += 1
    print(f"  Common genes: {common}/{n_features}")
    for i in range(n_features):
        if np.all(X[:, i] == 0):
            X[:, i] = mi["impute_stats"][i]
    Xc = X - mi["center_mean"][np.newaxis, :]
    preds = mi["intercept"] + np.dot(Xc, mi["coef"])
    preds_adj = preds * 48.0

    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    meta.to_csv(RESULTS_DIR / "gse230402_predictions.csv", index=False)

    # Summarize
    print("\nMean tAge by group:")
    print(meta.groupby(["sex", "diet"])["tAge_adj"].agg(["mean", "std", "count"]).to_string())

    # CR - AL differences
    print("\nCR - AL differences (months, negative=rejuvenation):")
    diff_rows = []
    for sex in ["Male", "Female"]:
        sub = meta[meta["sex"] == sex]
        al = sub[sub["diet"] == "AL"]["tAge_adj"]
        cr = sub[sub["diet"] == "CR"]["tAge_adj"]
        diff = cr.mean() - al.mean()
        diff_rows.append({"sex": sex, "n_cr": len(cr), "n_al": len(al), "cr_mean": cr.mean(), "al_mean": al.mean(), "difference": diff})
        print(f"  {sex}: Δ = {diff:+.1f} months (CR {cr.mean():.1f} vs AL {al.mean():.1f})")
    pd.DataFrame(diff_rows).to_csv(RESULTS_DIR / "gse230402_cr_al_differences.csv", index=False)

    # Plot
    plt.figure(figsize=(6, 5))
    sns.boxplot(data=meta, x="diet", y="tAge_adj", hue="sex")
    plt.title("GSE230402 Liver — tAge by Diet")
    plt.ylabel("Predicted tAge (months)")
    plt.tight_layout()
    out = FIGURES_DIR / "gse230402_liver_tage.png"
    plt.savefig(out)
    plt.close()
    print(f"\nSaved {out}")
    print("Done!")

if __name__ == "__main__":
    main()
