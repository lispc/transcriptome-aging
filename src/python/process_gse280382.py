#!/usr/bin/env python3
"""
Process GSE280382: GLP-1R agonist (exenatide) aging rejuvenation study.
AgedLT: 30 weeks, multi-tissue (Young_ctrl, Aged_ctrl, Aged_exenatide)
AgedST: 13 weeks (Aged_ctrl, Aged_exenatide, Aged_rapamycin, plus KO groups)
Young: young controls with/without exenatide
"""

import sys
import gzip
import pickle
import joblib
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

PROJECT_ROOT = Path(__file__).parent.parent.parent
MODEL_PATH = PROJECT_ROOT / "models" / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
FIGURES_DIR = PROJECT_ROOT / "figures" / "gse280382"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Load model
print("Loading model...")
model = joblib.load(MODEL_PATH)
imputer = model.named_steps["imputation"]
scaler = model.named_steps["scaler"]
estimator = model.named_steps["estimator"]
coef = estimator.coef_
intercept = float(estimator.intercept_)
model_features = list(model.feature_names_in_)
impute_stats = imputer.statistics_
center_mean = scaler.mean_

# mygene for symbol->Entrez
print("Setting up mygene...")
import mygene
mg = mygene.MyGeneInfo()

# ============================================
# Helper: map gene symbols to Entrez, aggregate
# ============================================
def map_symbols_to_entrez(symbols):
    """Map mouse gene symbols to Entrez Gene IDs, aggregate duplicates by mean."""
    print(f"  Querying mygene for {len(symbols)} symbols...")
    # Query in batches
    batch_size = 1000
    all_results = []
    for i in range(0, len(symbols), batch_size):
        batch = symbols[i:i+batch_size]
        res = mg.querymany(batch, scopes="symbol", species="mouse",
                           fields="entrezgene", as_dataframe=True)
        all_results.append(res)
    res_df = pd.concat(all_results)
    
    # Build mapping
    mapping = {}
    for sym in symbols:
        if sym in res_df.index:
            row = res_df.loc[sym]
            if isinstance(row, pd.DataFrame):
                # Multiple matches, take first with entrezgene
                for _, r in row.iterrows():
                    eg = r.get("entrezgene")
                    if pd.notna(eg):
                        mapping[sym] = str(int(eg))
                        break
            else:
                eg = row.get("entrezgene")
                if pd.notna(eg):
                    mapping[sym] = str(int(eg))
    
    print(f"  Mapped {len(mapping)}/{len(symbols)} symbols to Entrez")
    return mapping

# ============================================
# Helper: predict tAge from counts DataFrame
# ============================================
def predict_tage(counts_df, mapping, dataset_name):
    """
    counts_df: index=symbols, columns=samples
    Returns: DataFrame with tAge, tAge_adj per sample
    """
    # Filter to mapped genes
    mapped = counts_df.index.intersection(list(mapping.keys()))
    print(f"  {len(mapped)} genes have Entrez mapping")
    
    # Build Entrez-indexed matrix
    entrez_list = [mapping[s] for s in mapped]
    expr = counts_df.loc[mapped].copy()
    expr.index = entrez_list
    
    # Aggregate duplicates by mean
    expr = expr.groupby(level=0).mean()
    
    # Align to model features
    expr_aligned = pd.DataFrame(index=model_features, columns=expr.columns)
    common = expr.index.intersection(model_features)
    print(f"  {len(common)} features overlap with model")
    expr_aligned.loc[common] = expr.loc[common].astype(float)
    expr_aligned = expr_aligned.fillna(0.0)
    
    # TMM normalization via edgeR (using Rscript subprocess for stability)
    print("  Running edgeR TMM normalization via Rscript...")
    import subprocess
    import tempfile
    import os
    
    # Ensure float dtype
    expr_aligned = expr_aligned.astype(float)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = os.path.join(tmpdir, "counts.csv")
        output_csv = os.path.join(tmpdir, "logcpm.csv")
        expr_aligned.to_csv(input_csv)
        
        r_script = str(PROJECT_ROOT / "src" / "r" / "tmm_normalize.R")
        result = subprocess.run(
            ["Rscript", r_script, input_csv, output_csv],
            capture_output=True, text=True, check=True
        )
        print(result.stdout.strip())
        if result.stderr:
            print(result.stderr.strip())
        
        logcpm_df = pd.read_csv(output_csv, index_col=0)
        # CRITICAL: Ensure index and columns are strings for proper alignment
        logcpm_df.index = logcpm_df.index.astype(str)
        logcpm_df.columns = logcpm_df.columns.astype(str)
        # Align to model features and samples
        logcpm_df = logcpm_df.reindex(index=model_features, columns=expr_aligned.columns)
    
    # Impute (manual to avoid sklearn version incompatibility with pickled model)
    X_arr = logcpm_df.values.T  # samples x genes
    X_imp = np.where(np.isnan(X_arr), impute_stats, X_arr)
    
    # Center (manual to match model)
    X_centered = X_imp - center_mean
    
    # Predict
    y = intercept + X_centered @ coef
    y_adj = (y - 5.5) * 48
    
    preds = pd.DataFrame({
        "sample_id": expr_aligned.columns,
        "tAge": y,
        "tAge_adj": y_adj,
    })
    preds["dataset"] = dataset_name
    return preds


# ============================================
# Process AgedLT (Long Term)
# ============================================
def process_agedlt():
    print("\n=== Processing AgedLT ===")
    counts_path = PROJECT_ROOT / "data" / "geo" / "GSE280382_AgedLT_counts.csv.gz"
    meta_path = PROJECT_ROOT / "data" / "geo" / "GSE280382_AgedLT_metadata.csv.gz"
    
    print("  Reading counts...")
    counts = pd.read_csv(counts_path, index_col=0)
    print(f"  Shape: {counts.shape}")
    
    meta = pd.read_csv(meta_path)
    print(f"  Metadata: {len(meta)} samples")
    
    # Map symbols
    cache = RESULTS_DIR / "gse280382_symbol_to_entrez.pkl"
    if cache.exists():
        print("  Using cached symbol mapping")
        with open(cache, "rb") as f:
            mapping = pickle.load(f)
    else:
        mapping = map_symbols_to_entrez(list(counts.index))
        with open(cache, "wb") as f:
            pickle.dump(mapping, f)
    
    preds = predict_tage(counts, mapping, "AgedLT")
    
    # Merge with metadata
    meta_dict = meta.set_index("Name").to_dict("index")
    preds["Mouse"] = preds["sample_id"].map(lambda x: meta_dict.get(x, {}).get("Mouse", x.split("_")[0]))
    preds["Tissue"] = preds["sample_id"].map(lambda x: meta_dict.get(x, {}).get("Tissue", "Unknown"))
    preds["Group"] = preds["sample_id"].map(lambda x: meta_dict.get(x, {}).get("Group", "Unknown"))
    preds["Batch"] = preds["sample_id"].map(lambda x: meta_dict.get(x, {}).get("Batch", "Unknown"))
    
    preds.to_csv(RESULTS_DIR / "gse280382_AgedLT_predictions.csv", index=False)
    print(f"  Saved to {RESULTS_DIR / 'gse280382_AgedLT_predictions.csv'}")
    return preds


# ============================================
# Process AgedST (Short Term)
# ============================================
def process_agedst():
    print("\n=== Processing AgedST ===")
    counts_path = PROJECT_ROOT / "data" / "geo" / "GSE280382_AgedST_counts.csv.gz"
    meta_path = PROJECT_ROOT / "data" / "geo" / "GSE280382_AgedST_metadata.csv.gz"
    
    print("  Reading counts...")
    counts = pd.read_csv(counts_path, index_col=0)
    print(f"  Shape: {counts.shape}")
    
    meta = pd.read_csv(meta_path)
    print(f"  Metadata: {len(meta)} samples")
    
    # Use cached mapping
    cache = RESULTS_DIR / "gse280382_symbol_to_entrez.pkl"
    with open(cache, "rb") as f:
        mapping = pickle.load(f)
    
    preds = predict_tage(counts, mapping, "AgedST")
    
    # Merge with metadata (counts columns have tissue suffix, e.g. "E13_Adipose")
    meta_dict = meta.set_index("Animal").to_dict("index")
    preds["Mouse"] = preds["sample_id"].map(lambda x: x.split("_")[0])
    preds["Tissue"] = preds["sample_id"].map(lambda x: "_".join(x.split("_")[1:]) if "_" in x else "Unknown")
    preds["Group"] = preds["sample_id"].map(lambda x: meta_dict.get(x.split("_")[0], {}).get("Group", "Unknown"))
    
    preds.to_csv(RESULTS_DIR / "gse280382_AgedST_predictions.csv", index=False)
    print(f"  Saved to {RESULTS_DIR / 'gse280382_AgedST_predictions.csv'}")
    return preds


# ============================================
# Process Young
# ============================================
def process_young():
    print("\n=== Processing Young ===")
    counts_path = PROJECT_ROOT / "data" / "geo" / "GSE280382_young_counts.csv.gz"
    meta_path = PROJECT_ROOT / "data" / "geo" / "GSE280382_young_metadata.csv.gz"
    
    print("  Reading counts...")
    counts = pd.read_csv(counts_path, index_col=0)
    print(f"  Shape: {counts.shape}")
    
    meta = pd.read_csv(meta_path)
    print(f"  Metadata: {len(meta)} samples")
    
    cache = RESULTS_DIR / "gse280382_symbol_to_entrez.pkl"
    with open(cache, "rb") as f:
        mapping = pickle.load(f)
    
    preds = predict_tage(counts, mapping, "Young")
    
    meta_dict = meta.set_index("Animal").to_dict("index")
    preds["Mouse"] = preds["sample_id"].map(lambda x: x.split("_")[0])
    preds["Tissue"] = preds["sample_id"].map(lambda x: "_".join(x.split("_")[1:]) if "_" in x else "Unknown")
    preds["Group"] = preds["sample_id"].map(lambda x: meta_dict.get(x.split("_")[0], {}).get("Group", "Unknown"))
    
    preds.to_csv(RESULTS_DIR / "gse280382_young_predictions.csv", index=False)
    print(f"  Saved to {RESULTS_DIR / 'gse280382_young_predictions.csv'}")
    return preds


# ============================================
# Generate comparison tables
# ============================================
def generate_comparisons(agedlt, agedst, young):
    print("\n=== Generating Comparisons ===")
    
    # AgedLT: per-tissue group differences
    print("\n-- AgedLT per-tissue differences (Aged_exenatide - Aged_ctrl) --")
    tissues = agedlt["Tissue"].unique()
    agedlt_diffs = []
    for tissue in tissues:
        tdf = agedlt[agedlt["Tissue"] == tissue]
        ctrl = tdf[tdf["Group"] == "Aged_ctrl"]["tAge_adj"]
        exen = tdf[tdf["Group"] == "Aged_exenatide"]["tAge_adj"]
        if len(ctrl) > 0 and len(exen) > 0:
            diff = exen.mean() - ctrl.mean()
            agedlt_diffs.append({
                "tissue": tissue,
                "n_ctrl": len(ctrl),
                "n_exenatide": len(exen),
                "ctrl_mean": ctrl.mean(),
                "exenatide_mean": exen.mean(),
                "difference": diff,
            })
            print(f"  {tissue}: ctrl={ctrl.mean():.1f}, exen={exen.mean():.1f}, diff={diff:.1f} (n={len(ctrl)}/{len(exen)})")
    
    pd.DataFrame(agedlt_diffs).to_csv(
        RESULTS_DIR / "gse280382_AgedLT_differences.csv", index=False
    )
    
    # AgedLT: pooled across tissues
    print("\n-- AgedLT pooled --")
    ctrl_all = agedlt[agedlt["Group"] == "Aged_ctrl"]["tAge_adj"]
    exen_all = agedlt[agedlt["Group"] == "Aged_exenatide"]["tAge_adj"]
    young_all = agedlt[agedlt["Group"] == "Young_ctrl"]["tAge_adj"]
    print(f"  Young_ctrl: {young_all.mean():.1f} (n={len(young_all)})")
    print(f"  Aged_ctrl:  {ctrl_all.mean():.1f} (n={len(ctrl_all)})")
    print(f"  Aged_exen:  {exen_all.mean():.1f} (n={len(exen_all)})")
    print(f"  Aging effect (Aged - Young): {ctrl_all.mean() - young_all.mean():.1f}")
    print(f"  Rejuvenation (Exen - Aged):  {exen_all.mean() - ctrl_all.mean():.1f}")
    
    # AgedST: Exenatide vs Rapamycin
    print("\n-- AgedST: Exenatide vs Rapamycin vs Control --")
    groups = ["Aged_ctrl", "Aged_exenatide", "Aged_rapamycin"]
    for g in groups:
        gdf = agedst[agedst["Group"] == g]["tAge_adj"]
        print(f"  {g}: {gdf.mean():.1f} (n={len(gdf)})")
    
    ctrl_st = agedst[agedst["Group"] == "Aged_ctrl"]["tAge_adj"]
    exen_st = agedst[agedst["Group"] == "Aged_exenatide"]["tAge_adj"]
    rapa_st = agedst[agedst["Group"] == "Aged_rapamycin"]["tAge_adj"]
    
    agedst_diffs = [
        {"group1": "Aged_exenatide", "group2": "Aged_ctrl",
         "n1": len(exen_st), "n2": len(ctrl_st),
         "mean1": exen_st.mean(), "mean2": ctrl_st.mean(),
         "difference": exen_st.mean() - ctrl_st.mean()},
        {"group1": "Aged_rapamycin", "group2": "Aged_ctrl",
         "n1": len(rapa_st), "n2": len(ctrl_st),
         "mean1": rapa_st.mean(), "mean2": ctrl_st.mean(),
         "difference": rapa_st.mean() - ctrl_st.mean()},
        {"group1": "Aged_exenatide", "group2": "Aged_rapamycin",
         "n1": len(exen_st), "n2": len(rapa_st),
         "mean1": exen_st.mean(), "mean2": rapa_st.mean(),
         "difference": exen_st.mean() - rapa_st.mean()},
    ]
    pd.DataFrame(agedst_diffs).to_csv(
        RESULTS_DIR / "gse280382_AgedST_differences.csv", index=False
    )
    
    # AgedST: KO groups
    print("\n-- AgedST: KO effect (hypothalamic GLP-1R KO) --")
    for g in agedst["Group"].unique():
        if "KD" in g:
            gdf = agedst[agedst["Group"] == g]["tAge_adj"]
            print(f"  {g}: {gdf.mean():.1f} (n={len(gdf)})")
    
    # Young: exenatide effect in young mice
    print("\n-- Young: Exenatide effect in young mice --")
    yc = young[young["Group"] == "Young_ctrl"]["tAge_adj"]
    ye = young[young["Group"] == "Young_exenatide"]["tAge_adj"]
    print(f"  Young_ctrl: {yc.mean():.1f} (n={len(yc)})")
    print(f"  Young_exen: {ye.mean():.1f} (n={len(ye)})")
    print(f"  Difference: {ye.mean() - yc.mean():.1f}")


# ============================================
# Visualization
# ============================================
def make_figures(agedlt, agedst, young):
    print("\n=== Making Figures ===")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Figure 1: AgedLT per-tissue boxplot
    fig, ax = plt.subplots(figsize=(16, 8))
    tissues_order = sorted(agedlt["Tissue"].unique())
    plot_data = []
    for tissue in tissues_order:
        tdf = agedlt[agedlt["Tissue"] == tissue]
        for _, row in tdf.iterrows():
            if row["Group"] in ["Aged_ctrl", "Aged_exenatide", "Young_ctrl"]:
                plot_data.append({
                    "Tissue": tissue,
                    "Group": row["Group"],
                    "tAge_adj": row["tAge_adj"],
                })
    pdf = pd.DataFrame(plot_data)
    
    # Relabel
    pdf["Group"] = pdf["Group"].map({
        "Young_ctrl": "Young Ctrl",
        "Aged_ctrl": "Aged Ctrl",
        "Aged_exenatide": "Aged + Exenatide",
    })
    
    sns.boxplot(data=pdf, x="Tissue", y="tAge_adj", hue="Group", ax=ax,
                palette={"Young Ctrl": "#2ecc71", "Aged Ctrl": "#e74c3c", "Aged + Exenatide": "#3498db"})
    ax.set_ylabel("Transcriptomic Age (months)")
    ax.set_xlabel("")
    ax.set_title("GSE280382 AgedLT: GLP-1RA Rejuvenation Across Tissues (30 weeks)")
    ax.tick_params(axis="x", rotation=45)
    ax.axhline(y=0, color="black", linestyle="--", alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "gse280382_AgedLT_tissue_boxplot.png", dpi=150)
    plt.close()
    print(f"  Saved {FIGURES_DIR / 'gse280382_AgedLT_tissue_boxplot.png'}")
    
    # Figure 2: AgedST group comparison (Exenatide vs Rapamycin)
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_data2 = []
    for _, row in agedst.iterrows():
        if row["Group"] in ["Aged_ctrl", "Aged_exenatide", "Aged_rapamycin"]:
            plot_data2.append({
                "Group": row["Group"],
                "tAge_adj": row["tAge_adj"],
            })
    pdf2 = pd.DataFrame(plot_data2)
    pdf2["Group"] = pdf2["Group"].map({
        "Aged_ctrl": "Aged Control",
        "Aged_exenatide": "Aged + Exenatide",
        "Aged_rapamycin": "Aged + Rapamycin",
    })
    
    sns.boxplot(data=pdf2, x="Group", y="tAge_adj", ax=ax,
                palette={"Aged Control": "#e74c3c", "Aged + Exenatide": "#3498db", "Aged + Rapamycin": "#9b59b6"})
    ax.set_ylabel("Transcriptomic Age (months)")
    ax.set_title("GSE280382 AgedST: Exenatide vs Rapamycin (13 weeks, whole blood)")
    ax.axhline(y=0, color="black", linestyle="--", alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "gse280382_AgedST_exen_vs_rapa.png", dpi=150)
    plt.close()
    print(f"  Saved {FIGURES_DIR / 'gse280382_AgedST_exen_vs_rapa.png'}")
    
    # Figure 3: Rejuvenation effect per tissue (barplot)
    fig, ax = plt.subplots(figsize=(14, 6))
    diff_df = pd.read_csv(RESULTS_DIR / "gse280382_AgedLT_differences.csv")
    diff_df = diff_df.sort_values("difference")
    colors = ["#2ecc71" if d < 0 else "#e74c3c" for d in diff_df["difference"]]
    ax.barh(range(len(diff_df)), diff_df["difference"], color=colors)
    ax.set_yticks(range(len(diff_df)))
    ax.set_yticklabels(diff_df["tissue"])
    ax.axvline(x=0, color="black", linestyle="-", alpha=0.5)
    ax.set_xlabel("Rejuvenation Effect (Exenatide - Control, months)")
    ax.set_title("GLP-1RA (Exenatide) Rejuvenation Effect by Tissue\nNegative = Rejuvenation (lower tAge)")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "gse280382_rejuvenation_by_tissue.png", dpi=150)
    plt.close()
    print(f"  Saved {FIGURES_DIR / 'gse280382_rejuvenation_by_tissue.png'}")


# ============================================
# Main
# ============================================
def main():
    print("=" * 60)
    print("GSE280382 Processing Pipeline")
    print("=" * 60)
    
    agedlt = process_agedlt()
    agedst = process_agedst()
    young = process_young()
    
    generate_comparisons(agedlt, agedst, young)
    make_figures(agedlt, agedst, young)
    
    print("\n" + "=" * 60)
    print("DONE. Results in:", RESULTS_DIR)
    print("Figures in:", FIGURES_DIR)
    print("=" * 60)


if __name__ == "__main__":
    main()
