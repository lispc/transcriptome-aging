#!/usr/bin/env python3
"""
Extract selected samples from ARCHS4 H5, run TMM normalization (edgeR),
align features, run tAge prediction, and save results.
"""
import os, sys, subprocess, tempfile, json
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
import joblib
import warnings
from sklearn.exceptions import InconsistentVersionWarning

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings("ignore", category=UserWarning)

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
H5_PATH = PROJECT_ROOT / "data" / "archs4" / "mouse_gene_v2.5.h5"
MODEL_PATH = PROJECT_ROOT / "models" / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def load_model_info():
    print("Loading tAge model...")
    model = joblib.load(MODEL_PATH)
    imputer = model.named_steps["imputation"]
    scaler = model.named_steps["scaler"]
    estimator = model.named_steps["estimator"]
    feature_names = list(model.feature_names_in_)
    coef = estimator.coef_
    intercept = float(estimator.intercept_)
    impute_stats = imputer.statistics_
    center_mean = scaler.mean_
    print(f"  Features: {len(feature_names)}")
    print(f"  Non-zero coefs: {(coef != 0).sum()}")
    return {
        "model": model,
        "feature_names": feature_names,
        "coef": coef,
        "intercept": intercept,
        "impute_stats": impute_stats,
        "center_mean": center_mean,
    }

def get_gene_symbols(h5_file):
    """Read gene symbols (Entrez IDs) from H5."""
    with h5py.File(h5_file, "r") as f:
        genes = f["meta"]["genes"][:]
    # Decode bytes
    genes = [g.decode("utf-8", errors="replace") if isinstance(g, bytes) else str(g) for g in genes]
    return genes

def extract_counts_for_samples(h5_file, sample_indices, out_tsv):
    """
    Extract count matrix for given sample indices from H5.
    Reads in gene chunks to keep memory low.
    """
    print(f"Extracting {len(sample_indices)} samples from {h5_file}...")
    with h5py.File(h5_file, "r") as f:
        data = f["data"]["expression"]
        n_genes, n_samples_h5 = data.shape
        print(f"  H5 shape: {n_genes} x {n_samples_h5}")

        # Read gene symbols
        genes = f["meta"]["genes"][:]
        genes = [g.decode("utf-8", errors="replace") if isinstance(g, bytes) else str(g) for g in genes]

        # Read sample metadata for column names
        sample_geo = f["meta"]["Sample_geo_accession"][:]
        sample_geo = [s.decode("utf-8", errors="replace") if isinstance(s, bytes) else str(s) for s in sample_geo]
        col_names = [sample_geo[i] for i in sample_indices]

        # Extract in chunks of genes to avoid memory spike
        chunk_size = 5000
        chunks = []
        for start in range(0, n_genes, chunk_size):
            end = min(start + chunk_size, n_genes)
            chunk = data[start:end, sample_indices]
            # h5py returns array; ensure it's dense
            chunk = np.array(chunk, dtype=np.float64)
            chunks.append(chunk)
            if (start // chunk_size) % 10 == 0:
                print(f"  Read genes {start}:{end}")

    mat = np.vstack(chunks)
    df = pd.DataFrame(mat, index=genes, columns=col_names)
    df.index.name = "gene_id"
    df.to_csv(out_tsv, sep="\t")
    print(f"  Saved counts to {out_tsv} ({df.shape})")
    return df

def run_edger_tmm(counts_file, output_prefix):
    logcpm_file = RESULTS_DIR / f"{output_prefix}_logcpm.tsv"
    r_script = f"""
suppressPackageStartupMessages(library(edgeR))
counts <- read.delim("{counts_file}", row.names=1, check.names=FALSE)
counts <- as.matrix(counts)
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge, method="TMM")
logcpm <- cpm(dge, log=TRUE, prior.count=3)
write.table(logcpm, file="{logcpm_file}", sep="\\t", quote=FALSE)
cat("logCPM saved to {logcpm_file}\\n")
cat("Dimensions:", nrow(logcpm), "x", ncol(logcpm), "\\n")
"""
    r_script_file = RESULTS_DIR / f"{output_prefix}_edger.R"
    with open(r_script_file, "w") as f:
        f.write(r_script)
    result = subprocess.run(["Rscript", str(r_script_file)], capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print("stderr:", result.stderr)
    if result.returncode != 0:
        raise RuntimeError("edgeR normalization failed")
    return logcpm_file

def tage_predict(logcpm_df, model_info):
    feature_names = model_info["feature_names"]
    impute_stats = model_info["impute_stats"]
    center_mean = model_info["center_mean"]
    coef = model_info["coef"]
    intercept = model_info["intercept"]

    logcpm_df.index = logcpm_df.index.astype(str)
    common_genes = list(logcpm_df.index.intersection(feature_names))
    n_samples = logcpm_df.shape[1]
    n_features = len(feature_names)
    X = np.zeros((n_samples, n_features), dtype=float)
    feat_to_idx = {name: i for i, name in enumerate(feature_names)}

    for gene in common_genes:
        X[:, feat_to_idx[gene]] = logcpm_df.loc[gene].values.astype(float)

    missing = 0
    for i, gene in enumerate(feature_names):
        if gene not in common_genes:
            X[:, i] = impute_stats[i]
            missing += 1
    print(f"  Common genes: {len(common_genes)}, Missing imputed: {missing}")

    X_centered = X - center_mean[np.newaxis, :]
    predictions = intercept + np.dot(X_centered, coef)
    return predictions

def process_sample_set(name, sample_indices, metadata_df):
    print(f"\n{'='*60}")
    print(f"Processing: {name} ({len(sample_indices)} samples)")
    print(f"{'='*60}")

    counts_file = RESULTS_DIR / f"{name}_counts.tsv"
    if not counts_file.exists():
        extract_counts_for_samples(H5_PATH, sample_indices, counts_file)
    else:
        print(f"  Using existing counts: {counts_file}")

    logcpm_file = RESULTS_DIR / f"{name}_logcpm.tsv"
    if not logcpm_file.exists():
        logcpm_file = run_edger_tmm(counts_file, name)
    else:
        print(f"  Using existing logCPM: {logcpm_file}")

    logcpm_df = pd.read_csv(logcpm_file, sep="\t", index_col=0)
    model_info = load_model_info()
    preds = tage_predict(logcpm_df, model_info)

    meta_subset = metadata_df.iloc[sample_indices].copy()
    meta_subset["tAge"] = preds
    # Apply mouse species adjustment (48 months)
    meta_subset["tAge_adj"] = preds * 48.0

    out = RESULTS_DIR / f"{name}_predictions.csv"
    meta_subset.to_csv(out, index=False)
    print(f"  Saved predictions to {out}")
    return meta_subset

def main():
    if not H5_PATH.exists():
        print(f"H5 file not found: {H5_PATH}")
        sys.exit(1)

    # Load metadata
    meta_file = RESULTS_DIR / "archs4_sample_metadata.csv.gz"
    if not meta_file.exists():
        print("Run explore_archs4_h5.py first to extract metadata.")
        sys.exit(1)

    print("Loading metadata...")
    metadata_df = pd.read_csv(meta_file, compression="gzip", low_memory=False)
    print(f"  Total samples: {len(metadata_df)}")

    # Define sample sets by keyword search on metadata
    text = (
        metadata_df.get("Sample_title", "") + " " +
        metadata_df.get("Sample_characteristics_ch1", "") + " " +
        metadata_df.get("Sample_source_name_ch1", "")
    )
    text = text.fillna("").str.lower()

    def get_indices(keywords, extra_mask=None):
        mask = pd.Series(False, index=metadata_df.index)
        for kw in keywords:
            mask |= text.str.contains(kw.lower(), na=False)
        if extra_mask is not None:
            mask &= extra_mask
        idx = metadata_df.index[mask].tolist()
        return idx

    # Limit sample sets to reasonable sizes to keep compute manageable
    max_n = 500

    sets = {
        "rapamycin": get_indices(["rapamycin"]),
        "cr_liver": get_indices(["caloric restriction", "calorie restriction"]),
        "metformin": get_indices(["metformin"]),
        "aging_liver": get_indices(["aging", "aged", "old", "elderly", "young"]),
        "dexamethasone": get_indices(["dexamethasone"]),
        "doxorubicin": get_indices(["doxorubicin"]),
    }

    all_results = []
    for name, idx in sets.items():
        if not idx:
            print(f"  {name}: no samples found")
            continue
        idx = idx[:max_n]
        print(f"  {name}: {len(idx)} samples (limited to {max_n})")
        res = process_sample_set(name, idx, metadata_df)
        all_results.append(res)

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(RESULTS_DIR / "archs4_all_predictions.csv", index=False)
        print(f"\nSaved combined predictions ({len(combined)} rows)")

if __name__ == "__main__":
    main()
