#!/usr/bin/env python3
"""Process only remaining ARCHS4 liver chunks (15, 16, 17). 
Skips known-bad H5 columns (corrupted filter data)."""

import sys
import pickle
import subprocess
import tempfile
import os
import gc
import time
import numpy as np
import pandas as pd
import h5py
import joblib
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
MODEL_PATH = PROJECT_ROOT / "models" / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
H5_PATH = PROJECT_ROOT / "data" / "archs4" / "mouse_gene_v2.5.h5"

CHUNK_SIZE = 2000
N_GENES = 53511
# Known corrupted H5 columns (filter failure during read)
BAD_H5_COLS = {930747, 930749, 930750, 930753}

def log(msg):
    print(msg, flush=True)

# Load model
log("Loading model...")
model = joblib.load(MODEL_PATH)
imputer = model.named_steps["imputation"]
scaler = model.named_steps["scaler"]
estimator = model.named_steps["estimator"]
coef = estimator.coef_
intercept = float(estimator.intercept_)
model_features = list(model.feature_names_in_)
impute_stats = imputer.statistics_
center_mean = scaler.mean_

# Load liver indices
log("Loading liver indices...")
with open(RESULTS_DIR / "archs4_liver_indices.pkl", "rb") as f:
    liver_indices = pickle.load(f)
log(f"Total liver samples: {len(liver_indices)}")

# Load gene symbol mapping
log("Loading cached symbol→Entrez mapping...")
with open(RESULTS_DIR / "gse280382_symbol_to_entrez.pkl", "rb") as f:
    symbol_to_entrez = pickle.load(f)

# Get ARCHS4 gene symbols
log("Reading ARCHS4 gene symbols...")
with h5py.File(H5_PATH, "r") as f:
    archs4_symbols = [s.decode("utf-8", errors="replace") for s in f["meta"]["genes"]["symbol"][:]]

# Build mapping
entrez_ids = [symbol_to_entrez.get(sym) for sym in archs4_symbols]
mapped_count = sum(1 for e in entrez_ids if e is not None)
log(f"Mapped {mapped_count}/{len(archs4_symbols)} symbols")

from collections import defaultdict
entrez_to_rows = defaultdict(list)
for i, eid in enumerate(entrez_ids):
    if eid is not None:
        entrez_to_rows[eid].append(i)
log(f"Unique Entrez IDs: {len(entrez_to_rows)}")

n_chunks = (len(liver_indices) + CHUNK_SIZE - 1) // CHUNK_SIZE
remaining = [i for i in range(n_chunks) if not (RESULTS_DIR / f"archs4_liver_chunk_{i:03d}.csv").exists()]
log(f"Remaining chunks: {remaining}")

row_lists = list(entrez_to_rows.values())
entrez_list = list(entrez_to_rows.keys())

for chunk_i in remaining:
    start_idx = chunk_i * CHUNK_SIZE
    end_idx = min(start_idx + CHUNK_SIZE, len(liver_indices))
    chunk_sample_indices = liver_indices[start_idx:end_idx]
    n_samples = len(chunk_sample_indices)
    
    log(f"\n=== Chunk {chunk_i+1}/{n_chunks} ({n_samples} samples) ===")
    
    # Separate good and bad samples
    good_indices = []
    bad_indices = []
    for i, idx in enumerate(chunk_sample_indices):
        if idx in BAD_H5_COLS:
            bad_indices.append((i, idx))
        else:
            good_indices.append((i, idx))
    
    log(f"  Good samples: {len(good_indices)}, Bad samples: {len(bad_indices)}")
    
    # Read only good samples
    log("  Extracting counts from H5 (good samples only)...")
    expr_good = np.zeros((N_GENES, len(good_indices)), dtype=np.float64)
    with h5py.File(H5_PATH, "r") as f:
        ds = f["data"]["expression"]
        for i, (_, sample_idx) in enumerate(good_indices):
            if i % 200 == 0:
                log(f"    Sample {i}/{len(good_indices)} (h5_col={sample_idx})")
            for attempt in range(3):
                try:
                    expr_good[:, i] = ds[:, sample_idx]
                    break
                except OSError as e:
                    log(f"    H5 read failed (attempt {attempt+1}/3): col={sample_idx}, err={str(e)[:80]}")
                    time.sleep(1)
                    if attempt == 2:
                        raise
    
    log("  Aggregating by Entrez...")
    expr_agg = np.zeros((len(entrez_to_rows), len(good_indices)), dtype=np.float64)
    for i, rows in enumerate(row_lists):
        expr_agg[i, :] = expr_good[rows, :].mean(axis=0)
    
    good_sample_ids = [idx for _, idx in good_indices]
    df = pd.DataFrame(expr_agg, index=entrez_list, columns=[f"S{idx}" for idx in good_sample_ids])
    
    log("  Aligning to model features...")
    expr_aligned = pd.DataFrame(index=model_features, columns=df.columns)
    common = df.index.intersection(model_features)
    log(f"  {len(common)} features overlap")
    expr_aligned.loc[common] = df.loc[common].astype(float)
    expr_aligned = expr_aligned.fillna(0.0).astype(float)
    
    log("  TMM normalization via Rscript...")
    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = os.path.join(tmpdir, "counts.csv")
        output_csv = os.path.join(tmpdir, "logcpm.csv")
        expr_aligned.to_csv(input_csv)
        
        r_script = str(PROJECT_ROOT / "src" / "r" / "tmm_normalize.R")
        result = subprocess.run(
            ["Rscript", r_script, input_csv, output_csv],
            capture_output=True, text=True, check=True
        )
        log(f"  Rscript: {result.stdout.strip()}")
        
        logcpm_df = pd.read_csv(output_csv, index_col=0)
        logcpm_df.index = logcpm_df.index.astype(str)
        logcpm_df.columns = logcpm_df.columns.astype(str)
        logcpm_df = logcpm_df.reindex(index=model_features, columns=expr_aligned.columns)
    
    log("  Predicting tAge...")
    X_arr = logcpm_df.values.T
    X_imp = np.where(np.isnan(X_arr), impute_stats, X_arr)
    X_centered = X_imp - center_mean
    y_good = intercept + X_centered @ coef
    y_adj_good = (y_good - 5.5) * 48
    
    # Build full result with NaN for bad samples
    y = np.full(n_samples, np.nan)
    y_adj = np.full(n_samples, np.nan)
    for i, (_, _) in enumerate(good_indices):
        y[i] = y_good[i]
        y_adj[i] = y_adj_good[i]
    
    chunk_preds = pd.DataFrame({
        "sample_idx": chunk_sample_indices,
        "tAge": y,
        "tAge_adj": y_adj,
    })
    chunk_file = RESULTS_DIR / f"archs4_liver_chunk_{chunk_i:03d}.csv"
    chunk_preds.to_csv(chunk_file, index=False)
    valid = chunk_preds["tAge_adj"].notna()
    log(f"  Saved chunk {chunk_i} ({valid.sum()}/{n_samples} valid, range: {chunk_preds.loc[valid,'tAge_adj'].min():.1f} to {chunk_preds.loc[valid,'tAge_adj'].max():.1f})")
    
    del expr_good, expr_agg, df, expr_aligned, logcpm_df, X_arr, X_imp, X_centered
    gc.collect()

log("\nDone! All remaining chunks processed.")
