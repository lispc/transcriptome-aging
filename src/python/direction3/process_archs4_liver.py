#!/usr/bin/env python3
"""
Process ARCHS4 liver samples in chunks for tAge prediction.
35,422 liver samples × 53,511 genes.
Uses small-batch H5 reads to avoid filter failures on long runs.
"""

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
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHUNK_SIZE = 2000
H5_BATCH = 100  # read H5 in small batches to avoid filter errors

def log(msg):
    print(msg, flush=True)

def read_h5_small_batches(h5_path, sample_indices, n_genes=53511):
    """Read H5 expression data in small batches to avoid filter failures."""
    n_samples = len(sample_indices)
    expr = np.zeros((n_genes, n_samples), dtype=np.float64)
    with h5py.File(h5_path, "r") as f:
        ds = f["data"]["expression"]
        for start in range(0, n_samples, H5_BATCH):
            end = min(start + H5_BATCH, n_samples)
            batch_idx = sample_indices[start:end]
            for attempt in range(3):
                try:
                    expr[:, start:end] = ds[:, batch_idx]
                    break
                except OSError as e:
                    log(f"    H5 batch read failed (attempt {attempt+1}/3): {e}")
                    time.sleep(1)
                    if attempt == 2:
                        raise
    return expr

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
cache = RESULTS_DIR / "gse280382_symbol_to_entrez.pkl"
if not cache.exists():
    log("ERROR: No symbol mapping cache found!")
    sys.exit(1)
log("Loading cached symbol→Entrez mapping...")
with open(cache, "rb") as f:
    symbol_to_entrez = pickle.load(f)

# Get ARCHS4 gene symbols
log("Reading ARCHS4 gene symbols...")
with h5py.File(PROJECT_ROOT / "data" / "archs4" / "mouse_gene_v2.5.h5", "r") as f:
    archs4_symbols = [s.decode("utf-8", errors="replace") for s in f["meta"]["genes"]["symbol"][:]]

# Build mapping
log(f"Mapping {len(archs4_symbols)} ARCHS4 symbols to Entrez...")
entrez_ids = []
for sym in archs4_symbols:
    entrez_ids.append(symbol_to_entrez.get(sym))

mapped_count = sum(1 for e in entrez_ids if e is not None)
log(f"  Mapped {mapped_count}/{len(archs4_symbols)} symbols")

from collections import defaultdict
entrez_to_rows = defaultdict(list)
for i, eid in enumerate(entrez_ids):
    if eid is not None:
        entrez_to_rows[eid].append(i)

log(f"  Unique Entrez IDs: {len(entrez_to_rows)}")

# Check for resume
n_chunks = (len(liver_indices) + CHUNK_SIZE - 1) // CHUNK_SIZE
completed_chunks = set()
for i in range(n_chunks):
    chunk_file = RESULTS_DIR / f"archs4_liver_chunk_{i:03d}.csv"
    if chunk_file.exists():
        completed_chunks.add(i)
        log(f"  Chunk {i} already completed (resuming)")

log(f"Completed chunks: {len(completed_chunks)}/{n_chunks}")

# Load completed predictions
all_preds = []
for i in sorted(completed_chunks):
    chunk_file = RESULTS_DIR / f"archs4_liver_chunk_{i:03d}.csv"
    all_preds.append(pd.read_csv(chunk_file))

H5_PATH = PROJECT_ROOT / "data" / "archs4" / "mouse_gene_v2.5.h5"

# Process remaining chunks
for chunk_i in range(n_chunks):
    if chunk_i in completed_chunks:
        continue
    
    start_idx = chunk_i * CHUNK_SIZE
    end_idx = min(start_idx + CHUNK_SIZE, len(liver_indices))
    chunk_sample_indices = liver_indices[start_idx:end_idx]
    n_samples = len(chunk_sample_indices)
    
    log(f"\n=== Chunk {chunk_i+1}/{n_chunks} ({n_samples} samples) ===")
    
    log("  Extracting counts from H5 (small batches)...")
    expr = read_h5_small_batches(H5_PATH, chunk_sample_indices)
    log(f"  Extracted shape: {expr.shape}")
    
    log("  Aggregating by Entrez...")
    expr_agg = np.zeros((len(entrez_to_rows), n_samples), dtype=np.float64)
    agg_index = []
    for i, (eid, rows) in enumerate(entrez_to_rows.items()):
        expr_agg[i, :] = expr[rows, :].mean(axis=0)
        agg_index.append(eid)
    
    df = pd.DataFrame(expr_agg, index=agg_index, columns=[f"S{idx}" for idx in chunk_sample_indices])
    
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
    y = intercept + X_centered @ coef
    y_adj = (y - 5.5) * 48
    
    chunk_preds = pd.DataFrame({
        "sample_idx": chunk_sample_indices,
        "tAge": y,
        "tAge_adj": y_adj,
    })
    chunk_preds.to_csv(RESULTS_DIR / f"archs4_liver_chunk_{chunk_i:03d}.csv", index=False)
    log(f"  Saved chunk {chunk_i} ({y_adj.min():.1f} to {y_adj.max():.1f})")
    
    all_preds.append(chunk_preds)
    
    del expr, expr_agg, df, expr_aligned, logcpm_df, X_arr, X_imp, X_centered
    gc.collect()

# Combine
log("\n=== Combining all predictions ===")
final_preds = pd.concat(all_preds, ignore_index=True)
final_preds.to_csv(RESULTS_DIR / "archs4_liver_all_predictions.csv", index=False)
log(f"Saved {len(final_preds)} predictions")
log(f"tAge_adj range: {final_preds['tAge_adj'].min():.1f} to {final_preds['tAge_adj'].max():.1f}")
log(f"Mean tAge_adj: {final_preds['tAge_adj'].mean():.1f}")
