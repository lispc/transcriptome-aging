#!/usr/bin/env python3
"""Quick analysis of GSE299228 using group-averaged DESeq2 normalized counts."""

import pandas as pd
import numpy as np
import joblib
import sys
from pathlib import Path

DATA_DIR = Path('/home/scroll/zzhang/transcriptome-aging/data/geo')
RESULTS_DIR = Path('/home/scroll/zzhang/transcriptome-aging/results/geo_longterm')
MODEL_PATH = Path('/home/scroll/zzhang/transcriptome-aging/models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl')

# Load model
print("Loading model...")
model = joblib.load(MODEL_PATH)
imputer = model.named_steps["imputation"]
scaler = model.named_steps["scaler"]
estimator = model.named_steps["estimator"]
feature_names = list(model.feature_names_in_)
coef = estimator.coef_
intercept = float(estimator.intercept_)
impute_stats = imputer.statistics_
center_mean = scaler.mean_

# Load GSE299228 comparative study (group averages)
print("\nLoading GSE299228 group averages...")
df = pd.read_csv(DATA_DIR / 'GSE299228_Comparative_Study_DESeq2.csv.gz', compression='gzip')
print(f"  Shape: {df.shape}")
print(f"  Columns: {df.columns.tolist()[:7]}")

# Gene IDs are in first column (Gene Symbol with Ensembl version)
genes = df.iloc[:, 0].values
print(f"  First 5 gene IDs: {genes[:5]}")

# Extract group average counts
group_cols = ['Cntrl', 'Rapamycin', 'Metformin', 'TM5614', 'CalRestrn']
counts = df[group_cols].values

# Map Gene Symbols to Entrez IDs using mygene
print("\nMapping Gene Symbols to Entrez IDs...")
from mygene import MyGeneInfo
mg = MyGeneInfo()

# Clean gene IDs (remove Ensembl version)
clean_genes = [str(g).split('.')[0] for g in genes]

mapping = {}
batch_size = 1000
for i in range(0, len(clean_genes), batch_size):
    batch = clean_genes[i:i+batch_size]
    result = mg.querymany(batch, scopes='symbol,ensembl.gene', 
                         fields='entrezgene', species='mouse',
                         verbose=False, as_dataframe=True)
    for idx, row in result.iterrows():
        if 'entrezgene' in row and pd.notna(row['entrezgene']):
            mapping[idx] = int(row['entrezgene'])

print(f"  Mapped: {len(mapping)}/{len(clean_genes)}")

# Build expression matrix with Entrez IDs
entrez_list = []
expr_list = []
for i, gene in enumerate(clean_genes):
    if gene in mapping:
        entrez_list.append(str(mapping[gene]))
        expr_list.append(counts[i])

expr_df = pd.DataFrame(np.array(expr_list), columns=group_cols, index=entrez_list)
expr_df = expr_df.groupby(level=0).mean()  # aggregate duplicates
print(f"  After aggregation: {len(expr_df)} unique Entrez IDs")

# Note: These are DESeq2 size-factor normalized counts (not raw counts)
# For consistency with our pipeline, we'll log2-transform them
# DESeq2 normalized counts are approximately on raw count scale
log_expr = np.log2(expr_df + 1)  # log2(x+1) pseudo-count

# Align with model features
n_features = len(feature_names)
X = np.zeros((len(group_cols), n_features), dtype=float)
feat_to_idx = {name: i for i, name in enumerate(feature_names)}

common = 0
for gene in log_expr.index:
    if gene in feat_to_idx:
        X[:, feat_to_idx[gene]] = log_expr.loc[gene].values
        common += 1

print(f"\nCommon genes with model: {common}/{n_features}")

# Impute missing
missing = n_features - common
for i in range(n_features):
    if np.all(X[:, i] == 0):
        X[:, i] = impute_stats[i]

# Center and predict
X_centered = X - center_mean[np.newaxis, :]
predictions = intercept + np.dot(X_centered, coef)

print("\n" + "=" * 60)
print("GSE299228 tAge Predictions (from group-averaged counts)")
print("=" * 60)
for i, group in enumerate(group_cols):
    print(f"  {group:12s}: {predictions[i]:.4f}")

# Compute differences from Control
ctrl_pred = predictions[0]  # Cntrl is first
print(f"\n  Control baseline: {ctrl_pred:.4f}")
print(f"\n  Drug - Control (negative = rejuvenation):")
for i, group in enumerate(group_cols[1:], 1):
    diff = predictions[i] - ctrl_pred
    print(f"  {group:12s}: {diff:+.4f}")

# Save
result_df = pd.DataFrame({
    'group': group_cols,
    'tage_prediction': predictions,
    'difference_from_control': predictions - ctrl_pred
})
result_df.to_csv(RESULTS_DIR / 'gse299228_predictions.csv', index=False)
print(f"\nSaved to {RESULTS_DIR / 'gse299228_predictions.csv'}")
