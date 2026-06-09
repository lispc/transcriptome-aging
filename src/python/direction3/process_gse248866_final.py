#!/usr/bin/env python3
"""
Process GSE248866 ZT12 subset: merge Kallisto abundances -> TMM -> tAge prediction.
"""

import sys
import pickle
import subprocess
import tempfile
import os
import numpy as np
import pandas as pd
import joblib
from pathlib import Path
from collections import defaultdict

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
KALLISTO_DIR = PROJECT_ROOT / "data" / "gse248866_kallisto"

def log(msg):
    print(msg, flush=True)

# Load model
log("Loading tAge model...")
model = joblib.load(PROJECT_ROOT / "models" / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
imputer = model.named_steps["imputation"]
scaler = model.named_steps["scaler"]
estimator = model.named_steps["estimator"]
coef = estimator.coef_
intercept = float(estimator.intercept_)
model_features = list(model.feature_names_in_)
impute_stats = imputer.statistics_
center_mean = scaler.mean_

# Load Ensembl->Entrez mapping
log("Loading Ensembl->Entrez mapping...")
with open(RESULTS_DIR / "gse248866_ensembl_to_entrez.pkl", "rb") as f:
    cache = pickle.load(f)
ensembl_to_entrez = cache["entrez"]

# Sample mapping
samples = {
    "AL_1": "SRR26975750",
    "AL_2": "SRR26975749",
    "AL_3": "SRR26975748",
    "CR_1": "SRR26975747",
    "CR_2": "SRR26975746",
    "CR_3": "SRR26975745",
    "NR_1": "SRR33542637",
    "NR_2": "SRR33542636",
    "NR_3": "SRR33542635",
    "Cort_1": "SRR33542634",
    "Cort_2": "SRR33542633",
    "Cort_3": "SRR33542632",
}

# ============================================================
# Step 1: Merge Kallisto est_counts (for TMM)
# ============================================================
log("\n=== Step 1: Merging Kallisto est_counts ===")

counts_matrix = None
tx_ids = None

for name in samples:
    tsv = KALLISTO_DIR / name / "abundance.tsv"
    df = pd.read_csv(tsv, sep="\t", usecols=["target_id", "est_counts"])
    
    if tx_ids is None:
        tx_ids = df["target_id"].values
        counts_matrix = pd.DataFrame(index=tx_ids)
    
    counts_matrix[name] = df["est_counts"].values
    log(f"  {name}: {len(df)} transcripts")

log(f"Merged counts matrix: {counts_matrix.shape}")

# ============================================================
# Step 2: Transcript -> Gene aggregation (sum counts)
# ============================================================
log("\n=== Step 2: Transcript -> Gene aggregation ===")

# Load transcript->gene mapping
tx2gene = pd.read_csv(
    PROJECT_ROOT / "data" / "ref" / "mus_musculus" / "transcripts_to_genes.txt",
    sep="\t", header=None, names=["transcript_id", "gene_id", "gene_name"]
)

tx_to_gene = dict(zip(tx2gene["transcript_id"], tx2gene["gene_id"]))

# Map transcripts to genes (sum counts)
gene_counts = defaultdict(lambda: pd.Series(0.0, index=counts_matrix.columns))
for tx in counts_matrix.index:
    gene = tx_to_gene.get(tx)
    if gene:
        gene_counts[gene] += counts_matrix.loc[tx]

gene_df = pd.DataFrame(gene_counts).T
log(f"Gene-level counts: {gene_df.shape}")

# ============================================================
# Step 3: Ensembl gene ID -> Entrez mapping
# ============================================================
log("\n=== Step 3: Ensembl -> Entrez mapping ===")

gene_df["ensembl_clean"] = gene_df.index.str.split("\.").str[0]
gene_df["entrez"] = gene_df["ensembl_clean"].map(ensembl_to_entrez)

mapped = gene_df["entrez"].notna().sum()
log(f"Mapped {mapped}/{len(gene_df)} genes to Entrez")

# Aggregate by Entrez (sum counts for duplicate Entrez IDs)
entrez_counts = gene_df.groupby("entrez")[list(samples.keys())].sum()
entrez_counts = entrez_counts.loc[entrez_counts.index.notna()]
entrez_counts.index = entrez_counts.index.astype(str)
log(f"Entrez-level counts: {entrez_counts.shape}")

# ============================================================
# Step 4: TMM normalization via Rscript
# ============================================================
log("\n=== Step 4: TMM normalization ===")

# Align to model features
expr_aligned = pd.DataFrame(index=model_features, columns=entrez_counts.columns)
common = entrez_counts.index.intersection(model_features)
log(f"{len(common)} features overlap with model")
expr_aligned.loc[common] = entrez_counts.loc[common].astype(float)
expr_aligned = expr_aligned.fillna(0.0).astype(float)

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

# ============================================================
# Step 5: tAge prediction
# ============================================================
log("\n=== Step 5: tAge prediction ===")

X_arr = logcpm_df.values.T
X_imp = np.where(np.isnan(X_arr), impute_stats, X_arr)
X_centered = X_imp - center_mean
y = intercept + X_centered @ coef
y_adj = (y - 5.5) * 48

# Map sample names to conditions
condition_map = {
    "AL_1": "AL", "AL_2": "AL", "AL_3": "AL",
    "CR_1": "CR", "CR_2": "CR", "CR_3": "CR",
    "NR_1": "NR", "NR_2": "NR", "NR_3": "NR",
    "Cort_1": "Cort", "Cort_2": "Cort", "Cort_3": "Cort",
}

preds = pd.DataFrame({
    "sample": list(samples.keys()),
    "srr": list(samples.values()),
    "condition": [condition_map[s] for s in samples.keys()],
    "tAge": y,
    "tAge_adj": y_adj,
})

preds.to_csv(RESULTS_DIR / "gse248866_predictions.csv", index=False)
log(f"\nPredictions saved to gse248866_predictions.csv")

# Summary statistics
summary = preds.groupby("condition").agg(
    n=("tAge_adj", "count"),
    mean_tAge=("tAge_adj", "mean"),
    std_tAge=("tAge_adj", "std"),
    median_tAge=("tAge_adj", "median"),
    min_tAge=("tAge_adj", "min"),
    max_tAge=("tAge_adj", "max"),
).reset_index()

log("\n=== Summary ===")
print(summary.to_string(index=False))

# Pairwise comparisons (relative to AL)
log("\n=== Pairwise comparisons vs AL ===")
al_mean = preds[preds["condition"] == "AL"]["tAge_adj"].mean()
al_values = preds[preds["condition"] == "AL"]["tAge_adj"].values
for cond in ["CR", "NR", "Cort"]:
    sub = preds[preds["condition"] == cond]["tAge_adj"]
    diff = sub.mean() - al_mean
    # Simple t-test
    from scipy import stats
    t, p = stats.ttest_ind(al_values, sub.values)
    log(f"  {cond} vs AL: {diff:+.1f} months, t={t:.2f}, p={p:.4f}")
    log(f"    {cond}: {sub.mean():.1f}±{sub.std():.1f}, AL: {al_mean:.1f}")

log("\nDone!")
