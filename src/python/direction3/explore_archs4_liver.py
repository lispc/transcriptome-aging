#!/usr/bin/env python3
"""
Explore ARCHS4 liver tAge predictions across studies.
35,422 liver samples from hundreds of GEO series.
"""

import sys
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
from collections import Counter, defaultdict

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"

def log(msg):
    print(msg, flush=True)

# Load predictions
log("Loading tAge predictions...")
preds = pd.read_csv(RESULTS_DIR / "archs4_liver_all_predictions.csv")
preds = preds[preds["tAge_adj"].notna()].copy()
log(f"Valid predictions: {len(preds)}")

# Load ARCHS4 metadata for liver samples
log("Loading ARCHS4 metadata...")
sample_indices = preds["sample_idx"].values

with h5py.File(PROJECT_ROOT / "data" / "archs4" / "mouse_gene_v2.5.h5", "r") as f:
    meta = f["meta"]["samples"]
    
    series_list = []
    title_list = []
    source_list = []
    char_list = []
    geo_list = []
    
    for idx in sample_indices:
        s = meta["series_id"][idx]
        series_list.append(s.decode("utf-8", errors="replace") if isinstance(s, bytes) else str(s) if s is not None else "")
        
        t = meta["title"][idx]
        title_list.append(t.decode("utf-8", errors="replace") if isinstance(t, bytes) else str(t) if t is not None else "")
        
        src = meta["source_name_ch1"][idx]
        source_list.append(src.decode("utf-8", errors="replace") if isinstance(src, bytes) else str(src) if src is not None else "")
        
        c = meta["characteristics_ch1"][idx]
        char_list.append(c.decode("utf-8", errors="replace") if isinstance(c, bytes) else str(c) if c is not None else "")
        
        g = meta["geo_accession"][idx]
        geo_list.append(g.decode("utf-8", errors="replace") if isinstance(g, bytes) else str(g) if g is not None else "")

meta_df = pd.DataFrame({
    "sample_idx": sample_indices,
    "geo_accession": geo_list,
    "series_id": series_list,
    "title": title_list,
    "source_name": source_list,
    "characteristics": char_list,
})

# Merge with predictions
merged = pd.merge(preds, meta_df, on="sample_idx", how="left")
log(f"Merged: {len(merged)} samples")

# Save full merged data
merged.to_csv(RESULTS_DIR / "archs4_liver_merged_metadata.csv", index=False)
log("Saved merged metadata")

# ============================================================
# Analysis 1: Per-series statistics
# ============================================================
log("\n=== Analysis 1: Per-series statistics ===")

# A sample may belong to multiple series (comma-separated)
series_stats = defaultdict(lambda: {"count": 0, "tAge_sum": 0.0, "tAge_sq_sum": 0.0, "samples": []})

for _, row in merged.iterrows():
    tAge = row["tAge_adj"]
    series_ids = [s.strip() for s in str(row["series_id"]).split(",") if s.strip()]
    for sid in series_ids:
        if sid.startswith("GSE") or sid.startswith("GDS"):
            series_stats[sid]["count"] += 1
            series_stats[sid]["tAge_sum"] += tAge
            series_stats[sid]["tAge_sq_sum"] += tAge ** 2
            series_stats[sid]["samples"].append(row["sample_idx"])

series_rows = []
for sid, stats in series_stats.items():
    n = stats["count"]
    if n >= 3:  # at least 3 samples
        mean_tAge = stats["tAge_sum"] / n
        var_tAge = stats["tAge_sq_sum"] / n - mean_tAge ** 2
        std_tAge = np.sqrt(max(0, var_tAge))
        series_rows.append({
            "series": sid,
            "n_samples": n,
            "mean_tAge": mean_tAge,
            "std_tAge": std_tAge,
            "min_tAge": min(stats["samples"]),
            "max_tAge": max(stats["samples"]),
        })

series_df = pd.DataFrame(series_rows)
series_df = series_df.sort_values("mean_tAge")
series_df.to_csv(RESULTS_DIR / "archs4_liver_per_series.csv", index=False)
log(f"Series with >=3 samples: {len(series_df)}")
log(f"Lowest mean tAge: {series_df.iloc[0]['series']} = {series_df.iloc[0]['mean_tAge']:.1f} (n={series_df.iloc[0]['n_samples']})")
log(f"Highest mean tAge: {series_df.iloc[-1]['series']} = {series_df.iloc[-1]['mean_tAge']:.1f} (n={series_df.iloc[-1]['n_samples']})")

# ============================================================
# Analysis 2: Top/Bottom series with sample titles
# ============================================================
log("\n=== Analysis 2: Youngest series (lowest tAge) ===")
for _, row in series_df.head(20).iterrows():
    sid = row["series"]
    samples = merged[merged["series_id"].str.contains(sid, na=False)]
    titles = samples["title"].tolist()
    log(f"{sid}: mean={row['mean_tAge']:.1f}, n={row['n_samples']}, titles={titles[:3]}")

log("\n=== Analysis 3: Oldest series (highest tAge) ===")
for _, row in series_df.tail(20).iterrows():
    sid = row["series"]
    samples = merged[merged["series_id"].str.contains(sid, na=False)]
    titles = samples["title"].tolist()
    log(f"{sid}: mean={row['mean_tAge']:.1f}, n={row['n_samples']}, titles={titles[:3]}")

# ============================================================
# Analysis 3: Keyword-based intervention classification
# ============================================================
log("\n=== Analysis 4: Keyword-based intervention classification ===")

def classify_sample(title, chars, source):
    text = f"{title} {chars} {source}".lower()
    
    # Diet interventions
    if any(k in text for k in ["caloric restriction", "calorie restriction", "cr diet", "restricted diet"]):
        return "Diet_CR"
    if any(k in text for k in ["ketogenic", "keto diet", "high fat"]):
        return "Diet_Keto"
    if any(k in text for k in ["fasting", "starvation", "food deprivation"]):
        return "Diet_Fasting"
    if any(k in text for k in ["high-fat", "high fat diet", "hfd", "western diet", "obesogenic"]):
        return "Diet_HighFat"
    if any(k in text for k in ["methionine restriction", "mr diet", "methionine deficient"]):
        return "Diet_MethionineRestriction"
    
    # Drugs
    if any(k in text for k in ["rapamycin", "sirolimus", "mTOR"]):
        return "Drug_Rapamycin"
    if any(k in text for k in ["metformin", "biguanide"]):
        return "Drug_Metformin"
    if any(k in text for k in ["resveratrol"]):
        return "Drug_Resveratrol"
    if any(k in text for k in ["trametinib", "mek inhibitor", "mek1/2"]):
        return "Drug_Trametinib"
    if any(k in text for k in ["dexamethasone", "corticosterone", "glucocorticoid", "prednisolone", "hydrocortisone"]):
        return "Drug_Glucocorticoid"
    if any(k in text for k in ["doxorubicin", "cisplatin", "bleomycin", "chemotherapy"]):
        return "Drug_Chemotherapy"
    if any(k in text for k in ["acetaminophen", "apap", "paracetamol", "toxic"]):
        return "Drug_Toxicant"
    if any(k in text for k in ["nmn", "nicotinamide", "nad+", "nr supplement"]):
        return "Drug_NAD"
    if any(k in text for k in ["spermidine", "spermine"]):
        return "Drug_Spermidine"
    
    # Genetic
    if any(k in text for k in ["knockout", "ko mouse", "-/-", "null", "deficient", "deletion", "cKO", "flox"]):
        return "Genetic_KO"
    if any(k in text for k in ["overexpression", "transgenic", "tg mouse", "oe-"]):
        return "Genetic_OE"
    
    # Disease
    if any(k in text for k in ["nash", "steatosis", "steatohepatitis", "fatty liver", "mash"]):
        return "Disease_NASH"
    if any(k in text for k in ["hcc", "hepatocellular carcinoma", "liver cancer", "tumor", "carcinoma", "oncogene"]):
        return "Disease_LiverCancer"
    if any(k in text for k in ["fibrosis", "ccl4", "carbon tetrachloride", "thioacetamide", "taa"]):
        return "Disease_Fibrosis"
    if any(k in text for k in ["cirrhosis", "alcoholic", "ald", "alcohol"]):
        return "Disease_ALD"
    if any(k in text for k in ["diabetes", "stz", "streptozotocin", "insulin resistance", "db/db", "ob/ob"]):
        return "Disease_Diabetes"
    
    # Age
    if any(k in text for k in ["old", "aged", "aging", "elderly", "senescent", "senescence", "24 month", "28 month", "30 month", "18 month"]):
        return "Age_Old"
    if any(k in text for k in ["young", "adult", "2 month", "3 month", "4 month", "6 month", "8 week", "12 week"]):
        return "Age_Young"
    
    return "Other"

merged["intervention"] = merged.apply(lambda r: classify_sample(r["title"], r["characteristics"], r["source_name"]), axis=1)

intervention_stats = merged.groupby("intervention").agg(
    n_samples=("tAge_adj", "count"),
    mean_tAge=("tAge_adj", "mean"),
    std_tAge=("tAge_adj", "std"),
    median_tAge=("tAge_adj", "median"),
    min_tAge=("tAge_adj", "min"),
    max_tAge=("tAge_adj", "max"),
).reset_index().sort_values("mean_tAge")

intervention_stats.to_csv(RESULTS_DIR / "archs4_liver_intervention_stats.csv", index=False)
log(f"\nIntervention groups (n>=5):")
for _, row in intervention_stats[intervention_stats["n_samples"] >= 5].iterrows():
    log(f"  {row['intervention']}: n={row['n_samples']}, mean={row['mean_tAge']:.1f}±{row['std_tAge']:.1f}, median={row['median_tAge']:.1f}")

# ============================================================
# Analysis 5: Compare known interventions
# ============================================================
log("\n=== Analysis 5: Known interventions in ARCHS4 ===")

known_series = {
    "GSE288795": "Trametinib+Rapamycin combo",
    "GSE280382": "GLP-1RA Exenatide",
    "GSE230402": "Caloric Restriction",
    "GSE11923": "Rapamycin",
    "GSE55140": "Rapamycin",
}

for gse, name in known_series.items():
    sub = series_df[series_df["series"] == gse]
    if len(sub) > 0:
        row = sub.iloc[0]
        log(f"  {gse} ({name}): mean={row['mean_tAge']:.1f}, n={row['n_samples']}")
    else:
        log(f"  {gse} ({name}): NOT FOUND in liver samples")

# ============================================================
# Analysis 6: Individual outlier samples
# ============================================================
log("\n=== Analysis 6: Youngest individual samples ===")
youngest = merged.nsmallest(20, "tAge_adj")
for _, row in youngest.iterrows():
    log(f"  tAge={row['tAge_adj']:.1f}: {row['geo_accession']} | {row['title'][:60]} | {row['series_id'][:40]}")

log("\n=== Analysis 7: Oldest individual samples ===")
oldest = merged.nlargest(20, "tAge_adj")
for _, row in oldest.iterrows():
    log(f"  tAge={row['tAge_adj']:.1f}: {row['geo_accession']} | {row['title'][:60]} | {row['series_id'][:40]}")

log("\nDone! Results saved to results/direction3/")
