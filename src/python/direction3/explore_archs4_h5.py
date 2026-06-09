#!/usr/bin/env python3
"""
Explore ARCHS4 mouse_gene_v2.5.h5 structure and extract sample metadata.
Run after the H5 file is fully downloaded.
"""
import os, sys, json, time
from pathlib import Path
import h5py
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
H5_PATH = PROJECT_ROOT / "data" / "archs4" / "mouse_gene_v2.5.h5"
OUT_DIR = PROJECT_ROOT / "results" / "direction3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def explore_structure():
    print(f"Opening {H5_PATH}")
    with h5py.File(H5_PATH, "r") as f:
        print("Top-level keys:", list(f.keys()))
        if "meta" in f:
            print("Meta keys:", list(f["meta"].keys()))
        if "data" in f:
            print("Data keys:", list(f["data"].keys()))
            for k in f["data"].keys():
                d = f[f"data/{k}"]
                print(f"  data/{k}: shape={d.shape}, dtype={d.dtype}, chunks={d.chunks}")

def extract_metadata():
    print("\nExtracting sample metadata...")
    with h5py.File(H5_PATH, "r") as f:
        meta = f["meta"]
        n_samples = meta["Sample_geo_accession"].shape[0]
        print(f"Total samples: {n_samples}")

        # Decode byte strings safely
        def decode(idx, key):
            try:
                val = meta[key][idx]
                if isinstance(val, bytes):
                    return val.decode("utf-8", errors="replace")
                return str(val)
            except Exception:
                return ""

        records = []
        # Batch extraction to avoid too many individual reads
        keys = [
            "Sample_geo_accession",
            "Sample_title",
            "Series_geo_accession",
            "Sample_source_name_ch1",
            "Sample_characteristics_ch1",
            "Sample_library_strategy",
            "Sample_taxid_ch1",
            "Sample_organism_ch1",
        ]
        # Check which keys exist
        available_keys = [k for k in keys if k in meta]
        print(f"Available keys: {available_keys}")

        batch_size = 10000
        for start in range(0, n_samples, batch_size):
            end = min(start + batch_size, n_samples)
            print(f"  Reading batch {start}:{end}")
            batch = {k: meta[k][start:end] for k in available_keys}
            for i in range(end - start):
                rec = {}
                for k in available_keys:
                    v = batch[k][i]
                    if isinstance(v, bytes):
                        v = v.decode("utf-8", errors="replace")
                    rec[k] = str(v)
                records.append(rec)

    df = pd.DataFrame(records)
    out = OUT_DIR / "archs4_sample_metadata.csv.gz"
    df.to_csv(out, index=False, compression="gzip")
    print(f"Saved metadata to {out} ({len(df)} rows)")
    return df

def search_samples(metadata_df, keywords):
    """Search metadata for keywords in title and characteristics."""
    print(f"\nSearching for {keywords}...")
    mask = pd.Series(False, index=metadata_df.index)
    text_cols = [c for c in metadata_df.columns if "title" in c.lower() or "characteristics" in c.lower() or "source" in c.lower()]
    for col in text_cols:
        for kw in keywords:
            mask |= metadata_df[col].str.contains(kw, case=False, na=False)
    hits = metadata_df[mask]
    print(f"  Found {len(hits)} samples")
    return hits

def main():
    if not H5_PATH.exists():
        print(f"H5 file not found: {H5_PATH}")
        sys.exit(1)
    explore_structure()
    df = extract_metadata()

    # Quick keyword searches
    searches = {
        "rapamycin": ["rapamycin"],
        "caloric_restriction": ["caloric restriction", "calorie restriction", "cr diet"],
        "metformin": ["metformin"],
        "aging_liver": ["aging", "aged", "old", "elderly"],
        "dexamethasone": ["dexamethasone"],
        "doxorubicin_heart": ["doxorubicin"],
        "young_old_liver": ["young", "old", "liver"],
    }
    for name, kws in searches.items():
        hits = search_samples(df, kws)
        hits.to_csv(OUT_DIR / f"archs4_search_{name}.csv", index=False)

if __name__ == "__main__":
    main()
