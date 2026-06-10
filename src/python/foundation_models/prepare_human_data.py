"""
Preprocess human aging scRNA-seq data (AIDA + Kedlian) for foundation model training.
Uses backed mode for large files to avoid memory issues.
"""
import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

DATA_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data/scrna_aging/human")
OUTPUT_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def process_aida():
    """Process AIDA PBMC data."""
    print("=" * 60)
    print("Processing AIDA v1")
    print("=" * 60)
    
    # Use backed mode to avoid loading full matrix
    adata = sc.read_h5ad(DATA_DIR / "aida_v1.h5ad", backed='r')
    print(f"Loaded (backed): {adata.shape}")
    
    # Extract age from development_stage
    # Format: "33-year-old stage"
    def parse_age(stage_str):
        s = str(stage_str).lower()
        match = re.search(r'(\d+)-year-old', s)
        if match:
            return float(match.group(1))
        return np.nan
    
    ages = adata.obs['development_stage'].apply(parse_age)
    ages = pd.to_numeric(ages, errors='coerce')
    print(f"Age distribution (sample):")
    print(ages.value_counts().sort_index().head(20).to_string())
    print(f"... total unique ages: {ages.nunique()}")
    print(f"  Age range: {ages.min():.0f} - {ages.max():.0f}")
    print(f"  Non-null ages: {ages.notna().sum()} / {len(ages)}")
    
    # Create age group
    def age_group(age):
        if pd.isna(age):
            return 'unknown'
        if age < 30:
            return 'young'
        elif age < 55:
            return 'middle'
        else:
            return 'old'
    
    age_groups = ages.apply(age_group)
    print(f"\nAge group distribution:")
    print(age_groups.value_counts().to_string())
    
    # Cell type
    print(f"\nCell type distribution (top 10):")
    print(adata.obs['cell_type'].value_counts().head(10).to_string())
    
    # Donor info
    print(f"\nDonor count: {adata.obs['donor_id'].nunique()}")
    print(f"Sample donors: {list(adata.obs['donor_id'].unique()[:5])}")
    
    # Save metadata as CSV
    meta = adata.obs.copy()
    meta['age_numeric'] = ages
    meta['age_group'] = age_groups
    meta.to_csv(OUTPUT_DIR / "aida_v1_metadata.csv")
    print(f"\nMetadata saved to {OUTPUT_DIR / 'aida_v1_metadata.csv'}")
    
    adata.file.close()
    return meta


def process_kedlian():
    """Process Kedlian Muscle data (smaller, can load fully)."""
    print("\n" + "=" * 60)
    print("Processing Kedlian Muscle Atlas")
    print("=" * 60)
    
    adata = sc.read_h5ad(DATA_DIR / "kedlian_muscle.h5ad")
    print(f"Loaded: {adata.shape}")
    
    print(f"Age_group distribution:")
    print(adata.obs['Age_group'].value_counts().to_string())
    print(f"\nAge_bin distribution:")
    print(adata.obs['Age_bin'].value_counts().to_string())
    
    # Parse age
    def parse_age_range(age_str):
        s = str(age_str)
        if '-' in s:
            parts = s.split('-')
            return (float(parts[0]) + float(parts[1])) / 2
        return np.nan
    
    adata.obs['age_numeric'] = adata.obs['Age_group'].apply(parse_age_range)
    print(f"\nAge numeric distribution:")
    print(adata.obs['age_numeric'].value_counts().sort_index().to_string())
    
    print(f"\nCell type distribution (top 10):")
    print(adata.obs['cell_type'].value_counts().head(10).to_string())
    
    # Save
    output_path = OUTPUT_DIR / "kedlian_muscle_processed.h5ad"
    adata.write_h5ad(output_path)
    print(f"\nSaved to {output_path}")
    
    return adata


if __name__ == "__main__":
    aida_meta = process_aida()
    kedlian = process_kedlian()
    
    print("\n" + "=" * 60)
    print("Metadata extraction complete!")
    print("Next: tAge validation on AIDA pseudobulk")
    print("=" * 60)
