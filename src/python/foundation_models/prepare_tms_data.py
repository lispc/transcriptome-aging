"""
Prepare Tabula Muris Senis data for Geneformer and scGPT training.
- Load TMS h5ad
- Filter for aging-relevant metadata (age, tissue, cell type)
- Create pseudobulk for tAge validation
- Tokenize for Geneformer
- Preprocess for scGPT
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

DATA_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data/scrna_aging/tabula_muris_senis")
OUTPUT_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_tms_data(filepath):
    """Load TMS h5ad and inspect metadata."""
    print(f"Loading {filepath}...")
    adata = sc.read_h5ad(filepath)
    print(f"Shape: {adata.shape}")
    print(f"\nObs columns: {adata.obs.columns.tolist()}")
    print(f"\nSample obs:")
    print(adata.obs.head(3).to_string())
    return adata


def extract_age_info(adata):
    """Extract age information from TMS metadata.
    
    TMS uses 'age' column with values like '3m', '18m', '21m', '24m', '30m'
    """
    if 'age' in adata.obs.columns:
        ages = adata.obs['age'].astype(str)
        print(f"\nAge distribution:")
        print(ages.value_counts().to_string())
        
        # Categorize into age groups
        def categorize_age(age_str):
            try:
                age_num = float(age_str.replace('m', '').replace('M', ''))
                if age_num <= 6:
                    return 'young'
                elif age_num <= 12:
                    return 'middle'
                else:
                    return 'old'
            except:
                return 'unknown'
        
        adata.obs['age_group'] = ages.apply(categorize_age)
        print(f"\nAge group distribution:")
        print(adata.obs['age_group'].value_counts().to_string())
    else:
        print("Warning: No 'age' column found!")
        print(f"Available columns: {adata.obs.columns.tolist()}")
    
    return adata


def filter_qc(adata, min_genes=200, min_cells=3, max_mt_pct=20):
    """Basic QC filtering."""
    print(f"\nBefore QC: {adata.shape}")
    
    # Calculate mitochondrial percentage if not present
    if 'pct_counts_mt' not in adata.obs.columns:
        adata.var['mt'] = adata.var_names.str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter
    adata = adata[adata.obs.n_genes_by_counts >= min_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :]
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After QC: {adata.shape}")
    return adata


def create_pseudobulk(adata, groupby=['age_group', 'cell_type']):
    """Create pseudobulk profiles for tAge validation."""
    print(f"\nCreating pseudobulk by {groupby}...")
    
    # Use current matrix (after QC filtering)
    expr = pd.DataFrame(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                       index=adata.obs_names, columns=adata.var_names)
    
    meta = adata.obs[groupby].copy()
    expr.index = meta.index
    
    # Aggregate
    grouped = pd.concat([meta, expr], axis=1).groupby(groupby).mean()
    
    print(f"Pseudobulk shape: {grouped.shape}")
    return grouped


def save_for_geneformer(adata, output_prefix):
    """Save processed AnnData for Geneformer tokenization."""
    output_path = OUTPUT_DIR / f"{output_prefix}.processed.h5ad"
    adata.write_h5ad(output_path)
    print(f"\nSaved processed data to {output_path}")
    return output_path


if __name__ == "__main__":
    # Process Liver data
    liver_path = DATA_DIR / "tms_liver_10x.h5ad"
    if liver_path.exists():
        adata = load_tms_data(liver_path)
        adata = extract_age_info(adata)
        adata = filter_qc(adata)
        
        # Show cell type distribution by age
        if 'cell_ontology_class' in adata.obs.columns:
            print("\nCell type by age:")
            print(pd.crosstab(adata.obs['cell_ontology_class'], adata.obs['age_group']))
        
        pseudobulk = create_pseudobulk(adata)
        pseudobulk.to_csv(OUTPUT_DIR / "tms_liver_pseudobulk.csv")
        
        save_for_geneformer(adata, "tms_liver")
    else:
        print(f"File not found: {liver_path}")
