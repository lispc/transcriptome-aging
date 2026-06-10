"""
Prepare AIDA data for Geneformer tokenization.

Steps:
1. Stratified sampling by age group from AIDA v1 (backed mode)
2. Denormalize log1p(CPM) back to approximate raw counts
3. Add required columns: ensembl_id (var), n_counts (obs), age_group (obs)
4. Save as new h5ad for Geneformer tokenizer
"""

import argparse
import json
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path


def parse_age(development_stage: str) -> float:
    """Parse age from development_stage string like '33-year-old stage'."""
    import re
    m = re.search(r'(\d+)-year-old', str(development_stage))
    if m:
        return float(m.group(1))
    return np.nan


def assign_age_group(age: float) -> str:
    if pd.isna(age):
        return 'unknown'
    if age < 30:
        return 'young'
    elif age <= 55:
        return 'middle'
    else:
        return 'old'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='/home/scroll/zzhang/transcriptome-aging/data/scrna_aging/human/aida_v1.h5ad')
    parser.add_argument('--output', default='/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_for_geneformer.h5ad')
    parser.add_argument('--metadata', default='/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_metadata.csv')
    parser.add_argument('--cells-per-group', type=int, default=20000, help='Cells per age group')
    parser.add_argument('--chunk-size', type=int, default=5000)
    args = parser.parse_args()

    print("Loading metadata...")
    meta = pd.read_csv(args.metadata, index_col=0, low_memory=False)
    
    # Ensure age_group exists
    if 'age_group' not in meta.columns:
        meta['age_group'] = meta['age_numeric'].apply(assign_age_group)
    
    # Stratified sampling
    sampled_indices = []
    for group in ['young', 'middle', 'old']:
        group_meta = meta[meta['age_group'] == group]
        n = min(args.cells_per_group, len(group_meta))
        sampled = group_meta.sample(n=n, random_state=42)
        sampled_indices.extend(sampled.index.tolist())
        print(f"  {group}: {n} cells")
    
    print(f"Total sampled: {len(sampled_indices)}")
    
    # Open backed AIDA
    print("Opening AIDA (backed mode)...")
    adata_full = sc.read_h5ad(args.input, backed='r')
    
    # Create position map
    all_obs_names = list(adata_full.obs_names)
    pos_map = {name: idx for idx, name in enumerate(all_obs_names)}
    positions = [pos_map[idx] for idx in sampled_indices if idx in pos_map]
    print(f"Valid positions: {len(positions)}")
    
    # Process in chunks to avoid memory issues
    chunks = []
    obs_chunks = []
    
    for start in range(0, len(positions), args.chunk_size):
        end = min(start + args.chunk_size, len(positions))
        chunk_pos = positions[start:end]
        chunk_idx = sampled_indices[start:end]
        
        print(f"  Processing chunk {start//args.chunk_size + 1}/{(len(positions)-1)//args.chunk_size + 1}: {len(chunk_pos)} cells")
        
        # Read expression
        expr = adata_full.X[chunk_pos, :]
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        
        # Denormalize: log1p(counts / n_counts * 10000) -> counts
        ncounts = adata_full.obs['nCount_RNA'].values[chunk_pos]
        ncounts = ncounts.reshape(-1, 1)
        
        # raw_counts = (exp(X) - 1) * n_counts / 10000
        raw_counts = np.expm1(expr) * ncounts / 10000.0
        # Round to nearest integer and clip
        raw_counts = np.clip(np.round(raw_counts), 0, None).astype(np.float32)
        
        chunks.append(raw_counts)
        
        # Build obs chunk
        chunk_obs = meta.loc[chunk_idx].copy()
        chunk_obs['n_counts'] = adata_full.obs['nCount_RNA'].values[chunk_pos]
        obs_chunks.append(chunk_obs)
    
    print("Concatenating...")
    X_all = np.vstack(chunks)
    obs_all = pd.concat(obs_chunks)
    
    # Build var
    var_df = pd.DataFrame(index=adata_full.var.index)
    var_df['ensembl_id'] = adata_full.var.index
    
    # Create AnnData
    adata_new = sc.AnnData(X=X_all, obs=obs_all, var=var_df)
    
    # Ensure ensembl_id is string
    adata_new.var['ensembl_id'] = adata_new.var['ensembl_id'].astype(str)
    adata_new.obs['n_counts'] = adata_new.obs['n_counts'].astype(float)
    
    print(f"New adata shape: {adata_new.shape}")
    print(f"X dtype: {adata_new.X.dtype}")
    print(f"X max: {adata_new.X.max():.0f}")
    print(f"X mean: {adata_new.X.mean():.2f}")
    
    # Save
    print(f"Saving to {args.output}...")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata_new.write_h5ad(args.output)
    print("Done!")
    
    adata_full.file.close()


if __name__ == '__main__':
    main()
