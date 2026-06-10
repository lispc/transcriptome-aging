"""
Prepare AIDA data for scGPT.

scGPT uses gene symbols as vocabulary (not ENSG).
We need to:
1. Map AIDA's ENSG IDs to gene symbols
2. Filter to genes present in scGPT vocab
3. Save in scGPT-compatible format
"""

import json
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path


def main():
    # Paths
    aida_path = "/home/scroll/zzhang/transcriptome-aging/data/scrna_aging/human/aida_v1.h5ad"
    meta_path = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_metadata.csv"
    scgpt_vocab_path = "/home/scroll/zzhang/transcriptome-aging/data/scgpt/weights/checkpoint-0623/scGPT_human/vocab.json"
    output_path = "/home/scroll/zzhang/transcriptome-aging/data/scgpt/aida_v1_for_scgpt.h5ad"

    # Load scGPT vocab
    print("Loading scGPT vocab...")
    with open(scgpt_vocab_path) as f:
        vocab = json.load(f)
    print(f"Vocab size: {len(vocab)}")

    # Load metadata (for cell sampling)
    print("Loading metadata...")
    meta = pd.read_csv(meta_path, index_col=0, low_memory=False)

    # Stratified sampling: 10K per age group = 30K total
    sampled_indices = []
    for group in ["young", "middle", "old"]:
        group_meta = meta[meta["age_group"] == group]
        n = min(10000, len(group_meta))
        sampled = group_meta.sample(n=n, random_state=42)
        sampled_indices.extend(sampled.index.tolist())
    print(f"Sampled {len(sampled_indices)} cells")

    # Open AIDA backed
    print("Opening AIDA (backed mode)...")
    adata_full = sc.read_h5ad(aida_path, backed="r")

    # Get gene symbols
    symbols = list(adata_full.var["feature_name"])
    ensg_ids = list(adata_full.var.index)

    # Find genes in scGPT vocab
    mask = [s in vocab for s in symbols]
    n_in_vocab = sum(mask)
    print(f"Genes in scGPT vocab: {n_in_vocab} / {len(symbols)}")

    # Get positions of sampled cells
    all_obs_names = list(adata_full.obs_names)
    pos_map = {name: idx for idx, name in enumerate(all_obs_names)}
    positions = [pos_map[idx] for idx in sampled_indices if idx in pos_map]

    # Read expression for sampled cells, filtered genes
    print("Reading expression data...")
    expr = adata_full.X[positions, :]
    if hasattr(expr, "toarray"):
        expr = expr.toarray()

    # Filter to genes in vocab
    expr_filtered = expr[:, mask]
    symbols_filtered = [s for s, m in zip(symbols, mask) if m]

    print(f"Filtered expression shape: {expr_filtered.shape}")

    # Build obs
    obs_df = meta.loc[sampled_indices].copy()

    # Build var
    var_df = pd.DataFrame(index=symbols_filtered)

    # Create AnnData
    adata_new = sc.AnnData(X=expr_filtered, obs=obs_df, var=var_df)

    # Save
    print(f"Saving to {output_path}...")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata_new.write_h5ad(output_path)
    print("Done!")

    adata_full.file.close()


if __name__ == "__main__":
    main()
