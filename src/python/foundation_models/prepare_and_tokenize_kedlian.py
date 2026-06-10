"""
Prepare Kedlian Muscle Atlas for Geneformer and run tokenization.
"""

import json
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import sys
sys.path.insert(0, '/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2')
from geneformer import TranscriptomeTokenizer


def main():
    input_path = "/home/scroll/zzhang/transcriptome-aging/data/scrna_aging/human/kedlian_muscle.h5ad"
    output_h5ad = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/kedlian_muscle_for_geneformer.h5ad"
    output_dataset = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/kedlian_muscle_tokenized"

    print("Loading Kedlian Muscle data (backed)...")
    adata_full = sc.read_h5ad(input_path, backed='r')
    print(f"Full shape: {adata_full.shape}")

    # Stratified sampling by age bin
    print("Sampling cells...")
    n_per_group = 15000
    sampled_indices = []
    for age_bin in ["young", "old"]:
        mask = adata_full.obs["Age_bin"] == age_bin
        indices = np.where(mask)[0]
        n = min(n_per_group, len(indices))
        chosen = np.random.choice(indices, n, replace=False)
        sampled_indices.extend(chosen.tolist())
        print(f"  {age_bin}: {n} cells")

    print(f"Total sampled: {len(sampled_indices)}")

    # Read expression
    print("Reading expression...")
    expr = adata_full.X[sampled_indices, :]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()

    # Denormalize
    ncounts = adata_full.obs["n_counts"].values[sampled_indices]
    ncounts = ncounts.reshape(-1, 1)
    raw_counts = np.expm1(expr) * ncounts / 10000.0
    raw_counts = np.clip(np.round(raw_counts), 0, None).astype(np.float32)

    # Build obs
    obs_df = pd.DataFrame(index=adata_full.obs_names[sampled_indices])
    for col in ["Age_group", "Age_bin", "cell_type", "donor_id", "n_counts"]:
        if col in adata_full.obs.columns:
            obs_df[col] = adata_full.obs[col].values[sampled_indices]

    # Build var
    var_df = pd.DataFrame(index=adata_full.var.index)
    var_df["ensembl_id"] = adata_full.var.index.astype(str)

    # Create AnnData
    adata_new = sc.AnnData(X=raw_counts, obs=obs_df, var=var_df)
    adata_new.obs["n_counts"] = adata_new.obs["n_counts"].astype(float)

    print(f"New adata shape: {adata_new.shape}")
    print(f"X max: {adata_new.X.max():.0f}")

    # Save intermediate h5ad
    print(f"Saving to {output_h5ad}...")
    Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata_new.write_h5ad(output_h5ad)
    print("Saved.")

    adata_full.file.close()

    # Convert to sparse for tokenizer
    print("Converting to sparse...")
    import scipy.sparse as sp
    adata_tok = sc.read_h5ad(output_h5ad)
    adata_tok.X = sp.csr_matrix(adata_tok.X)

    # Remove ensembl_id_collapsed if exists
    if "ensembl_id_collapsed" in adata_tok.var.columns:
        del adata_tok.var["ensembl_id_collapsed"]

    adata_tok.write_h5ad(output_h5ad)
    print("Converted to sparse and saved.")

    # Tokenize
    print("Tokenizing for Geneformer...")
    custom_attrs = {
        "Age_bin": "age_bin",
        "Age_group": "age_group",
        "cell_type": "cell_type",
        "donor_id": "donor_id",
    }

    tk = TranscriptomeTokenizer(
        custom_attr_name_dict=custom_attrs,
        nproc=4,
        model_version="V2",
        use_h5ad_index=True,
    )

    tk.tokenize_data(
        data_directory=Path(output_h5ad).parent,
        output_directory=Path(output_dataset).parent,
        output_prefix="kedlian_muscle_tokenized",
        file_format="h5ad",
        input_identifier="kedlian_muscle_for_geneformer",
    )

    print(f"Tokenized dataset saved to {output_dataset}.dataset")


if __name__ == "__main__":
    np.random.seed(42)
    main()
