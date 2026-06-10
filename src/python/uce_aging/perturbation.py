"""
In-silico gene perturbation on UCE embeddings.

For each candidate gene, zero out its expression and re-run UCE to measure
changes in aging score.
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import subprocess
import os
from pathlib import Path


def prepare_perturbed_adata(adata, gene_name, output_path):
    """Create perturbed anndata with gene expression zeroed out."""
    adata_pert = adata.copy()
    if gene_name in adata_pert.var_names:
        idx = list(adata_pert.var_names).index(gene_name)
        if hasattr(adata_pert.X, 'toarray'):
            adata_pert.X = adata_pert.X.toarray()
        adata_pert.X[:, idx] = 0
        print(f"  Zeroed out {gene_name} (idx={idx})")
    else:
        print(f"  WARNING: {gene_name} not found in var_names")
    adata_pert.write(output_path)
    return output_path


def run_uce_inference(adata_path, output_dir, model_loc, nlayers=4, output_dim=1280, batch_size=100):
    """Run UCE inference via subprocess."""
    script_dir = Path("data/uce-repo").resolve()
    cmd = [
        "python", str(script_dir / "eval_single_anndata.py"),
        "--adata_path", str(Path(adata_path).resolve()),
        "--dir", str(Path(output_dir).resolve()),
        "--species", "mouse",
        "--model_loc", str(Path(model_loc).resolve()),
        "--batch_size", str(batch_size),
        "--nlayers", str(nlayers),
        "--output_dim", str(output_dim),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(script_dir))
    if result.returncode != 0:
        print(f"  UCE error: {result.stderr[:500]}")
        return None
    
    # Find output file
    out_files = list(Path(output_dir).glob("*uce_adata.h5ad"))
    if out_files:
        return str(out_files[0])
    return None


def compute_aging_score_shift(original_adata, perturbed_path, aging_axis):
    """Compute how aging scores shift after perturbation."""
    pert_adata = sc.read_h5ad(perturbed_path)
    X_pert = pert_adata.obsm["X_uce"]
    aging_scores_pert = X_pert @ aging_axis
    aging_scores_orig = original_adata.obs["aging_score"].values
    
    delta = aging_scores_pert - aging_scores_orig
    return {
        "mean_delta": delta.mean(),
        "std_delta": delta.std(),
        "median_delta": np.median(delta),
        "young_delta": delta[original_adata.obs["age_group"] == "young"].mean(),
        "old_delta": delta[original_adata.obs["age_group"] == "old"].mean(),
        "n_cells": len(delta),
    }


def main(args):
    print(f"Loading original data: {args.input}")
    adata = sc.read_h5ad(args.input)
    
    # Add age_group
    y_age = adata.obs["age"].astype(str).str.replace("m", "").astype(int).values
    adata.obs["age_group"] = ["young" if a <= 3 else "old" for a in y_age]
    
    # Load pre-computed aging axis
    aging_axis = np.load(args.aging_axis)
    
    # Read gene list
    with open(args.genes, "r") as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"Perturbing {len(genes)} genes")
    
    results = []
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    
    for gene in genes:
        print(f"\n=== Perturbing {gene} ===")
        
        # 1. Prepare perturbed data
        pert_path = work_dir / f"perturbed_{gene}.h5ad"
        prepare_perturbed_adata(adata, gene, str(pert_path))
        
        # 2. Run UCE
        out_dir = work_dir / f"uce_{gene}"
        out_dir.mkdir(exist_ok=True)
        uce_out = run_uce_inference(
            str(pert_path), str(out_dir),
            args.model_loc, args.nlayers, args.output_dim, args.batch_size
        )
        
        if uce_out is None:
            print(f"  SKIPPED {gene} (UCE failed)")
            continue
        
        # 3. Compute shift
        shift = compute_aging_score_shift(adata, uce_out, aging_axis)
        shift["gene"] = gene
        results.append(shift)
        print(f"  Mean delta: {shift['mean_delta']:.4f} (young: {shift['young_delta']:.4f}, old: {shift['old_delta']:.4f})")
    
    # Save
    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False)
    print(f"\nResults saved to {args.output}")
    
    # Summary
    if len(df) > 0:
        print("\n=== Summary ===")
        df_sorted = df.sort_values("mean_delta")
        print("Most rejuvenating (negative delta = reverse aging):")
        print(df_sorted.head(5)[["gene", "mean_delta", "young_delta", "old_delta"]].to_string(index=False))
        print("\nMost aging-accelerating (positive delta = accelerate aging):")
        print(df_sorted.tail(5)[["gene", "mean_delta", "young_delta", "old_delta"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Original UCE-embedded h5ad")
    parser.add_argument("--aging-axis", required=True, help="Path to aging axis .npy file")
    parser.add_argument("--genes", required=True, help="File with one gene per line")
    parser.add_argument("--model-loc", default="data/uce-repo/model_files/4layer_model.torch")
    parser.add_argument("--nlayers", type=int, default=4)
    parser.add_argument("--output-dim", type=int, default=1280)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--work-dir", default="results/perturbation_work")
    parser.add_argument("--output", default="results/perturbation_results.csv")
    args = parser.parse_args()
    main(args)
