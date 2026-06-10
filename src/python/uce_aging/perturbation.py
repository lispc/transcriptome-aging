"""
In-silico perturbation on UCE embeddings.

For candidate genes, zero out expression and re-embed, measure aging score change.
NOTE: UCE issue #57 means perturbation affects other cells' embeddings.
Results should be interpreted conservatively.
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import subprocess
import os


def perturb_and_reembed(adata_path, gene, uce_dir, species="mouse", 
                        model_loc="data/uce-repo/model_files/4layer_model.torch",
                        nlayers=4, batch_size=1):
    """
    Zero out a gene's expression and re-run UCE embedding.
    Uses batch_size=1 to minimize cross-cell effects (issue #57).
    """
    adata = sc.read_h5ad(adata_path)
    
    # Check gene exists
    if gene not in adata.var_names:
        print(f"Gene {gene} not found in adata")
        return None
    
    # Create perturbed adata
    adata_pert = adata.copy()
    adata_pert[:, gene].X = 0
    
    # Save perturbed data
    pert_path = adata_path.replace(".h5ad", f"_pert_{gene}.h5ad")
    adata_pert.write(pert_path)
    
    # Run UCE on perturbed data
    out_dir = os.path.dirname(adata_path)
    cmd = [
        "python", "data/uce-repo/eval_single_anndata.py",
        "--adata_path", pert_path,
        "--dir", out_dir,
        "--species", species,
        "--model_loc", model_loc,
        "--batch_size", str(batch_size),
        "--nlayers", str(nlayers),
        "--output_dim", "1280",
    ]
    
    print(f"Running UCE on perturbed data: {pert_path}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Load perturbed embeddings
    pert_out = pert_path.replace(".h5ad", "_uce_adata.h5ad")
    if os.path.exists(pert_out):
        adata_pert_emb = sc.read_h5ad(pert_out)
        return adata_pert_emb.obsm["X_uce"]
    else:
        print(f"UCE output not found: {pert_out}")
        print(f"STDERR: {result.stderr[-500:]}")
        return None


def compute_perturbation_effect(baseline_emb, perturbed_emb, aging_axis):
    """Compute aging score change after perturbation."""
    baseline_scores = baseline_emb @ aging_axis
    perturbed_scores = perturbed_emb @ aging_axis
    
    delta = perturbed_scores - baseline_scores
    return {
        "mean_delta": delta.mean(),
        "std_delta": delta.std(),
        "n_cells": len(delta),
        "fraction_reversed": (delta < 0).mean(),  # negative = toward young
    }


def main(args):
    print(f"Loading baseline UCE embeddings from {args.input}")
    adata = sc.read_h5ad(args.input)
    baseline_emb = adata.obsm["X_uce"]
    
    # Load aging axis
    from aging_axis import compute_aging_axis
    aging_axis, _, _ = compute_aging_axis(adata)
    
    results = []
    for gene in args.genes:
        print(f"\n=== Perturbing {gene} ===")
        perturbed_emb = perturb_and_reembed(
            args.input, gene, 
            species=args.species,
            model_loc=args.model_loc,
            nlayers=args.nlayers,
            batch_size=args.batch_size
        )
        
        if perturbed_emb is not None:
            effect = compute_perturbation_effect(baseline_emb, perturbed_emb, aging_axis)
            effect["gene"] = gene
            results.append(effect)
            print(f"  Mean delta: {effect['mean_delta']:.4f}")
            print(f"  Fraction reversed: {effect['fraction_reversed']:.3f}")
    
    if results:
        df = pd.DataFrame(results)
        df.to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")
        print(df[["gene", "mean_delta", "fraction_reversed"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="UCE-embedded h5ad")
    parser.add_argument("--genes", nargs="+", required=True, help="Genes to perturb")
    parser.add_argument("--output", default="results/perturbation_results.csv")
    parser.add_argument("--species", default="mouse")
    parser.add_argument("--model-loc", default="data/uce-repo/model_files/4layer_model.torch")
    parser.add_argument("--nlayers", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=1, help="Use 1 to minimize issue #57")
    args = parser.parse_args()
    main(args)
