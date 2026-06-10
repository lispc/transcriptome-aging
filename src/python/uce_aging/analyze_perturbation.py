"""
Analyze perturbation results: compute aging score shifts for all perturbed genes.
"""
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path


def main():
    # Load original subsampled data with aging scores from full data
    orig_full = sc.read_h5ad("results/uce_spleen/tms_spleen_uce_adata_aging_axis.h5ad")
    orig_sub = sc.read_h5ad("results/uce_spleen/tms_spleen_subsampled_1k.h5ad")
    
    # Align by cell ID
    cell_ids = orig_sub.obs["cell"].values
    mask = np.isin(orig_full.obs["cell"].values, cell_ids)
    orig_scores = orig_full.obs.loc[mask, "aging_score"].values
    orig_cell_order = {c: i for i, c in enumerate(orig_full.obs.loc[mask, "cell"].values)}
    
    # Age group
    y_age = orig_sub.obs["age"].astype(str).str.replace("m", "").astype(int).values
    
    axis = np.load("results/uce_spleen/aging_axis.npy")
    
    work_dir = Path("results/perturbation_work")
    gene_list = "results/uce_spleen/perturbation_genes.txt"
    
    with open(gene_list) as f:
        genes = [l.strip() for l in f if l.strip()]
    
    results = []
    for gene in genes:
        # Try both path formats (UCE sometimes writes to subdir, sometimes flat)
        uce_file = work_dir / f"uce_{gene}" / f"perturbed_{gene}_uce_adata.h5ad"
        if not uce_file.exists():
            uce_file = work_dir / f"uce_{gene}perturbed_{gene}_uce_adata.h5ad"
        if not uce_file.exists():
            print(f"SKIP {gene}: file not found")
            continue
        
        pert = sc.read_h5ad(uce_file)
        pert_scores = pert.obsm["X_uce"] @ axis
        
        # Align
        pert_order = [orig_cell_order[c] for c in pert.obs["cell"].values]
        orig_scores_aligned = orig_scores[pert_order]
        
        delta = pert_scores - orig_scores_aligned
        
        results.append({
            "gene": gene,
            "mean_delta": delta.mean(),
            "median_delta": np.median(delta),
            "std_delta": delta.std(),
            "young_delta": delta[:500].mean(),
            "old_delta": delta[500:].mean(),
            "n_cells": len(delta),
        })
        
        print(f"{gene}: mean_delta={delta.mean():.4f} (young={delta[:500].mean():.4f}, old={delta[500:].mean():.4f})")
    
    df = pd.DataFrame(results)
    df.to_csv("results/perturbation_results.csv", index=False)
    
    print("\n=== SUMMARY ===")
    df_sorted = df.sort_values("mean_delta")
    print("\nMost rejuvenating (KO reduces aging score):")
    print(df_sorted.head(10)[["gene", "mean_delta", "young_delta", "old_delta"]].to_string(index=False))
    print("\nMost aging-accelerating (KO increases aging score):")
    print(df_sorted.tail(10)[["gene", "mean_delta", "young_delta", "old_delta"]].to_string(index=False))
    
    # Direction validation
    spleen_corr = pd.read_csv("results/uce_spleen/gene_aging_correlation.csv")
    merged = df.merge(spleen_corr, on="gene")
    concordant = ((merged["corr_with_aging_score"] > 0) & (merged["mean_delta"] < 0)).sum() + \
                ((merged["corr_with_aging_score"] < 0) & (merged["mean_delta"] > 0)).sum()
    print(f"\nDirection concordance: {concordant}/{len(merged)} ({concordant/len(merged)*100:.0f}%)")
    print("(Pro-aging genes should have negative delta when KO'd; anti-aging genes should have positive delta)")


if __name__ == "__main__":
    main()
