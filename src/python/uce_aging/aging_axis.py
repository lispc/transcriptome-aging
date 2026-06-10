"""
Aging axis analysis on UCE embeddings.

Compute young vs old mean embedding difference, project cells onto aging axis.
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def compute_aging_axis(adata, young_cutoff=3, old_cutoff=18):
    """Compute aging axis = mean(old) - mean(young)."""
    X = adata.obsm["X_uce"]
    y_age = adata.obs["age"].astype(str).str.replace("m", "").astype(int).values
    
    young_mask = y_age <= young_cutoff
    old_mask = y_age >= old_cutoff
    
    young_mean = X[young_mask].mean(axis=0)
    old_mean = X[old_mask].mean(axis=0)
    
    aging_axis = old_mean - young_mean
    aging_axis = aging_axis / np.linalg.norm(aging_axis)
    
    # Project all cells onto aging axis
    aging_scores = X @ aging_axis
    adata.obs["aging_score"] = aging_scores
    
    return aging_axis, young_mean, old_mean


def plot_aging_distributions(adata, outdir="results/figures"):
    """Plot aging score distributions by age, tissue, cell type."""
    import os
    os.makedirs(outdir, exist_ok=True)
    
    tissue = adata.obs["tissue"].iloc[0] if "tissue" in adata.obs.columns else "unknown"
    
    # 1. Aging score by age
    fig, ax = plt.subplots(figsize=(10, 6))
    ages = sorted(adata.obs["age"].unique())
    data_to_plot = [adata.obs[adata.obs["age"] == a]["aging_score"].values for a in ages]
    bp = ax.boxplot(data_to_plot, labels=ages, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    ax.set_xlabel("Age")
    ax.set_ylabel("Aging Score")
    ax.set_title(f"Aging Score Distribution by Age - {tissue}")
    plt.tight_layout()
    plt.savefig(f"{outdir}/aging_score_by_age_{tissue}.png", dpi=150)
    plt.close()
    
    # 2. Aging score by cell type (top 10)
    if "cell_ontology_class" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 8))
        top_cts = adata.obs["cell_ontology_class"].value_counts().head(10).index
        sub = adata[adata.obs["cell_ontology_class"].isin(top_cts)]
        
        df_plot = sub.obs[["cell_ontology_class", "aging_score", "age"]].copy()
        df_plot["age_group"] = df_plot["age"].apply(lambda x: "young" if int(x.replace("m","")) <= 3 else "old")
        
        sns.boxplot(data=df_plot, x="cell_ontology_class", y="aging_score", hue="age_group", ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_title(f"Aging Score by Cell Type - {tissue}")
        plt.tight_layout()
        plt.savefig(f"{outdir}/aging_score_by_celltype_{tissue}.png", dpi=150)
        plt.close()
    
    print(f"Plots saved to {outdir}/")


def main(args):
    print(f"Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    
    aging_axis, young_mean, old_mean = compute_aging_axis(adata)
    
    tissue = adata.obs["tissue"].iloc[0] if "tissue" in adata.obs.columns else "unknown"
    print(f"\n=== Aging Axis Analysis: {tissue} ===")
    print(f"Aging axis norm: {np.linalg.norm(aging_axis):.4f}")
    
    # Statistics by age
    ages = sorted(adata.obs["age"].unique())
    for age in ages:
        mask = adata.obs["age"] == age
        scores = adata.obs.loc[mask, "aging_score"]
        print(f"  {age}: mean={scores.mean():.3f}, std={scores.std():.3f}, n={mask.sum()}")
    
    # Save annotated adata
    if args.output:
        adata.write(args.output)
        print(f"\nSaved to {args.output}")
    
    # Plot
    if args.plot:
        plot_aging_distributions(adata, outdir=args.plot_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="UCE-embedded h5ad")
    parser.add_argument("--output", default=None, help="Output path")
    parser.add_argument("--plot", action="store_true", help="Generate plots")
    parser.add_argument("--plot-dir", default="results/figures", help="Plot output dir")
    args = parser.parse_args()
    main(args)
