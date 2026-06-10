"""
Analyze Geneformer ISP results from CSV checkpoint.
- Load ISP checkpoint CSV
- Rank genes by age-classifier logit change
- Compare with known longevity databases
- Visualize top hits
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import pickle


def load_gene_mapping():
    """Load ENSG → gene symbol mapping."""
    map_path = Path("/home/scroll/zzhang/transcriptome-aging/data/human_ensembl_to_entrez.pkl")
    if map_path.exists():
        with open(map_path, "rb") as f:
            mapping = pickle.load(f)
        # mapping is ensembl_id -> (entrez_id, gene_symbol)
        return {k: v[1] for k, v in mapping.items() if isinstance(v, tuple) and len(v) > 1}
    return {}


def load_isp_checkpoint(csv_path):
    """Load ISP checkpoint CSV."""
    df = pd.read_csv(csv_path)
    # Combine young and old effects into a single score
    # Weighted average of both groups
    df["mean_change"] = (df["delta_young_mean"] + df["delta_old_mean"]) / 2
    df["std_change"] = np.sqrt((df["delta_young_std"]**2 + df["delta_old_std"]**2) / 2)
    df["n_cells"] = df["n_cells_with_gene"]
    
    # Map ENSG to gene symbol
    mapping = load_gene_mapping()
    df["gene_symbol"] = df["ensembl_id"].map(mapping)
    df["gene"] = df["gene_symbol"].fillna(df["ensembl_id"])
    
    # Statistical tests
    df["z_score"] = df["mean_change"] / (df["std_change"] / np.sqrt(df["n_cells"]))
    df["p_value"] = 2 * (1 - stats.norm.cdf(np.abs(df["z_score"])))
    
    # Benjamini-Hochberg FDR
    pvals = df["p_value"].values
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    n = len(pvals)
    fdr_corrected = np.empty(n)
    for i in range(n - 1, -1, -1):
        if i == n - 1:
            fdr_corrected[sorted_idx[i]] = sorted_pvals[i]
        else:
            fdr_corrected[sorted_idx[i]] = min(sorted_pvals[i] * n / (i + 1), fdr_corrected[sorted_idx[i + 1]])
    df["fdr"] = fdr_corrected
    
    return df.sort_values("mean_change", ascending=False)


def load_longevity_genes():
    """Load known longevity-related gene lists."""
    genage_human = [
        "IGF1R", "FOXO3", "APOE", "SIRT1", "SIRT6", "MTOR", "AMPK", "TFEB",
        "CISd2", "CLOCK", "BMAL1", "SIRT2", "SIRT3", "SIRT5", "SIRT7",
        "PARP1", "NADK2", "CARM1", "DOT1L", "KAT7", "SETD2", "MLL4",
        "EP300", "HDAC1", "HDAC2", "HDAC3", "HDAC6", "KDM6B", "KDM4B",
        "ATM", "BRCA1", "BRCA2", "TP53", "CDKN2A", "RB1", "PTEN",
        "TERT", "TERC", "WRN", "LMNA", "ZMPSTE24", "ERCC1", "ERCC2",
        "ERCC3", "ERCC4", "ERCC5", "ERCC6", "XPA", "XPC", "CSA", "CSB",
        "GDF11", "GDF15", "FGF21", "KL", "FGF23", "GHR", "GHRH", "GH1",
        "IRS1", "IRS2", "PIK3CA", "AKT1", "TSC1", "TSC2", "RHEB", "RPTOR",
        "RPS6KB1", "EIF4E", "ULK1", "ATG5", "ATG7", "BECN1", "LC3B",
        "SQSTM1", "TFEB", "TFE3", "PPARGC1A", "PPARG", "PPARA", "SIRT1",
        "NAMPT", "NMNAT1", "NMNAT2", "NMNAT3", "NADSYN1", "NADK",
        "CD38", "SARM1", "PARP1", "PARP2", "ACMSD", "IDO1", "KYNU",
        "HAAO", "NADK2", "NMRK1", "NMRK2", "NRK1", "NRK2", "NPT1",
    ]
    gwas_longevity = [
        "FOXO3", "APOE", "CDKN2B-AS1", "SH2B3", "ATXN1", "CYP2U1",
        "FBXO25", "MED1", "STK39", "ULK4", "USP35", "WRN",
    ]
    all_known = set(genage_human + gwas_longevity)
    return all_known


def annotate_known_genes(df, known_genes):
    df = df.copy()
    df["known_longevity"] = df["gene"].isin(known_genes)
    return df


def plot_top_hits(df, top_n=50, output_path="results/isp_top_hits.png"):
    fig, axes = plt.subplots(1, 2, figsize=(16, 10))
    
    pro = df.head(top_n).sort_values("mean_change")
    colors = ["red" if k else "gray" for k in pro["known_longevity"]]
    axes[0].barh(range(len(pro)), pro["mean_change"], color=colors)
    axes[0].set_yticks(range(len(pro)))
    axes[0].set_yticklabels(pro["gene"], fontsize=8)
    axes[0].set_xlabel("Mean Δ Age Score (perturbed - original)")
    axes[0].set_title(f"Top {top_n} Pro-Aging Genes (↑ age)")
    axes[0].axvline(0, color="black", linestyle="--", alpha=0.5)
    
    anti = df.tail(top_n).sort_values("mean_change", ascending=False)
    colors = ["green" if k else "gray" for k in anti["known_longevity"]]
    axes[1].barh(range(len(anti)), anti["mean_change"], color=colors)
    axes[1].set_yticks(range(len(anti)))
    axes[1].set_yticklabels(anti["gene"], fontsize=8)
    axes[1].set_xlabel("Mean Δ Age Score (perturbed - original)")
    axes[1].set_title(f"Top {top_n} Anti-Aging Genes (↓ age)")
    axes[1].axvline(0, color="black", linestyle="--", alpha=0.5)
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="red", label="Known longevity"),
        Patch(facecolor="green", label="Known longevity"),
        Patch(facecolor="gray", label="Novel"),
    ]
    fig.legend(handles=legend_elements, loc="upper center", ncol=3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot to {output_path}")


def plot_distribution(df, output_path="results/isp_distribution.png"):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    known = df[df["known_longevity"]]
    novel = df[~df["known_longevity"]]
    
    ax.hist(novel["mean_change"], bins=100, alpha=0.5, label="Novel genes", density=True)
    if len(known) > 0:
        ax.hist(known["mean_change"], bins=max(10, len(known)//2), alpha=0.8, label="Known longevity genes", density=True)
    
    ax.axvline(0, color="black", linestyle="--", alpha=0.5)
    ax.set_xlabel("Mean Δ Age Score")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of ISP Effects on Age Prediction")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Saved distribution plot to {output_path}")


def main():
    result_dir = Path("/home/scroll/zzhang/transcriptome-aging/results")
    checkpoint_path = result_dir / "isp_finetuned_geneformer_age.checkpoint.csv"
    
    if not checkpoint_path.exists():
        print("No ISP checkpoint found yet.")
        return
    
    print(f"Loading checkpoint: {checkpoint_path}")
    df = load_isp_checkpoint(checkpoint_path)
    print(f"Loaded {len(df)} genes")
    
    print("Loading known longevity genes...")
    known_genes = load_longevity_genes()
    df = annotate_known_genes(df, known_genes)
    
    print(f"\nTop 10 pro-aging genes:")
    print(df.head(10)[["gene", "mean_change", "z_score", "p_value", "fdr", "known_longevity"]].to_string())
    
    print(f"\nTop 10 anti-aging genes:")
    print(df.tail(10)[["gene", "mean_change", "z_score", "p_value", "fdr", "known_longevity"]].to_string())
    
    top10pct = len(df) // 10
    print(f"\nKnown longevity genes in top 10% (pro-aging): {df.head(top10pct)['known_longevity'].sum()}")
    print(f"Known longevity genes in bottom 10% (anti-aging): {df.tail(top10pct)['known_longevity'].sum()}")
    
    # Save ranked list
    output_csv = result_dir / "isp_ranked_genes.csv"
    df.to_csv(output_csv, index=False)
    print(f"\nSaved ranked list to {output_csv}")
    
    # Plots
    print("\nGenerating plots...")
    plot_top_hits(df, top_n=30, output_path=result_dir / "isp_top_hits.png")
    plot_distribution(df, output_path=result_dir / "isp_distribution.png")
    
    print("Done!")


if __name__ == "__main__":
    main()
