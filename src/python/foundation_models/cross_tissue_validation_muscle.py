"""
Cross-tissue validation: Compare Geneformer ISP hits from PBMC
with age-associated expression in Kedlian Muscle Atlas.

Hypothesis: Genes identified as pro-aging in PBMC ISP should be
upregulated in old muscle; anti-aging genes should be upregulated in young muscle.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_isp_results(csv_path, top_pct=10):
    """Load ISP results and return top/bottom genes."""
    df = pd.read_csv(csv_path)
    df["mean_change"] = (df["delta_young_mean"] + df["delta_old_mean"]) / 2
    df = df.sort_values("mean_change", ascending=False).reset_index(drop=True)
    
    n = len(df) * top_pct // 100
    pro_genes = df.head(n)["ensembl_id"].tolist()
    anti_genes = df.tail(n)["ensembl_id"].tolist()
    return df, pro_genes, anti_genes


def compute_muscle_logfc(adata, gene_list):
    """Compute log fold change (old / young) for given genes in muscle."""
    young = adata[adata.obs["Age_bin"] == "young", :]
    old = adata[adata.obs["Age_bin"] == "old", :]
    
    results = []
    for gene in gene_list:
        if gene not in adata.var_names:
            continue
        
        y_expr = young[:, gene].X
        o_expr = old[:, gene].X
        
        if hasattr(y_expr, "toarray"):
            y_expr = y_expr.toarray().flatten()
            o_expr = o_expr.toarray().flatten()
        else:
            y_expr = np.array(y_expr).flatten()
            o_expr = np.array(o_expr).flatten()
        
        # Log1p transform if not already
        y_mean = np.mean(y_expr)
        o_mean = np.mean(o_expr)
        
        # Pseudocount for log ratio
        logfc = np.log2((o_mean + 1) / (y_mean + 1))
        
        # Wilcoxon rank-sum test
        if len(y_expr) > 0 and len(o_expr) > 0:
            _, pval = stats.ranksums(o_expr, y_expr)
        else:
            pval = 1.0
        
        results.append({
            "ensembl_id": gene,
            "young_mean": y_mean,
            "old_mean": o_mean,
            "logFC_old_vs_young": logfc,
            "p_value": pval,
        })
    
    return pd.DataFrame(results)


def plot_validation(df_pro, df_anti, output_path):
    """Plot logFC distributions for pro- and anti-aging genes in muscle."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(df_pro["logFC_old_vs_young"], bins=30, alpha=0.6, label="PBMC pro-aging genes", color="red", density=True)
    ax.hist(df_anti["logFC_old_vs_young"], bins=30, alpha=0.6, label="PBMC anti-aging genes", color="green", density=True)
    ax.axvline(0, color="black", linestyle="--", alpha=0.5)
    ax.axvline(df_pro["logFC_old_vs_young"].median(), color="red", linestyle="-", alpha=0.8, label=f"Pro-aging median={df_pro['logFC_old_vs_young'].median():.3f}")
    ax.axvline(df_anti["logFC_old_vs_young"].median(), color="green", linestyle="-", alpha=0.8, label=f"Anti-aging median={df_anti['logFC_old_vs_young'].median():.3f}")
    
    ax.set_xlabel("Muscle log2FC (old / young)")
    ax.set_ylabel("Density")
    ax.set_title("Cross-Tissue Validation: PBMC ISP Hits in Muscle")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Saved plot to {output_path}")


def main():
    base_dir = Path("/home/scroll/zzhang/transcriptome-aging")
    isp_path = base_dir / "results/isp_finetuned_geneformer_age.csv"
    muscle_path = base_dir / "data/scrna_aging/human/kedlian_muscle.h5ad"
    output_dir = base_dir / "results"
    output_dir.mkdir(exist_ok=True)
    
    print("Loading ISP results...")
    df_isp, pro_genes, anti_genes = load_isp_results(isp_path, top_pct=10)
    print(f"Top 10% pro-aging: {len(pro_genes)} genes")
    print(f"Top 10% anti-aging: {len(anti_genes)} genes")
    
    print("\nLoading Kedlian Muscle Atlas...")
    adata = sc.read_h5ad(muscle_path)
    print(f"Muscle data: {adata.shape}")
    print(f"  Young cells: {(adata.obs['Age_bin'] == 'young').sum()}")
    print(f"  Old cells: {(adata.obs['Age_bin'] == 'old').sum()}")
    
    # Ensure gene IDs match (ISP uses ENSG, muscle may use symbols)
    print(f"\nMuscle var names (first 5): {list(adata.var_names[:5])}")
    
    # Check overlap
    muscle_genes = set(adata.var_names)
    pro_in_muscle = [g for g in pro_genes if g in muscle_genes]
    anti_in_muscle = [g for g in anti_genes if g in muscle_genes]
    print(f"Pro-aging genes found in muscle: {len(pro_in_muscle)}/{len(pro_genes)}")
    print(f"Anti-aging genes found in muscle: {len(anti_in_muscle)}/{len(anti_genes)}")
    
    print("\nComputing muscle logFC for pro-aging genes...")
    df_pro_muscle = compute_muscle_logfc(adata, pro_in_muscle)
    
    print("Computing muscle logFC for anti-aging genes...")
    df_anti_muscle = compute_muscle_logfc(adata, anti_in_muscle)
    
    # Add gene symbols
    sym_map = dict(zip(df_isp["ensembl_id"], df_isp.get("symbol", df_isp["ensembl_id"])))
    df_pro_muscle["symbol"] = df_pro_muscle["ensembl_id"].map(sym_map)
    df_anti_muscle["symbol"] = df_anti_muscle["ensembl_id"].map(sym_map)
    
    # Summary statistics
    print(f"\n{'='*80}")
    print("VALIDATION SUMMARY")
    print(f"{'='*80}")
    
    pro_median = df_pro_muscle["logFC_old_vs_young"].median()
    anti_median = df_anti_muscle["logFC_old_vs_young"].median()
    
    print(f"\nPro-aging genes (PBMC) in muscle:")
    print(f"  Median logFC (old/young): {pro_median:+.4f}")
    print(f"  % up in old (logFC>0): {(df_pro_muscle['logFC_old_vs_young'] > 0).mean()*100:.1f}%")
    print(f"  Wilcoxon p (pro vs 0): {stats.wilcoxon(df_pro_muscle['logFC_old_vs_young']).pvalue:.4f}")
    
    print(f"\nAnti-aging genes (PBMC) in muscle:")
    print(f"  Median logFC (old/young): {anti_median:+.4f}")
    print(f"  % up in young (logFC<0): {(df_anti_muscle['logFC_old_vs_young'] < 0).mean()*100:.1f}%")
    print(f"  Wilcoxon p (anti vs 0): {stats.wilcoxon(df_anti_muscle['logFC_old_vs_young']).pvalue:.4f}")
    
    # Mann-Whitney U test between pro and anti
    u_stat, u_pval = stats.mannwhitneyu(df_pro_muscle["logFC_old_vs_young"], df_anti_muscle["logFC_old_vs_young"], alternative="two-sided")
    print(f"\nMann-Whitney U (pro vs anti): p={u_pval:.4f}")
    
    # Consistency score
    # Pro-aging genes should have positive logFC (up in old)
    # Anti-aging genes should have negative logFC (up in young)
    pro_consistent = (df_pro_muscle["logFC_old_vs_young"] > 0).sum()
    anti_consistent = (df_anti_muscle["logFC_old_vs_young"] < 0).sum()
    total_consistent = pro_consistent + anti_consistent
    total = len(df_pro_muscle) + len(df_anti_muscle)
    print(f"\nConsistency rate: {total_consistent}/{total} = {total_consistent/total*100:.1f}%")
    print(f"  (Pro-aging up in old: {pro_consistent}/{len(df_pro_muscle)})")
    print(f"  (Anti-aging up in young: {anti_consistent}/{len(df_anti_muscle)})")
    
    # Print top consistent genes
    print(f"\n{'='*80}")
    print("Top pro-aging genes UP in old muscle")
    print(f"{'='*80}")
    top_pro = df_pro_muscle.nlargest(10, "logFC_old_vs_young")
    for _, r in top_pro.iterrows():
        print(f"  {r['symbol']:12s} | logFC={r['logFC_old_vs_young']:+.3f} | p={r['p_value']:.3f}")
    
    print(f"\n{'='*80}")
    print("Top anti-aging genes UP in young muscle")
    print(f"{'='*80}")
    top_anti = df_anti_muscle.nsmallest(10, "logFC_old_vs_young")
    for _, r in top_anti.iterrows():
        print(f"  {r['symbol']:12s} | logFC={r['logFC_old_vs_young']:+.3f} | p={r['p_value']:.3f}")
    
    # Save results
    df_pro_muscle.to_csv(output_dir / "muscle_validation_pro_aging.csv", index=False)
    df_anti_muscle.to_csv(output_dir / "muscle_validation_anti_aging.csv", index=False)
    
    # Plot
    plot_validation(df_pro_muscle, df_anti_muscle, output_dir / "muscle_validation.png")
    
    print(f"\nDone! Files saved to {output_dir}")


if __name__ == "__main__":
    main()
