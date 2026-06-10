"""
Pathway enrichment analysis for Geneformer ISP top hits.
Uses Enrichr API for GO/KEGG/Reactome enrichment.
"""

import pandas as pd
import numpy as np
import requests
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


ENRICHR_URL = "https://maayanlab.cloud/Enrichr"


def enrichr_submit(gene_list, description=""):
    """Submit gene list to Enrichr."""
    genes_str = "\n".join(gene_list)
    response = requests.post(
        f"{ENRICHR_URL}/addList",
        files={"list": (None, genes_str), "description": (None, description)}
    )
    if response.ok:
        return response.json()
    else:
        raise RuntimeError(f"Enrichr submit failed: {response.status_code}")


def enrichr_get_results(user_list_id, database="GO_Biological_Process_2023"):
    """Get enrichment results from Enrichr.
    Returns list of [rank, term, pvalue, oddsratio, combinedscore, genes, adj_pval, old_pval, old_adj_pval]
    """
    response = requests.get(
        f"{ENRICHR_URL}/enrich",
        params={"userListId": user_list_id, "backgroundType": database}
    )
    if response.ok:
        data = response.json()
        # Enrichr returns {database_name: [results_list]}
        if isinstance(data, dict):
            return data.get(database, [])
        return data
    else:
        raise RuntimeError(f"Enrichr get failed: {response.status_code}")


def run_enrichment(gene_list, label, databases=None):
    """Run enrichment for a gene list across multiple databases."""
    if databases is None:
        databases = {
            "GO_BP": "GO_Biological_Process_2023",
            "GO_MF": "GO_Molecular_Function_2023",
            "GO_CC": "GO_Cellular_Component_2023",
            "KEGG": "KEGG_2021_Human",
            "Reactome": "Reactome_2022",
            "WikiPathways": "WikiPathway_2023_Human",
        }
    
    print(f"\nSubmitting {label} ({len(gene_list)} genes) to Enrichr...")
    submit_res = enrichr_submit(gene_list, description=label)
    user_list_id = submit_res["userListId"]
    print(f"  List ID: {user_list_id}")
    time.sleep(2)
    
    all_results = {}
    for db_name, db_code in databases.items():
        print(f"  Querying {db_name}...")
        try:
            res = enrichr_get_results(user_list_id, db_code)
            if res:
                df = pd.DataFrame(res, columns=["rank", "term", "pvalue", "oddsratio", "combinedscore", "genes", "adjusted_pvalue", "old_pvalue", "old_adjusted_pvalue"])
                df = df.sort_values("adjusted_pvalue")
                all_results[db_name] = df
                print(f"    {len(df)} terms returned, top FDR: {df['adjusted_pvalue'].min():.2e}")
            else:
                all_results[db_name] = pd.DataFrame()
                print(f"    No results")
        except Exception as e:
            print(f"    Error: {e}")
            all_results[db_name] = pd.DataFrame()
        time.sleep(1)
    
    return all_results


def plot_enrichment(results_dict, title_prefix, output_path, max_terms=15):
    """Plot bubble chart of enrichment results."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 8))
    dbs = ["GO_BP", "KEGG", "Reactome"]
    colors = ["#E74C3C", "#3498DB", "#2ECC71"]
    
    for ax, db, color in zip(axes, dbs, colors):
        if db not in results_dict or results_dict[db].empty:
            ax.set_title(f"{db}: No significant terms")
            ax.axis("off")
            continue
        
        df = results_dict[db].head(max_terms).copy()
        df = df.sort_values("combinedscore", ascending=True)
        
        # -log10(FDR) for size
        df["neg_log_fdr"] = -np.log10(df["adjusted_pvalue"].clip(lower=1e-300))
        
        # Plot
        y_pos = np.arange(len(df))
        sizes = (df["neg_log_fdr"] / df["neg_log_fdr"].max() * 300 + 50)
        ax.scatter(df["neg_log_fdr"], y_pos, s=sizes, c=color, alpha=0.6, edgecolors="black", linewidth=0.5)
        
        # Labels
        labels = [t[:60] + "..." if len(t) > 60 else t for t in df["term"]]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("-log10(FDR)", fontsize=10)
        ax.set_title(f"{db}", fontsize=12, fontweight="bold")
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
    
    fig.suptitle(f"{title_prefix} — Pathway Enrichment", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot to {output_path}")


def save_enrichment_tables(results_dict, output_dir, prefix):
    """Save all enrichment results to CSV."""
    for db_name, df in results_dict.items():
        if not df.empty:
            path = output_dir / f"{prefix}_{db_name}_enrichment.csv"
            df.to_csv(path, index=False)
            print(f"Saved {path}")


def main():
    base_dir = Path("/home/scroll/zzhang/transcriptome-aging")
    isp_path = base_dir / "results/isp_final_ranked.csv"
    output_dir = base_dir / "results" / "enrichment"
    output_dir.mkdir(exist_ok=True)
    
    # Load ISP results
    df = pd.read_csv(isp_path)
    df = df.sort_values("mean_change", ascending=False).reset_index(drop=True)
    
    n_top = len(df) // 10
    pro_genes = df.head(n_top)["symbol"].dropna().tolist()
    anti_genes = df.tail(n_top)["symbol"].dropna().tolist()
    
    print(f"Pro-aging genes: {len(pro_genes)}")
    print(f"Anti-aging genes: {len(anti_genes)}")
    
    # Run enrichment
    print("\n" + "="*80)
    print("PRO-AGING GENE ENRICHMENT")
    print("="*80)
    pro_results = run_enrichment(pro_genes, "PBMC_pro_aging_ISP")
    save_enrichment_tables(pro_results, output_dir, "pro_aging")
    plot_enrichment(pro_results, "Pro-Aging Genes", output_dir / "pro_aging_enrichment.png")
    
    print("\n" + "="*80)
    print("ANTI-AGING GENE ENRICHMENT")
    print("="*80)
    anti_results = run_enrichment(anti_genes, "PBMC_anti_aging_ISP")
    save_enrichment_tables(anti_results, output_dir, "anti_aging")
    plot_enrichment(anti_results, "Anti-Aging Genes", output_dir / "anti_aging_enrichment.png")
    
    # Print top terms
    print("\n" + "="*80)
    print("TOP ENRICHED TERMS — PRO-AGING")
    print("="*80)
    for db in ["GO_BP", "KEGG", "Reactome"]:
        if db in pro_results and not pro_results[db].empty:
            print(f"\n{db}:")
            for _, r in pro_results[db].head(5).iterrows():
                print(f"  {r['term']:60s} | FDR={r['adjusted_pvalue']:.2e} | OR={r['oddsratio']:.2f}")
    
    print("\n" + "="*80)
    print("TOP ENRICHED TERMS — ANTI-AGING")
    print("="*80)
    for db in ["GO_BP", "KEGG", "Reactome"]:
        if db in anti_results and not anti_results[db].empty:
            print(f"\n{db}:")
            for _, r in anti_results[db].head(5).iterrows():
                print(f"  {r['term']:60s} | FDR={r['adjusted_pvalue']:.2e} | OR={r['oddsratio']:.2f}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
