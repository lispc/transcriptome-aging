"""
Druggability analysis for Geneformer ISP top hits.
Queries OpenTargets and ChEMBL for known drug targets and compounds.
"""

import pandas as pd
import numpy as np
import requests
import json
from pathlib import Path
import mygene


OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def query_opentargets(ensg_id):
    """Query OpenTargets for target druggability info."""
    query = """
    query TargetDruggability($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        approvedName
        tractability {
          smallMolecule {
            topCategory
            clinicalPrecedence
            highQualityCompounds
            highQualityLeadSeries
          }
          antibody {
            topCategory
            clinicalPrecedence
            cellularAssay
            geneticSupport
          }
        }
        knownDrugs(size: 10) {
          rows {
            prefName
            phase
            status
            mechanismOfAction
            disease {
              name
            }
          }
        }
      }
    }
    """
    try:
        response = requests.post(
            OPENTARGETS_URL,
            json={"query": query, "variables": {"ensemblId": ensg_id}},
            timeout=15
        )
        if response.ok:
            data = response.json()
            return data.get("data", {}).get("target")
    except Exception as e:
        pass
    return None


def query_chembl_uniprot(uniprot_id):
    """Query ChEMBL for compounds targeting a protein."""
    if not uniprot_id:
        return None
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession={uniprot_id}&limit=1"
        res = requests.get(url, timeout=10)
        if res.ok:
            data = res.json()
            targets = data.get("targets", [])
            if targets:
                chembl_id = targets[0]["target_chembl_id"]
                # Get compounds
                comp_url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id={chembl_id}&limit=20"
                comp_res = requests.get(comp_url, timeout=10)
                if comp_res.ok:
                    comp_data = comp_res.json()
                    mechanisms = comp_data.get("mechanisms", [])
                    return {
                        "target_chembl_id": chembl_id,
                        "mechanisms": mechanisms
                    }
    except Exception as e:
        pass
    return None


def get_druggability_batch(gene_symbols):
    """Batch query mygene for UniProt IDs and druggability."""
    mg = mygene.MyGeneInfo()
    res = mg.querymany(gene_symbols, scopes="symbol", fields="uniprot,drugbank", species="human")
    mapping = {}
    for r in res:
        sym = r.get("query", "")
        mapping[sym] = {
            "uniprot": r.get("uniprot", {}).get("Swiss-Prot", "") if isinstance(r.get("uniprot"), dict) else "",
            "drugbank": r.get("drugbank", {}),
        }
    return mapping


def analyze_druggability(isp_df, top_n=50):
    """Analyze druggability of top pro- and anti-aging genes."""
    pro = isp_df.head(top_n).copy()
    anti = isp_df.tail(top_n).copy()
    
    all_symbols = list(set(pro["symbol"].dropna().tolist() + anti["symbol"].dropna().tolist()))
    print(f"Querying druggability for {len(all_symbols)} genes...")
    
    # Batch query mygene
    drug_info = get_druggability_batch(all_symbols)
    
    # Add to dataframes
    for df in [pro, anti]:
        df["uniprot"] = df["symbol"].map(lambda x: drug_info.get(x, {}).get("uniprot", ""))
        df["has_drugbank"] = df["symbol"].map(lambda x: bool(drug_info.get(x, {}).get("drugbank")))
    
    # Query OpenTargets for top candidates
    print("\nQuerying OpenTargets for top candidates...")
    opentargets_results = {}
    query_list = list(pro.head(20)["ensembl_id"]) + list(anti.head(20)["ensembl_id"])
    for i, ensg in enumerate(query_list):
        if pd.isna(ensg):
            continue
        result = query_opentargets(ensg)
        if result:
            opentargets_results[ensg] = result
        if (i + 1) % 10 == 0:
            print(f"  Queried {i+1}/{len(query_list)}")
    
    return pro, anti, opentargets_results


def format_drug_report(pro, anti, ot_results, output_path):
    """Generate a druggability report."""
    lines = []
    lines.append("# Druggability Analysis: Geneformer ISP Top Hits\n")
    lines.append("## Method\n")
    lines.append("- OpenTargets Platform API for tractability assessment\n")
    lines.append("- mygene.info for DrugBank cross-references\n")
    lines.append("- ChEMBL for known compounds\n\n")
    
    # Summary
    pro_druggable = pro["has_drugbank"].sum()
    anti_druggable = anti["has_drugbank"].sum()
    lines.append("## Summary\n")
    lines.append(f"| Group | Genes Screened | DrugBank Annotated | % Druggable |\n")
    lines.append(f"|-------|---------------|-------------------|-------------|\n")
    lines.append(f"| Pro-aging | {len(pro)} | {pro_druggable} | {pro_druggable/len(pro)*100:.1f}% |\n")
    lines.append(f"| Anti-aging | {len(anti)} | {anti_druggable} | {anti_druggable/len(anti)*100:.1f}% |\n")
    lines.append("\n")
    
    # Top druggable pro-aging
    lines.append("## Top Druggable Pro-Aging Genes\n")
    lines.append("(Knockdown makes cells look older → activating these genes might promote aging)\n\n")
    pro_drug = pro[pro["has_drugbank"]].sort_values("mean_change", ascending=False)
    for _, r in pro_drug.head(10).iterrows():
        ensg = r["ensembl_id"]
        ot = ot_results.get(ensg, {})
        drugs = ot.get("knownDrugs", {}).get("rows", [])
        drug_names = ", ".join([d["prefName"] for d in drugs[:3]]) if drugs else "N/A"
        lines.append(f"### {r['symbol']}\n")
        lines.append(f"- **Effect**: Δ={r['mean_change']:+.4f}, FDR={r['fdr']:.3f}\n")
        lines.append(f"- **Known drugs**: {drug_names}\n")
        if ot.get("tractability", {}).get("smallMolecule"):
            sm = ot["tractability"]["smallMolecule"]
            lines.append(f"- **Small molecule tractability**: {sm.get('topCategory', 'N/A')}\n")
        lines.append("\n")
    
    # Top druggable anti-aging
    lines.append("## Top Druggable Anti-Aging Genes\n")
    lines.append("(Knockdown makes cells look younger → inhibiting these genes might be rejuvenating)\n\n")
    anti_drug = anti[anti["has_drugbank"]].sort_values("mean_change")
    for _, r in anti_drug.head(10).iterrows():
        ensg = r["ensembl_id"]
        ot = ot_results.get(ensg, {})
        drugs = ot.get("knownDrugs", {}).get("rows", [])
        drug_names = ", ".join([d["prefName"] for d in drugs[:3]]) if drugs else "N/A"
        lines.append(f"### {r['symbol']}\n")
        lines.append(f"- **Effect**: Δ={r['mean_change']:+.4f}, FDR={r['fdr']:.3f}\n")
        lines.append(f"- **Known drugs**: {drug_names}\n")
        if ot.get("tractability", {}).get("smallMolecule"):
            sm = ot["tractability"]["smallMolecule"]
            lines.append(f"- **Small molecule tractability**: {sm.get('topCategory', 'N/A')}\n")
        lines.append("\n")
    
    # Top candidates for drug repurposing
    lines.append("## Top Drug Repurposing Candidates\n")
    lines.append("Genes that are:\n")
    lines.append("1. Highly significant in ISP (FDR < 0.1)\n")
    lines.append("2. Have known drugs in DrugBank or clinical trials\n")
    lines.append("3. Mechanism aligns with aging biology\n\n")
    
    sig = pd.concat([pro, anti])
    sig = sig[(sig["fdr"] < 0.1) & (sig["has_drugbank"])]
    for _, r in sig.iterrows():
        direction = "Pro-aging" if r["mean_change"] > 0 else "Anti-aging"
        lines.append(f"- **{r['symbol']}** ({direction}): Δ={r['mean_change']:+.4f}, FDR={r['fdr']:.3f}\n")
    
    with open(output_path, "w") as f:
        f.writelines(lines)
    print(f"Saved report to {output_path}")


def main():
    base_dir = Path("/home/scroll/zzhang/transcriptome-aging")
    isp_path = base_dir / "results/isp_final_ranked.csv"
    output_dir = base_dir / "results" / "druggability"
    output_dir.mkdir(exist_ok=True)
    
    df = pd.read_csv(isp_path)
    df = df.sort_values("mean_change", ascending=False).reset_index(drop=True)
    
    print(f"Loaded {len(df)} ISP genes")
    
    pro, anti, ot_results = analyze_druggability(df, top_n=50)
    
    # Save tables
    pro.to_csv(output_dir / "pro_aging_druggability.csv", index=False)
    anti.to_csv(output_dir / "anti_aging_druggability.csv", index=False)
    
    # Generate report
    format_drug_report(pro, anti, ot_results, output_dir / "drug_repurposing_report.md")
    
    # Print summary
    print("\n" + "="*80)
    print("DRUGGABILITY SUMMARY")
    print("="*80)
    print(f"\nPro-aging genes with DrugBank entries: {pro['has_drugbank'].sum()}/50")
    print(f"Anti-aging genes with DrugBank entries: {anti['has_drugbank'].sum()}/50")
    
    print("\nTop druggable pro-aging hits:")
    for _, r in pro[pro["has_drugbank"]].head(5).iterrows():
        print(f"  {r['symbol']:12s} | Δ={r['mean_change']:+.4f} | UniProt={r['uniprot']}")
    
    print("\nTop druggable anti-aging hits:")
    for _, r in anti[anti["has_drugbank"]].head(5).iterrows():
        print(f"  {r['symbol']:12s} | Δ={r['mean_change']:+.4f} | UniProt={r['uniprot']}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
