#!/usr/bin/env python3
"""
Drug Module Fingerprint Pipeline (Direction C)

Compute module-specific tAge scores for LINCS L1000 drug perturbations
using multi-species module mortality clocks from Tyshkovskiy et al. 2026.

Inputs:
  - supp_table_7.xlsx (Nature supplementary): multi-species module clock coefficients
  - lincs_consensi_pert_id.tsv.bz2: LINCS consensus signatures
  - lincs_genes.tsv: LINCS gene metadata (human Entrez IDs)
  - Table_of_orthologs.csv: human → mouse ortholog mapping

Outputs:
  - results/drug_module_fingerprint.csv: composite + module-specific tAge per drug
"""

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Target drugs (default list from research notes)
# ---------------------------------------------------------------------------
DEFAULT_DRUGS = {
    "metformin": ["BRD-K79602928"],
    "sirolimus": [
        "BRD-A23770159",
        "BRD-A50287119",
        "BRD-A79768653",
        "BRD-K84937637",
        "BRD-K89626439",
    ],
    "everolimus": ["BRD-K13514097"],
    "doxorubicin": ["BRD-A52530684", "BRD-A76941896", "BRD-K92093830"],
    "dexamethasone": [
        "BRD-A10188456",
        "BRD-A35108200",
        "BRD-A69951442",
        "BRD-K00835182",
        "BRD-K30207000",
    ],
}

# Multi-species mortality module columns in supp_table_7.xlsx Sheet (C)
# Columns 18-32 (0-indexed: 17-31), but we skip col 18 (composite) for modules
MORTALITY_MODULE_COLS = {
    "blue": 3,          # Col D in Excel (0-indexed 3) -> actually need to recalculate
}

# After inspecting the Excel, mortality modules start at column index 17 (0-based)
# where index 0 = Entrez ID, 1 = Gene symbol, 2-16 = Chrono modules, 17-31 = Mortality modules
MORTALITY_MODULES = [
    ("composite", "All module genes", 17),
    ("blue", "Muscle contraction / Cytoskeleton / Glycolysis", 18),
    ("brown4", "ECM organization / EMT", 19),
    ("darkgreen", "Adaptive immunity / T cell signaling", 20),
    ("darkmagenta", "Interferon signaling", 21),
    ("darkred", "mRNA splicing", 22),
    ("darkslateblue", "Protein processing in ER / UPR", 23),
    ("green", "Cell cycle / DNA replication", 24),
    ("ivory", "Fatty acid metabolism / Peroxisome", 25),
    ("orange", "Chromatin modification", 26),
    ("pink", "Mitochondrial translation / OxPhos", 27),
    ("plum1", "Protein folding / Translation", 28),
    ("sienna3", "VEGF signaling", 29),
    ("turquoise", "Innate immunity / Inflammation", 30),
    ("white", "OxPhos / Heme metabolism", 31),
]


def load_module_coefficients(excel_path: Path) -> dict:
    """
    Load multi-species module mortality clock coefficients from
    supp_table_7.xlsx Sheet '(C) Module multi-species clocks'.

    Returns dict: {module_name: pd.Series(index=mouse_entrez, values=coeff)}
    """
    print(f"[1/5] Loading module coefficients from {excel_path}")
    # Read the Excel sheet; skip first 4 header rows, use row 5 as header
    df = pd.read_excel(
        excel_path,
        sheet_name="(C) Module multi-species clocks",
        header=None,
        skiprows=4,  # Rows 1-4 are headers
    )

    # Column 0 = Entrez ID, Column 1 = Gene symbol, Columns 2+ = module coefficients
    df = df.rename(columns={0: "entrez_id", 1: "gene_symbol"})

    # Remove intercept row
    df = df[df["entrez_id"] != "Intercept"].copy()
    df["entrez_id"] = df["entrez_id"].astype(int)
    df = df.set_index("entrez_id")

    coeffs = {}
    for mod_name, mod_annot, col_idx in MORTALITY_MODULES:
        # col_idx is 0-based for the full sheet including Entrez + symbol cols
        # But after renaming, col 0 = entrez_id, col 1 = gene_symbol, col 2 = original col 2
        # Actually df has columns 0,1,2,... where 0=entrez, 1=symbol, 2=Chrono_All...
        # So the mortality modules are at df columns: 2+15 = 17 for composite, etc.
        # Wait, let me recalculate:
        # Original Excel: col A=Entrez, B=Symbol, C=Chrono_All, D=Chrono_blue, ...
        # So 0-based: 0=Entrez, 1=Symbol, 2=Chrono_All, 3=Chrono_blue...
        # After read_excel with header=None and skiprows=4, the data starts at row 5
        # Column 0 = Entrez ID, Column 1 = Gene symbol
        # Column 2 = Chrono All, Column 3 = Chrono blue, ... Column 16 = Chrono white
        # Column 17 = Mortality All, Column 18 = Mortality blue, ... Column 31 = Mortality white
        # But after rename(columns={0: "entrez_id", 1: "gene_symbol"}),
        # the remaining columns are 2,3,4,... which correspond to original 2,3,4,...
        # So mortality composite = column 17 (original) = still column 17 in df after rename
        # Because rename only renames 0 and 1, doesn't shift indices.
        
        # Actually let me just use the column index directly
        col_name = col_idx
        ser = df[col_name].dropna()
        ser = ser[ser != 0]  # Keep only non-zero coefficients
        coeffs[mod_name] = ser
        print(f"  {mod_name:15s} ({mod_annot:45s}): {len(ser):4d} non-zero features")

    return coeffs


def load_orthologs(ortholog_path: Path) -> dict:
    """Load human → mouse Entrez ID mapping."""
    print(f"[2/5] Loading orthologs from {ortholog_path}")
    df = pd.read_csv(ortholog_path)
    # The file has columns like: Human, Mouse, Rat, etc.
    # Let's inspect first
    print(f"  Columns: {df.columns.tolist()}")
    
    # Typical format: first column might be human, second mouse
    # Let's use the column names to identify
    human_col = None
    mouse_col = None
    for c in df.columns:
        if "human" in c.lower() or c.lower() == "human_entrez":
            human_col = c
        if "mouse" in c.lower() or c.lower() == "mouse_entrez":
            mouse_col = c
    
    if human_col is None or mouse_col is None:
        # Fallback: use first two numeric columns
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if len(numeric_cols) >= 2:
            human_col = numeric_cols[0]
            mouse_col = numeric_cols[1]
        else:
            raise ValueError(f"Cannot identify human/mouse columns in {ortholog_path}")
    
    mapping = (
        df[[human_col, mouse_col]]
        .dropna()
        .drop_duplicates(subset=[human_col])
        .set_index(human_col)[mouse_col]
        .to_dict()
    )
    print(f"  {len(mapping):,} human→mouse mappings loaded")
    return mapping


def load_lincs_consensus(consensus_path: Path, pert_ids: list) -> pd.DataFrame:
    """Load LINCS consensus signatures for specified pert_ids."""
    print(f"[3/5] Loading LINCS consensus signatures from {consensus_path}")
    
    # First, get the list of available perturbagens
    print("  Scanning file for target perturbagens...")
    found = {}
    import bz2, csv
    
    with bz2.open(consensus_path, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            pid = row[0]
            if pid in pert_ids:
                found[pid] = row
                if len(found) == len(pert_ids):
                    break
    
    if not found:
        raise ValueError("No target perturbagens found in consensus file")
    
    # Build DataFrame
    data = []
    for pid in pert_ids:
        if pid in found:
            data.append(found[pid])
        else:
            print(f"  WARNING: {pid} not found in consensus file")
    
    df = pd.DataFrame(data, columns=header)
    df = df.set_index("perturbagen")
    
    # Convert to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    
    print(f"  Loaded {len(df)} signatures × {len(df.columns)} genes")
    return df


def map_human_to_mouse(sig_df: pd.DataFrame, ortholog_map: dict) -> pd.DataFrame:
    """Map human Entrez columns to mouse Entrez, aggregate duplicates by mean."""
    print(f"[4/5] Mapping human→mouse Entrez IDs")
    
    # Convert column names to int where possible
    new_cols = {}
    for col in sig_df.columns:
        try:
            h_entrez = int(col)
            m_entrez = ortholog_map.get(h_entrez)
            if m_entrez is not None:
                new_cols[col] = int(m_entrez)
            else:
                new_cols[col] = None
        except (ValueError, TypeError):
            new_cols[col] = None
    
    mapped = sig_df.rename(columns=new_cols)
    # Drop unmapped columns
    mapped = mapped.loc[:, mapped.columns.notna()]
    # Aggregate duplicate mouse entrez by mean
    mapped = mapped.T.groupby(level=0).mean().T
    
    n_before = len(sig_df.columns)
    n_after = len(mapped.columns)
    print(f"  {n_before} → {n_after} genes after ortholog mapping ({n_after/n_before*100:.1f}%)")
    return mapped


def compute_module_scores(sig_df: pd.DataFrame, module_coeffs: dict) -> pd.DataFrame:
    """Compute composite + module-specific tAge scores."""
    print(f"[5/5] Computing module-specific tAge scores")
    
    results = []
    for pert_id, sig in sig_df.iterrows():
        row = {"pert_id": pert_id}
        
        for mod_name, coeff_series in module_coeffs.items():
            # Align: keep only genes present in both signature and module coefficients
            common_genes = sig.index.intersection(coeff_series.index)
            if len(common_genes) == 0:
                score = np.nan
            else:
                score = sig.loc[common_genes].values @ coeff_series.loc[common_genes].values
            row[mod_name] = score
            
        results.append(row)
    
    return pd.DataFrame(results).set_index("pert_id")


def aggregate_by_drug(scores_df: pd.DataFrame, drug_map: dict) -> pd.DataFrame:
    """Aggregate pert_id-level scores to drug-level (mean across pert_ids)."""
    records = []
    for drug_name, pert_ids in drug_map.items():
        available = [p for p in pert_ids if p in scores_df.index]
        if not available:
            print(f"  WARNING: No data for {drug_name}")
            continue
        
        drug_scores = scores_df.loc[available].mean()
        drug_scores["drug"] = drug_name
        drug_scores["n_signatures"] = len(available)
        records.append(drug_scores)
    
    df = pd.DataFrame(records)
    df = df.set_index("drug")
    return df


def main():
    parser = argparse.ArgumentParser(description="Drug Module Fingerprint Pipeline")
    parser.add_argument("--supp-table", default="data/raw/supp_table_7.xlsx")
    parser.add_argument("--consensus", default="data/raw/lincs_consensi_pert_id.tsv.bz2")
    parser.add_argument("--genes", default="data/raw/lincs_genes.tsv")
    parser.add_argument("--orthologs", default="data/raw/tage_metadata/Table_of_orthologs.csv")
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--drugs-json", default=None, help="JSON file with drug→pert_id mapping")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load drug mapping
    if args.drugs_json:
        with open(args.drugs_json) as f:
            drug_map = json.load(f)
    else:
        drug_map = DEFAULT_DRUGS
    
    all_pert_ids = [p for pids in drug_map.values() for p in pids]

    # 1. Load module coefficients
    module_coeffs = load_module_coefficients(Path(args.supp_table))

    # 2. Load orthologs
    ortholog_map = load_orthologs(Path(args.orthologs))

    # 3. Load LINCS consensus signatures
    sig_df = load_lincs_consensus(Path(args.consensus), all_pert_ids)

    # 4. Map human→mouse
    sig_mapped = map_human_to_mouse(sig_df, ortholog_map)

    # 5. Compute module scores
    scores_df = compute_module_scores(sig_mapped, module_coeffs)

    # 6. Aggregate by drug
    drug_scores = aggregate_by_drug(scores_df, drug_map)

    # 7. Save outputs
    # Per-perturbagen scores
    scores_df.to_csv(out_dir / "module_fingerprint_per_pert_id.csv")
    print(f"\nSaved per-pert_id scores to {out_dir}/module_fingerprint_per_pert_id.csv")

    # Per-drug aggregated scores
    drug_scores.to_csv(out_dir / "module_fingerprint_per_drug.csv")
    print(f"Saved per-drug scores to {out_dir}/module_fingerprint_per_drug.csv")

    # Print summary
    print("\n" + "=" * 70)
    print("DRUG MODULE FINGERPRINT SUMMARY")
    print("=" * 70)
    
    # Sort modules for display: composite first, then modules alphabetically
    module_cols = [c for c in drug_scores.columns if c != "n_signatures"]
    module_cols_sorted = ["composite"] + sorted([c for c in module_cols if c != "composite"])
    
    display_cols = [c for c in module_cols_sorted if c in drug_scores.columns]
    
    for drug in drug_scores.index:
        print(f"\n{drug.upper()} (n={int(drug_scores.loc[drug, 'n_signatures'])} signatures)")
        composite = drug_scores.loc[drug, "composite"]
        print(f"  Composite mortality tAge: {composite:+.4f}")
        print(f"  Module contributions:")
        for mod in display_cols:
            if mod == "composite":
                continue
            val = drug_scores.loc[drug, mod]
            rel = val - composite if not pd.isna(val) else np.nan
            print(f"    {mod:15s}: {val:+.4f}  (Δ vs composite: {rel:+.4f})")

    print("\n" + "=" * 70)
    print("Pipeline complete.")


if __name__ == "__main__":
    main()
