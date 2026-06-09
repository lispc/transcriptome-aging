#!/usr/bin/env python3
"""
Drug Module Fingerprint Pipeline v2 — Composite Clock Decomposition

Uses the validated multi-species mortality composite clock (Zenodo .pkl)
and decomposes its prediction into module-level contributions based on
WGCNA module gene lists from Tyshkovskiy et al. 2026 Nature.

Method:
  1. Load EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl
  2. Extract composite coefficients + preprocessing params (imputer, scaler)
  3. Load WGCNA module gene lists from supp_table_7.xlsx Sheet (C)
  4. Load LINCS L1000 consensus signatures
  5. Map human→mouse Entrez, align to clock feature space
  6. Apply median imputation + centering (using training-set params)
  7. Compute composite tAge and per-module contributions
  8. Aggregate by drug and output fingerprints

Outputs:
  - results/module_fingerprint_v2_per_pert_id.csv
  - results/module_fingerprint_v2_per_drug.csv
  - figures/module_fingerprint_heatmap.pdf
"""

import argparse
import json
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

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

# Mortality module columns in supp_table_7.xlsx Sheet (C)
# col 17 = composite (All module genes), 18-31 = individual modules
MORTALITY_MODULES = [
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


def load_clock_model(pkl_path: Path):
    """Load .pkl clock model and extract coefficients + preprocessing params."""
    print(f"[1/6] Loading clock model from {pkl_path}")
    model = joblib.load(pkl_path)

    imputer = model.named_steps["imputation"]
    scaler = model.named_steps["scaler"]
    estimator = model.named_steps["estimator"]

    params = {
        "feature_names": list(model.feature_names_in_),
        "coef": estimator.coef_,
        "intercept": float(estimator.intercept_),
        "impute_stats": imputer.statistics_,
        "center_mean": scaler.mean_,
    }
    print(
        f"  Features: {len(params['feature_names'])}, "
        f"non-zero coefs: {(params['coef'] != 0).sum()}, "
        f"intercept: {params['intercept']:.4f}"
    )
    return params


def load_module_gene_lists(excel_path: Path) -> dict:
    """Load module gene lists from supp_table_7.xlsx Sheet (C)."""
    print(f"[2/6] Loading module gene lists from {excel_path}")
    df = pd.read_excel(
        excel_path,
        sheet_name="(C) Module multi-species clocks",
        header=None,
        skiprows=4,
    )
    df = df[df[0] != "Intercept"].copy()
    df[0] = df[0].astype(int)
    df = df.set_index(0)

    modules = {}
    for mod_name, mod_annot, col_idx in MORTALITY_MODULES:
        genes = set(df[df[col_idx] != 0].index.astype(str))
        modules[mod_name] = genes
        print(f"  {mod_name:15s} ({mod_annot:45s}): {len(genes):4d} genes")
    return modules


def load_orthologs(ortholog_path: Path) -> dict:
    """Load human → mouse Entrez ID mapping."""
    print(f"[3/6] Loading orthologs from {ortholog_path}")
    df = pd.read_csv(ortholog_path)
    mapping = (
        df[["Entrez.Human", "Entrez.Mouse"]]
        .dropna()
        .drop_duplicates(subset=["Entrez.Human"])
        .set_index("Entrez.Human")["Entrez.Mouse"]
        .astype(int)
        .to_dict()
    )
    print(f"  {len(mapping):,} human→mouse mappings loaded")
    return mapping


def load_lincs_consensus(consensus_path: Path, pert_ids: list) -> pd.DataFrame:
    """Load LINCS consensus signatures for specified pert_ids."""
    print(f"[4/6] Loading LINCS consensus signatures from {consensus_path}")

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

    data = []
    missing = []
    for pid in pert_ids:
        if pid in found:
            data.append(found[pid])
        else:
            missing.append(pid)

    if missing:
        print(f"  WARNING: {len(missing)} perturbagen(s) not found: {missing}")

    df = pd.DataFrame(data, columns=header)
    df = df.set_index("perturbagen")
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    print(f"  Loaded {len(df)} signatures × {len(df.columns)} genes")
    return df


def align_and_preprocess(
    sig_df: pd.DataFrame,
    ortholog_map: dict,
    clock_params: dict,
) -> pd.DataFrame:
    """Map human→mouse, align to clock features, impute, and center."""
    print(f"[5/6] Aligning signatures to clock feature space")

    # Map columns (human Entrez → mouse Entrez)
    new_cols = {}
    for col in sig_df.columns:
        try:
            h = int(col)
            m = ortholog_map.get(h)
            new_cols[col] = str(int(m)) if m is not None else None
        except (ValueError, TypeError):
            new_cols[col] = None

    mapped = sig_df.rename(columns=new_cols)
    mapped = mapped.loc[:, mapped.columns.notna()]
    mapped = mapped.T.groupby(level=0).mean().T

    n_before = len(sig_df.columns)
    n_after = len(mapped.columns)
    print(f"  {n_before} → {n_after} genes after ortholog mapping")

    # Align to clock feature space
    feature_names = clock_params["feature_names"]
    aligned = pd.DataFrame(index=mapped.index, columns=feature_names, dtype=float)

    # Fill available genes
    common = mapped.columns.intersection(feature_names)
    aligned[common] = mapped[common]
    print(f"  {len(common)} / {len(feature_names)} clock features present in data")

    # Impute missing genes using training-set medians
    missing = aligned.columns[aligned.isna().any()]
    impute_stats = clock_params["impute_stats"]
    feature_to_idx = {f: i for i, f in enumerate(feature_names)}
    for col in missing:
        idx = feature_to_idx[col]
        aligned[col] = aligned[col].fillna(impute_stats[idx])

    # Center using training-set means
    center_mean = clock_params["center_mean"]
    centered = aligned.copy()
    for col in centered.columns:
        idx = feature_to_idx[col]
        centered[col] = centered[col] - center_mean[idx]

    print(f"  Preprocessing complete: median imputation + centering")
    return centered


def compute_fingerprint(
    centered_df: pd.DataFrame,
    clock_params: dict,
    module_genes: dict,
) -> pd.DataFrame:
    """Compute composite tAge and per-module contributions."""
    print(f"[6/6] Computing module fingerprints")

    coef = clock_params["coef"]
    intercept = clock_params["intercept"]
    feature_names = clock_params["feature_names"]
    feature_to_idx = {f: i for i, f in enumerate(feature_names)}

    results = []
    for pert_id, row in centered_df.iterrows():
        x = row.values

        # Composite prediction
        composite = intercept + np.dot(coef, x)

        # Module contributions
        module_scores = {"pert_id": pert_id, "composite": composite}

        # Non-module contribution (all genes not in any module)
        all_module_genes = set()
        for genes in module_genes.values():
            all_module_genes.update(genes)

        non_module_mask = np.ones(len(feature_names), dtype=bool)
        for g in all_module_genes:
            if g in feature_to_idx:
                non_module_mask[feature_to_idx[g]] = False
        non_module_contrib = np.dot(coef[non_module_mask], x[non_module_mask])
        module_scores["non_module"] = non_module_contrib

        # Per-module contributions
        for mod_name, genes in module_genes.items():
            indices = [feature_to_idx[g] for g in genes if g in feature_to_idx]
            if indices:
                contrib = np.dot(coef[indices], x[indices])
                module_scores[mod_name] = contrib
            else:
                module_scores[mod_name] = np.nan

        results.append(module_scores)

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


def plot_fingerprint(drug_scores: pd.DataFrame, out_dir: Path):
    """Create heatmap of module fingerprints across drugs."""
    print("\n[7/7] Generating figures")
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Prepare data for heatmap: module columns only (exclude composite, non_module, n_signatures)
    module_cols = [c for c in drug_scores.columns if c not in ("composite", "non_module", "n_signatures")]
    module_cols = sorted(module_cols)

    plot_df = drug_scores[module_cols].copy()

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(plot_df.values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)

    ax.set_xticks(range(len(module_cols)))
    ax.set_xticklabels(module_cols, rotation=45, ha="right")
    ax.set_yticks(range(len(plot_df.index)))
    ax.set_yticklabels(plot_df.index)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Module contribution to mortality tAge")

    # Add text annotations
    for i in range(len(plot_df.index)):
        for j in range(len(module_cols)):
            val = plot_df.iloc[i, j]
            text_color = "white" if abs(val) > 1 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=text_color, fontsize=8)

    ax.set_title("Drug Module Fingerprint (Mortality Clock Decomposition)")
    plt.tight_layout()
    fig_path = fig_dir / "module_fingerprint_heatmap.pdf"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"  Saved heatmap to {fig_path}")

    # Bar plot: composite tAge per drug
    fig, ax = plt.subplots(figsize=(8, 5))
    drugs = drug_scores.index
    composites = drug_scores["composite"]
    colors = ["green" if v < 0 else "red" for v in composites]
    ax.barh(drugs, composites, color=colors, edgecolor="black")
    ax.axvline(x=0, color="black", linestyle="-", linewidth=0.5)
    ax.set_xlabel("Composite Mortality tAge")
    ax.set_title("Drug-Induced Mortality tAge Shift (LINCS L1000)")
    plt.tight_layout()
    fig_path = fig_dir / "composite_tage_barplot.pdf"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"  Saved composite barplot to {fig_path}")


def main():
    parser = argparse.ArgumentParser(description="Drug Module Fingerprint Pipeline v2")
    parser.add_argument("--pkl-model", default="models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    parser.add_argument("--supp-table", default="data/raw/supp_table_7.xlsx")
    parser.add_argument("--consensus", default="data/raw/lincs_consensi_pert_id.tsv.bz2")
    parser.add_argument("--orthologs", default="data/raw/tage_metadata/Table_of_orthologs.csv")
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--drugs-json", default=None)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.drugs_json:
        with open(args.drugs_json) as f:
            drug_map = json.load(f)
    else:
        drug_map = DEFAULT_DRUGS

    all_pert_ids = [p for pids in drug_map.values() for p in pids]

    # 1. Load clock model
    clock_params = load_clock_model(Path(args.pkl_model))

    # 2. Load module gene lists
    module_genes = load_module_gene_lists(Path(args.supp_table))

    # 3. Load orthologs
    ortholog_map = load_orthologs(Path(args.orthologs))

    # 4. Load LINCS signatures
    sig_df = load_lincs_consensus(Path(args.consensus), all_pert_ids)

    # 5. Align and preprocess
    centered_df = align_and_preprocess(sig_df, ortholog_map, clock_params)

    # 6. Compute fingerprints
    scores_df = compute_fingerprint(centered_df, clock_params, module_genes)

    # 7. Aggregate by drug
    drug_scores = aggregate_by_drug(scores_df, drug_map)

    # 8. Save outputs
    scores_df.to_csv(out_dir / "module_fingerprint_v2_per_pert_id.csv")
    drug_scores.to_csv(out_dir / "module_fingerprint_v2_per_drug.csv")
    print(f"\nSaved per-pert_id scores to {out_dir}/module_fingerprint_v2_per_pert_id.csv")
    print(f"Saved per-drug scores to {out_dir}/module_fingerprint_v2_per_drug.csv")

    # 9. Plot
    plot_fingerprint(drug_scores, out_dir)

    # 10. Print summary
    print("\n" + "=" * 70)
    print("DRUG MODULE FINGERPRINT SUMMARY (v2 — Composite Decomposition)")
    print("=" * 70)

    for drug in drug_scores.index:
        composite = drug_scores.loc[drug, "composite"]
        n_sigs = int(drug_scores.loc[drug, "n_signatures"])
        print(f"\n{drug.upper()} (n={n_sigs} signatures)")
        print(f"  Composite mortality tAge: {composite:+.4f}")
        print(f"  Module contributions:")

        # Sort modules by absolute contribution
        mod_vals = {
            c: drug_scores.loc[drug, c]
            for c in drug_scores.columns
            if c not in ("composite", "non_module", "n_signatures")
        }
        sorted_mods = sorted(mod_vals.items(), key=lambda x: abs(x[1]), reverse=True)

        for mod, val in sorted_mods:
            pct = val / composite * 100 if composite != 0 else 0
            direction = "↑" if val > 0 else "↓" if val < 0 else "→"
            print(f"    {mod:15s}: {val:+.4f} ({pct:+5.1f}%) {direction}")

    print("\n" + "=" * 70)
    print("Pipeline complete.")


if __name__ == "__main__":
    main()
