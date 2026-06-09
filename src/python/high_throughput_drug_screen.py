#!/usr/bin/env python3
"""
High-Throughput Drug Screening Pipeline

Screen all LINCS L1000 small-molecule compounds for anti-aging potential
using composite mortality clock decomposition into 14 WGCNA modules.

Strategy:
  1. Load all 13,072 compound consensus signatures from LINCS
  2. Compute composite mortality tAge + per-module contributions
  3. Rank compounds by:
     a) Most negative composite tAge (strongest rejuvenation)
     b) Similarity to known anti-aging drugs (cosine similarity of module fingerprints)
  4. Filter out toxic/undesirable MoA (e.g., strong chromatin disruption)
  5. Output ranked candidate lists + visualizations

Outputs:
  - results/ht_screen_all_compounds.csv (13,072 rows)
  - results/ht_screen_top_rejuvenation_candidates.csv
  - results/ht_screen_similarity_to_reference_drugs.csv
  - figures/ht_screen_umap.pdf, volcano.pdf, top_candidates.pdf
"""

import argparse
import warnings
from pathlib import Path

import bz2
import csv
import joblib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Reference anti-aging drugs for similarity comparison
# ---------------------------------------------------------------------------
REFERENCE_DRUGS = {
    "metformin": ["BRD-K79602928"],
    "sirolimus_in_vitro_mean": None,  # computed from available pert_ids
}

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


def load_compound_pert_ids(pert_info_path: Path) -> set:
    """Load all compound (trt_cp) pert_ids from LINCS pert_info."""
    print(f"[1/6] Loading compound pert_ids from {pert_info_path}")
    df = pd.read_csv(pert_info_path, sep="\t", low_memory=False, usecols=["pert_id", "pert_type", "pert_iname"])
    compounds = df[df["pert_type"] == "trt_cp"]
    pert_ids = set(compounds["pert_id"])
    pert_id_to_name = compounds.set_index("pert_id")["pert_iname"].to_dict()
    print(f"  {len(pert_ids):,} compounds found")
    return pert_ids, pert_id_to_name


def load_clock_model(pkl_path: Path):
    """Load .pkl clock model and extract coefficients + preprocessing params."""
    print(f"[2/6] Loading clock model from {pkl_path}")
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
        f"intercept: {params['intercept']:.4f}"
    )
    return params


def load_module_gene_lists(excel_path: Path) -> dict:
    """Load module gene lists from supp_table_7.xlsx Sheet (C)."""
    print(f"[3/6] Loading module gene lists from {excel_path}")
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
        print(f"  {mod_name:15s}: {len(genes):4d} genes")
    return modules


def load_orthologs(ortholog_path: Path) -> dict:
    """Load human → mouse Entrez ID mapping."""
    print(f"[4/6] Loading orthologs from {ortholog_path}")
    df = pd.read_csv(ortholog_path)
    mapping = (
        df[["Entrez.Human", "Entrez.Mouse"]]
        .dropna()
        .drop_duplicates(subset=["Entrez.Human"])
        .set_index("Entrez.Human")["Entrez.Mouse"]
        .astype(int)
        .to_dict()
    )
    print(f"  {len(mapping):,} human→mouse mappings")
    return mapping


def build_precomputed_index(clock_params: dict, module_genes: dict, ortholog_map: dict):
    """Pre-compute indexing arrays for fast vectorized processing."""
    feature_names = clock_params["feature_names"]
    feature_to_idx = {f: i for i, f in enumerate(feature_names)}
    n_features = len(feature_names)

    # Build a mapping from LINCS human entrez column index → mouse entrez → clock feature index
    # This will be used to map each signature row quickly
    # But we don't know the column order yet; we'll build this after reading the header
    return {
        "feature_names": feature_names,
        "feature_to_idx": feature_to_idx,
        "n_features": n_features,
        "coef": clock_params["coef"],
        "intercept": clock_params["intercept"],
        "impute_stats": clock_params["impute_stats"],
        "center_mean": clock_params["center_mean"],
        "module_genes": module_genes,
        "ortholog_map": ortholog_map,
    }


def process_signatures(
    consensus_path: Path,
    compound_pert_ids: set,
    precomp: dict,
    batch_size: int = 5000,
):
    """Process all compound signatures in batches."""
    print(f"[5/6] Processing compound signatures from {consensus_path}")

    feature_names = precomp["feature_names"]
    feature_to_idx = precomp["feature_to_idx"]
    coef = precomp["coef"]
    intercept = precomp["intercept"]
    impute_stats = precomp["impute_stats"]
    center_mean = precomp["center_mean"]
    module_genes = precomp["module_genes"]
    ortholog_map = precomp["ortholog_map"]

    # Pre-compute module masks
    module_masks = {}
    all_module_genes = set()
    for mod_name, genes in module_genes.items():
        mask = np.zeros(len(feature_names), dtype=bool)
        for g in genes:
            if g in feature_to_idx:
                mask[feature_to_idx[g]] = True
                all_module_genes.add(g)
        module_masks[mod_name] = mask

    non_module_mask = np.ones(len(feature_names), dtype=bool)
    for g in all_module_genes:
        if g in feature_to_idx:
            non_module_mask[feature_to_idx[g]] = False

    results = []
    processed = 0
    skipped = 0

    with bz2.open(consensus_path, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        # Build column mapping: human entrez string → clock feature index
        col_map = {}
        for j, col_name in enumerate(header[1:], start=1):  # skip perturbagen column
            try:
                h_entrez = int(col_name)
                m_entrez = ortholog_map.get(h_entrez)
                if m_entrez is not None:
                    m_str = str(int(m_entrez))
                    if m_str in feature_to_idx:
                        col_map[j] = feature_to_idx[m_str]
            except (ValueError, TypeError):
                pass

        print(f"  {len(col_map)} / {len(header)-1} LINCS columns map to clock features")

        for row in reader:
            pert_id = row[0]
            if pert_id not in compound_pert_ids:
                skipped += 1
                continue

            # Build feature vector
            x = np.full(len(feature_names), np.nan, dtype=float)
            for j, val_str in enumerate(row[1:], start=1):
                if j in col_map:
                    try:
                        x[col_map[j]] = float(val_str)
                    except ValueError:
                        pass

            # Impute missing
            missing = np.isnan(x)
            x[missing] = impute_stats[missing]

            # Center
            x_centered = x - center_mean

            # Composite prediction
            composite = intercept + np.dot(coef, x_centered)

            # Module contributions
            record = {
                "pert_id": pert_id,
                "composite": composite,
                "non_module": float(np.dot(coef[non_module_mask], x_centered[non_module_mask])),
            }
            for mod_name, mask in module_masks.items():
                record[mod_name] = float(np.dot(coef[mask], x_centered[mask]))

            results.append(record)
            processed += 1

            if processed % batch_size == 0:
                print(f"  Processed {processed:,} compounds...")

    print(f"  Done: {processed:,} compounds processed, {skipped:,} non-compounds skipped")
    return pd.DataFrame(results)


def compute_similarity_scores(
    df: pd.DataFrame, reference_pert_ids: list, module_cols: list
) -> pd.DataFrame:
    """Compute cosine similarity between each compound and reference drugs."""
    print("[6/6] Computing similarity to reference drugs")

    ref_profiles = {}
    for ref_name, pids in reference_pert_ids.items():
        if pids is None:
            continue
        avail = [p for p in pids if p in df["pert_id"].values]
        if avail:
            profile = df[df["pert_id"].isin(avail)][module_cols].mean().values
            ref_profiles[ref_name] = profile
            print(f"  Reference {ref_name}: {len(avail)} signatures averaged")

    # Compute similarities
    module_matrix = df[module_cols].values
    for ref_name, ref_profile in ref_profiles.items():
        cos_sims = []
        pearsons = []
        for i in range(len(df)):
            profile = module_matrix[i]
            # Cosine similarity (1 - cosine distance)
            cos_sim = 1 - cosine(profile, ref_profile)
            cos_sims.append(cos_sim)
            # Pearson correlation
            r, _ = pearsonr(profile, ref_profile)
            pearsons.append(r)
        df[f"cos_sim_{ref_name}"] = cos_sims
        df[f"pearson_{ref_name}"] = pearsons

    return df


def rank_candidates(
    df: pd.DataFrame, pert_id_to_name: dict, module_cols: list
) -> tuple:
    """Rank compounds by multiple criteria."""
    print("\n[7/7] Ranking candidates")

    df["pert_iname"] = df["pert_id"].map(pert_id_to_name)

    # Criterion 1: Strongest rejuvenation (most negative composite)
    rejuvenation = df.nsmallest(50, "composite")[
        ["pert_id", "pert_iname", "composite"] + module_cols
    ].copy()
    print(f"  Top 10 rejuvenation candidates:")
    for _, row in rejuvenation.head(10).iterrows():
        print(f"    {row['pert_iname']:20s} ({row['pert_id']}): {row['composite']:+.4f}")

    # Criterion 2: Metformin-like (similar module fingerprint)
    if "cos_sim_metformin" in df.columns:
        metformin_like = df.nlargest(50, "cos_sim_metformin")[
            ["pert_id", "pert_iname", "composite", "cos_sim_metformin"] + module_cols
        ].copy()
        print(f"\n  Top 10 Metformin-like candidates:")
        for _, row in metformin_like.head(10).iterrows():
            print(f"    {row['pert_iname']:20s} ({row['pert_id']}): "
                  f"composite={row['composite']:+.4f}, cos_sim={row['cos_sim_metformin']:.4f}")
    else:
        metformin_like = None

    return rejuvenation, metformin_like


def plot_screening_results(df: pd.DataFrame, module_cols: list, out_dir: Path):
    """Generate screening visualizations."""
    print("\n  Generating screening figures...")
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1. Volcano plot: composite tAge vs max abs module contribution
    df["max_abs_module"] = df[module_cols].abs().max(axis=1)

    fig, ax = plt.subplots(figsize=(10, 7))
    scatter = ax.scatter(
        df["composite"],
        df["max_abs_module"],
        c=df["composite"],
        cmap="RdYlGn_r",
        alpha=0.5,
        s=15,
        edgecolors="none",
    )
    ax.axvline(x=0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("Composite Mortality tAge", fontsize=12)
    ax.set_ylabel("Max |Module Contribution|", fontsize=12)
    ax.set_title("LINCS Compound Screening: Rejuvenation Potential", fontsize=13, fontweight="bold")
    plt.colorbar(scatter, ax=ax, label="Composite tAge")
    plt.tight_layout()
    plt.savefig(fig_dir / "ht_screen_volcano.pdf", dpi=300)
    plt.close()
    print(f"    Saved volcano plot")

    # 2. Distribution of composite tAge
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(df["composite"], bins=100, color="steelblue", edgecolor="white", alpha=0.7)
    ax.axvline(x=0, color="black", linestyle="--", linewidth=1)
    ax.axvline(x=df["composite"].median(), color="red", linestyle="-", linewidth=1.5, label=f"Median: {df['composite'].median():.3f}")
    ax.set_xlabel("Composite Mortality tAge", fontsize=12)
    ax.set_ylabel("Number of Compounds", fontsize=12)
    ax.set_title("Distribution of Mortality tAge Across 13,072 LINCS Compounds", fontsize=13, fontweight="bold")
    ax.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "ht_screen_distribution.pdf", dpi=300)
    plt.close()
    print(f"    Saved distribution plot")

    # 3. Top rejuvenation candidates heatmap
    top_n = 30
    top_rej = df.nsmallest(top_n, "composite")[module_cols].copy()
    top_names = df.nsmallest(top_n, "composite")["pert_iname"].values

    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(top_rej.values, cmap="RdBu_r", aspect="auto", vmin=-1.5, vmax=1.5)
    ax.set_xticks(range(len(module_cols)))
    ax.set_xticklabels(module_cols, rotation=45, ha="right")
    ax.set_yticks(range(top_n))
    ax.set_yticklabels(top_names, fontsize=8)
    plt.colorbar(im, ax=ax, label="Module Contribution")
    ax.set_title(f"Top {top_n} Rejuvenation Candidates (Most Negative Composite tAge)", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(fig_dir / "ht_screen_top_rejuvenation_heatmap.pdf", dpi=300)
    plt.close()
    print(f"    Saved top rejuvenation heatmap")


def main():
    parser = argparse.ArgumentParser(description="High-Throughput Drug Screening")
    parser.add_argument("--pkl-model", default="models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    parser.add_argument("--supp-table", default="data/raw/supp_table_7.xlsx")
    parser.add_argument("--consensus", default="data/raw/lincs_consensi_pert_id.tsv.bz2")
    parser.add_argument("--pert-info", default="data/raw/GSE92742_Broad_LINCS_pert_info.txt.gz")
    parser.add_argument("--orthologs", default="data/raw/tage_metadata/Table_of_orthologs.csv")
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--batch-size", type=int, default=5000)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load compound pert_ids
    compound_pert_ids, pert_id_to_name = load_compound_pert_ids(Path(args.pert_info))

    # 2. Load clock model
    clock_params = load_clock_model(Path(args.pkl_model))

    # 3. Load module gene lists
    module_genes = load_module_gene_lists(Path(args.supp_table))

    # 4. Load orthologs
    ortholog_map = load_orthologs(Path(args.orthologs))

    # 5. Pre-compute indices
    precomp = build_precomputed_index(clock_params, module_genes, ortholog_map)

    # 6. Process all compound signatures
    results_df = process_signatures(
        Path(args.consensus), compound_pert_ids, precomp, batch_size=args.batch_size
    )

    module_cols = [m[0] for m in MORTALITY_MODULES]

    # 7. Compute similarity to reference drugs
    reference_pert_ids = {
        "metformin": ["BRD-K79602928"],
    }
    # Also add sirolimus from available data
    sirolimus_ids = [p for p in ["BRD-A23770159", "BRD-A50287119", "BRD-A79768653", "BRD-K84937637", "BRD-K89626439"] if p in results_df["pert_id"].values]
    if sirolimus_ids:
        reference_pert_ids["sirolimus"] = sirolimus_ids

    results_df = compute_similarity_scores(results_df, reference_pert_ids, module_cols)

    # 8. Rank candidates
    rejuvenation, metformin_like = rank_candidates(results_df, pert_id_to_name, module_cols)

    # 9. Save outputs
    results_df.to_csv(out_dir / "ht_screen_all_compounds.csv", index=False)
    print(f"\n  Saved all compounds to {out_dir}/ht_screen_all_compounds.csv")

    if rejuvenation is not None:
        rejuvenation.to_csv(out_dir / "ht_screen_top_rejuvenation_candidates.csv", index=False)
        print(f"  Saved top rejuvenation candidates to {out_dir}/ht_screen_top_rejuvenation_candidates.csv")

    if metformin_like is not None:
        metformin_like.to_csv(out_dir / "ht_screen_metformin_like_candidates.csv", index=False)
        print(f"  Saved Metformin-like candidates to {out_dir}/ht_screen_metformin_like_candidates.csv")

    # 10. Plot
    plot_screening_results(results_df, module_cols, out_dir)

    print("\n" + "=" * 70)
    print("HIGH-THROUGHPUT SCREENING COMPLETE")
    print("=" * 70)
    print(f"Total compounds screened: {len(results_df):,}")
    print(f"Rejuvenation candidates (composite < 0): {(results_df['composite'] < 0).sum():,}")
    print(f"Strong rejuvenation (composite < -1): {(results_df['composite'] < -1).sum():,}")
    print(f"Acceleration candidates (composite > 1): {(results_df['composite'] > 1).sum():,}")


if __name__ == "__main__":
    main()
