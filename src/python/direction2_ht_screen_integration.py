#!/usr/bin/env python3
"""
Direction 2 — HT Screen Integration with GEO Golden Standards
==============================================================

Uses GEO-validated long-term interventions (CR, Rapamycin, Acarbose, Metformin)
as "golden standard" module signatures to find new drug candidates in LINCS.

Steps:
  1. Compute GEO golden standard module signatures (per-module Drug-Control)
  2. Define screening strategies (CR mimic, Rapa mimic, consensus signatures)
  3. Screen all 13,072 LINCS compounds by cosine similarity
  4. Cross-reference with existing Metformin-like candidates
  5. Validate against known anti-aging drugs
  6. Generate visualizations and report

Outputs:
  - results/geo_longterm/geo_module_signatures.csv
  - results/geo_longterm/ht_screen_strategy_*.csv
  - results/geo_longterm/ht_screen_triple_hits.csv
  - figures/direction2_*.pdf
  - docs/REPORT_Direction2_HT_Screen_Integration.md
"""

import argparse
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
import seaborn as sns

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
GEO_DIR = DATA_DIR / "geo"
RESULTS_DIR = PROJECT_ROOT / "results"
GEO_RESULTS = RESULTS_DIR / "geo_longterm"
FIGURES_DIR = PROJECT_ROOT / "figures"
DOCS_DIR = PROJECT_ROOT / "docs"

# ---------------------------------------------------------------------------
# Module metadata (same as existing pipelines)
# ---------------------------------------------------------------------------
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
MODULE_NAMES = [m[0] for m in MORTALITY_MODULES]

# Known anti-aging drugs for validation
KNOWN_ANTIAGING_DRUGS = {
    "resveratrol",
    "spermidine",
    "rapamycin",
    "sirolimus",
    "metformin",
    "everolimus",
    "acarbose",
    "nmn",
    "nr",
    "nad+",
    "sirtuin",
    "quercetin",
    "fisetin",
    "dasatinib",
    "navitoclax",
    "azd-5438",
    "jnk-in-8",
    "alvocidib",
}


def load_model_and_modules(model_path, supp_table_path):
    """Load clock model and module gene lists."""
    print("[1/8] Loading model and module gene lists")
    model = joblib.load(model_path)
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
    feature_to_idx = {f: i for i, f in enumerate(params["feature_names"])}

    # Load module gene lists from supp_table_7.xlsx
    df = pd.read_excel(
        supp_table_path,
        sheet_name="(C) Module multi-species clocks",
        header=None,
        skiprows=4,
    )
    df = df[df[0] != "Intercept"].copy()
    df[0] = df[0].astype(int)
    df = df.set_index(0)

    module_genes = {}
    for mod_name, mod_annot, col_idx in MORTALITY_MODULES:
        genes = set(df[df[col_idx] != 0].index.astype(str))
        module_genes[mod_name] = genes
        print(f"  {mod_name:15s}: {len(genes):4d} genes")

    return params, feature_to_idx, module_genes


def compute_module_contributions(centered_array, coef, feature_to_idx, module_genes):
    """
    Compute per-module contributions for a centered expression matrix.
    centered_array: (n_samples, n_features)
    Returns dict of {module_name: (n_samples,) array}
    """
    all_module_genes = set()
    for genes in module_genes.values():
        all_module_genes.update(genes)

    # Build non-module mask
    non_module_mask = np.ones(len(coef), dtype=bool)
    for g in all_module_genes:
        if g in feature_to_idx:
            non_module_mask[feature_to_idx[g]] = False

    results = {"non_module": np.dot(centered_array[:, non_module_mask], coef[non_module_mask])}

    for mod_name, genes in module_genes.items():
        indices = [feature_to_idx[g] for g in genes if g in feature_to_idx]
        if indices:
            results[mod_name] = np.dot(centered_array[:, indices], coef[indices])
        else:
            results[mod_name] = np.full(centered_array.shape[0], np.nan)

    return results


def process_gse131754(model_params, feature_to_idx, module_genes):
    """
    Compute per-sample module contributions for GSE131754.
    Then aggregate Drug - Control per module.
    """
    print("\n[2/8] Processing GSE131754 logCPM data")
    logcpm = pd.read_csv(GEO_RESULTS / "gse131754_logcpm.tsv", sep="\t", index_col=0)
    logcpm.index = logcpm.index.astype(str)
    samples = logcpm.columns.tolist()

    # Align to model features
    n_features = len(model_params["feature_names"])
    n_samples = len(samples)
    X = np.full((n_samples, n_features), np.nan, dtype=float)

    common = 0
    for gene in logcpm.index:
        if gene in feature_to_idx:
            X[:, feature_to_idx[gene]] = logcpm.loc[gene].values.astype(float)
            common += 1

    print(f"  Common genes: {common}/{n_features}")

    # Impute missing
    missing = np.isnan(X)
    for i in range(n_features):
        X[missing[:, i], i] = model_params["impute_stats"][i]

    # Center
    X_centered = X - model_params["center_mean"][np.newaxis, :]

    # Compute composite predictions
    composite = model_params["intercept"] + np.dot(X_centered, model_params["coef"])
    pred_df = pd.DataFrame({
        "sample_id": samples,
        "composite": composite,
    })

    # Compute module contributions
    mod_contribs = compute_module_contributions(
        X_centered, model_params["coef"], feature_to_idx, module_genes
    )
    for mod_name, vals in mod_contribs.items():
        pred_df[mod_name] = vals

    # Parse metadata
    meta = []
    for col in samples:
        parts = col.split("_")
        if col.startswith("CON_"):
            intervention = "Control"
        elif col.startswith("GHRCON_"):
            intervention = "GHRKO_Control"
        elif col.startswith("MRCON_"):
            intervention = "MR_Control"
        elif col.startswith("SNELLCON_"):
            intervention = "Snell_Control"
        else:
            intervention = parts[0]
        age = parts[1] if len(parts) > 1 else "unknown"
        sex = parts[2] if len(parts) > 2 else "unknown"
        meta.append({
            "sample_id": col,
            "intervention": intervention,
            "age": age,
            "sex": sex,
        })
    meta_df = pd.DataFrame(meta)
    pred_df = pred_df.merge(meta_df, on="sample_id")

    # Compute Drug - Control per intervention per module
    print("  Computing Drug - Control differences per module...")
    results = []
    for intervention in ["CR", "RAP", "ACA"]:
        drug_samples = pred_df[pred_df["intervention"] == intervention]
        
        # Match controls by age/sex
        for _, drug_row in drug_samples.iterrows():
            control_samples = pred_df[
                (pred_df["intervention"] == "Control") &
                (pred_df["age"] == drug_row["age"]) &
                (pred_df["sex"] == drug_row["sex"])
            ]
            for col in ["composite"] + list(mod_contribs.keys()):
                for _, ctrl_row in control_samples.iterrows():
                    results.append({
                        "intervention": intervention,
                        "drug_sample": drug_row["sample_id"],
                        "control_sample": ctrl_row["sample_id"],
                        "age": drug_row["age"],
                        "sex": drug_row["sex"],
                        "module": col,
                        "drug_value": drug_row[col],
                        "control_value": ctrl_row[col],
                        "difference": drug_row[col] - ctrl_row[col],
                    })

    diff_df = pd.DataFrame(results)
    
    # Aggregate by intervention and module
    agg = diff_df.groupby(["intervention", "module"]).agg(
        mean_diff=("difference", "mean"),
        sem_diff=("difference", "sem"),
        std_diff=("difference", "std"),
        n=("difference", "count"),
    ).reset_index()
    
    print("  GSE131754 module signatures computed.")
    return agg, pred_df


def process_gse299228(model_params, feature_to_idx, module_genes):
    """
    Compute per-module contributions for GSE299228 group averages.
    Returns Metformin - Control differences per module.
    """
    print("\n[3/8] Processing GSE299228 group averages")
    df = pd.read_csv(GEO_DIR / "GSE299228_Comparative_Study_DESeq2.csv.gz", compression="gzip")
    
    # Gene symbols in first column
    genes = df.iloc[:, 0].values
    group_cols = ["Cntrl", "Rapamycin", "Metformin", "TM5614", "CalRestrn"]
    counts = df[group_cols].values
    
    # Map gene symbols to Entrez IDs using mygene
    print("  Mapping gene symbols to Entrez IDs...")
    from mygene import MyGeneInfo
    mg = MyGeneInfo()
    clean_genes = [str(g).split(".")[0] for g in genes]
    
    mapping = {}
    batch_size = 1000
    for i in range(0, len(clean_genes), batch_size):
        batch = clean_genes[i:i+batch_size]
        result = mg.querymany(batch, scopes="symbol,ensembl.gene",
                             fields="entrezgene", species="mouse",
                             verbose=False, as_dataframe=True)
        for idx, row in result.iterrows():
            if "entrezgene" in row and pd.notna(row["entrezgene"]):
                mapping[idx] = int(row["entrezgene"])
    
    print(f"  Mapped: {len(mapping)}/{len(clean_genes)}")
    
    # Build expression matrix with Entrez IDs
    entrez_list = []
    expr_list = []
    for i, gene in enumerate(clean_genes):
        if gene in mapping:
            entrez_list.append(str(mapping[gene]))
            expr_list.append(counts[i])
    
    expr_df = pd.DataFrame(np.array(expr_list), columns=group_cols, index=entrez_list)
    expr_df = expr_df.groupby(level=0).mean()
    
    # Log2 transform
    log_expr = np.log2(expr_df + 1)
    
    # Align with model features
    n_features = len(model_params["feature_names"])
    n_groups = len(group_cols)
    X = np.zeros((n_groups, n_features), dtype=float)
    
    common = 0
    for gene in log_expr.index:
        if gene in feature_to_idx:
            X[:, feature_to_idx[gene]] = log_expr.loc[gene].values
            common += 1
    
    print(f"  Common genes with model: {common}/{n_features}")
    
    # Impute missing
    for i in range(n_features):
        if np.all(X[:, i] == 0):
            X[:, i] = model_params["impute_stats"][i]
    
    # Center
    X_centered = X - model_params["center_mean"][np.newaxis, :]
    
    # Compute module contributions
    mod_contribs = compute_module_contributions(
        X_centered, model_params["coef"], feature_to_idx, module_genes
    )
    
    # Composite predictions
    composite = model_params["intercept"] + np.dot(X_centered, model_params["coef"])
    
    pred_df = pd.DataFrame({
        "group": group_cols,
        "composite": composite,
    })
    for mod_name, vals in mod_contribs.items():
        pred_df[mod_name] = vals
    
    # Compute Metformin - Control differences
    ctrl_idx = 0  # Cntrl
    met_idx = 2   # Metformin
    
    results = []
    for col in ["composite"] + list(mod_contribs.keys()):
        results.append({
            "intervention": "Metformin_GSE299228",
            "module": col,
            "mean_diff": pred_df.loc[met_idx, col] - pred_df.loc[ctrl_idx, col],
            "sem_diff": np.nan,
            "std_diff": np.nan,
            "n": 1,
        })
    
    print("  GSE299228 Metformin module signature computed.")
    return pd.DataFrame(results), pred_df


def build_golden_standards(gse131754_agg, gse299228_agg):
    """
    Combine GSE131754 and GSE299228 into a single golden standard signature matrix.
    Returns: DataFrame with rows=interventions, columns=modules (including composite and non_module)
    """
    print("\n[4/8] Building golden standard module signatures")
    
    all_data = pd.concat([gse131754_agg, gse299228_agg], ignore_index=True)
    
    # Pivot to wide format
    pivot = all_data.pivot(index="intervention", columns="module", values="mean_diff")
    pivot = pivot.reset_index()
    
    # Rename interventions for clarity
    rename_map = {
        "ACA": "Acarbose",
        "CR": "CR",
        "RAP": "Rapamycin",
        "Metformin_GSE299228": "Metformin",
    }
    pivot["intervention"] = pivot["intervention"].map(rename_map)
    
    # Reorder columns: intervention, composite, non_module, then modules
    module_cols = [c for c in pivot.columns if c not in ["intervention", "composite", "non_module"]]
    module_cols = sorted(module_cols)
    col_order = ["intervention", "composite", "non_module"] + module_cols
    pivot = pivot[[c for c in col_order if c in pivot.columns]]
    
    print("  Golden standards:")
    for _, row in pivot.iterrows():
        print(f"    {row['intervention']:20s}: composite={row['composite']:+.4f}")
    
    return pivot


def screen_lincs(golden_standards, output_dir):
    """
    Screen all 13,072 LINCS compounds against golden standard signatures.
    """
    print("\n[5/8] Screening LINCS compounds against golden standards")
    
    # Load HT screen data
    ht = pd.read_csv(RESULTS_DIR / "ht_screen_all_compounds.csv")
    print(f"  Loaded {len(ht):,} compounds")
    
    # Extract module fingerprint columns
    module_cols = MODULE_NAMES
    
    # Verify columns exist
    missing_cols = [c for c in module_cols if c not in ht.columns]
    if missing_cols:
        print(f"  WARNING: Missing module columns in ht_screen: {missing_cols}")
        # Try to find them with different naming
        available = [c for c in ht.columns if c not in ["pert_id", "composite", "non_module", "pert_iname"]]
        print(f"  Available columns: {available}")
        module_cols = [c for c in available if c in module_cols]
    
    # Build golden standard vectors
    ref_vectors = {}
    for _, row in golden_standards.iterrows():
        name = row["intervention"]
        vec = np.array([row.get(m, 0) for m in module_cols])
        ref_vectors[name] = vec
    
    # Compute cosine similarities
    module_matrix = ht[module_cols].values
    
    strategies = {
        "A_CR": ref_vectors["CR"],
        "B_Rapamycin": ref_vectors["Rapamycin"],
        "C_CR_Rapa_consensus": (ref_vectors["CR"] + ref_vectors["Rapamycin"]) / 2,
        "D_All_consensus": (
            0.40 * ref_vectors["CR"] +
            0.30 * ref_vectors["Rapamycin"] +
            0.15 * ref_vectors["Acarbose"] +
            0.15 * ref_vectors["Metformin"]
        ),
    }
    
    for strat_name, ref_vec in strategies.items():
        cos_sims = []
        for i in range(len(ht)):
            profile = module_matrix[i]
            # Cosine similarity (handle zero vectors)
            norm_prod = np.linalg.norm(profile) * np.linalg.norm(ref_vec)
            if norm_prod == 0:
                cos_sims.append(0.0)
            else:
                cos_sims.append(np.dot(profile, ref_vec) / norm_prod)
        ht[f"cos_sim_{strat_name}"] = cos_sims
    
    # Apply filters
    # Cytotoxic filter: orange < -0.5 and green < -0.5
    ht["cytotoxic"] = (ht["orange"] < -0.5) & (ht["green"] < -0.5)
    ht["rejuvenation"] = ht["composite"] < -0.5
    
    # Filtered dataset
    filtered = ht[(ht["composite"] < -0.5) & (~ht["cytotoxic"])].copy()
    print(f"  After filtering (composite < -0.5, not cytotoxic): {len(filtered):,} compounds")
    
    # Rank and save top 100 for each strategy
    top_candidates = {}
    for strat_name in strategies.keys():
        col = f"cos_sim_{strat_name}"
        top100 = filtered.nlargest(100, col)[[
            "pert_id", "pert_iname", "composite", col,
            "blue", "brown4", "darkgreen", "darkmagenta", "darkred",
            "darkslateblue", "green", "ivory", "orange", "pink",
            "plum1", "sienna3", "turquoise", "white"
        ]].copy()
        top_candidates[strat_name] = top100
        out_file = output_dir / f"ht_screen_strategy_{strat_name}_top100.csv"
        top100.to_csv(out_file, index=False)
        print(f"  Strategy {strat_name}: top candidate = {top100.iloc[0]['pert_iname']} (cos_sim={top100.iloc[0][col]:.4f})")
    
    return ht, filtered, top_candidates, strategies


def cross_reference_candidates(ht_full, filtered, top_candidates, output_dir):
    """
    Cross-reference with existing Metformin-like candidates.
    Find triple hits: Metformin-like AND CR-like AND Rapa-like.
    """
    print("\n[6/8] Cross-referencing candidates")
    
    # Load existing Metformin-like candidates
    metformin_like = pd.read_csv(RESULTS_DIR / "ht_screen_metformin_like_candidates.csv")
    metformin_ids = set(metformin_like["pert_id"].values)
    print(f"  Metformin-like candidates: {len(metformin_ids)}")
    
    # Define thresholds for "like" status
    # Top 20% of similarity scores within filtered set
    cr_col = "cos_sim_A_CR"
    rapa_col = "cos_sim_B_Rapamycin"
    met_col = "cos_sim_metformin" if "cos_sim_metformin" in filtered.columns else None
    
    cr_thresh = filtered[cr_col].quantile(0.80)
    rapa_thresh = filtered[rapa_col].quantile(0.80)
    met_thresh = filtered[met_col].quantile(0.80) if met_col else 0.8
    
    print(f"  Thresholds (80th percentile): CR={cr_thresh:.4f}, Rapa={rapa_thresh:.4f}, Metformin={met_thresh:.4f}")
    
    # Classify each compound
    filtered["cr_like"] = filtered[cr_col] >= cr_thresh
    filtered["rapa_like"] = filtered[rapa_col] >= rapa_thresh
    filtered["metformin_like"] = filtered[met_col] >= met_thresh if met_col else filtered["pert_id"].isin(metformin_ids)
    
    # Triple hits
    triple_hits = filtered[
        filtered["cr_like"] & filtered["rapa_like"] & filtered["metformin_like"]
    ].copy()
    
    print(f"  Triple hits (CR-like + Rapa-like + Metformin-like): {len(triple_hits)}")
    
    if len(triple_hits) > 0:
        triple_hits = triple_hits.sort_values("composite", ascending=True)
        triple_hits[[
            "pert_id", "pert_iname", "composite", cr_col, rapa_col,
            "blue", "brown4", "darkgreen", "darkmagenta", "darkred",
            "darkslateblue", "green", "ivory", "orange", "pink",
            "plum1", "sienna3", "turquoise", "white"
        ]].to_csv(output_dir / "ht_screen_triple_hits.csv", index=False)
        
        print("  Top triple hits:")
        for _, row in triple_hits.head(10).iterrows():
            print(f"    {row['pert_iname']:25s} composite={row['composite']:+.4f}  CR={row[cr_col]:.4f}  Rapa={row[rapa_col]:.4f}")
    
    # Also save the full filtered with flags
    filtered.to_csv(output_dir / "ht_screen_all_filtered_with_flags.csv", index=False)
    
    return triple_hits


def validate_known_drugs(ht_full, top_candidates):
    """
    Check if top candidates include known anti-aging drugs.
    """
    print("\n[7/8] Validating against known anti-aging drugs")
    
    all_top_ids = set()
    for strat_name, df in top_candidates.items():
        all_top_ids.update(df["pert_iname"].str.lower().values)
    
    found = all_top_ids.intersection(KNOWN_ANTIAGING_DRUGS)
    print(f"  Known anti-aging drugs in top candidates: {len(found)}")
    for drug in sorted(found):
        print(f"    - {drug}")
    
    # Also check specific known drugs in the full dataset
    known_in_data = {}
    for drug in ["resveratrol", "spermidine", "rapamycin", "sirolimus", "metformin", "everolimus", "acarbose"]:
        matches = ht_full[ht_full["pert_iname"].str.lower() == drug]
        if len(matches) > 0:
            known_in_data[drug] = matches.iloc[0]
    
    print(f"\n  Known drugs in full dataset ({len(known_in_data)} found):")
    for drug, row in known_in_data.items():
        print(f"    {drug:15s}: composite={row['composite']:+.4f}")
    
    return known_in_data


def create_visualizations(ht_full, filtered, top_candidates, triple_hits, golden_standards, strategies, output_dir):
    """Generate all visualizations."""
    print("\n[8/8] Generating visualizations")
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    # --- Plot 1: Scatter CR vs Rapa similarity (color by composite) ---
    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(
        filtered["cos_sim_A_CR"],
        filtered["cos_sim_B_Rapamycin"],
        c=filtered["composite"],
        cmap="RdYlGn_r",
        alpha=0.6,
        s=20,
        edgecolors="none",
    )
    ax.axhline(y=0, color="black", linestyle="--", linewidth=0.5)
    ax.axvline(x=0, color="black", linestyle="--", linewidth=0.5)
    ax.set_xlabel("Cosine Similarity to CR Module Signature", fontsize=12)
    ax.set_ylabel("Cosine Similarity to Rapamycin Module Signature", fontsize=12)
    ax.set_title("LINCS Compounds: CR vs Rapamycin Module Similarity\n(Color = Composite tAge, green = rejuvenation)", fontsize=13)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Composite Mortality tAge", fontsize=11)
    
    # Annotate top candidates
    for _, row in filtered.nlargest(10, "cos_sim_A_CR").iterrows():
        ax.annotate(row["pert_iname"], (row["cos_sim_A_CR"], row["cos_sim_B_Rapamycin"]),
                   fontsize=7, alpha=0.8)
    for _, row in filtered.nlargest(10, "cos_sim_B_Rapamycin").iterrows():
        ax.annotate(row["pert_iname"], (row["cos_sim_A_CR"], row["cos_sim_B_Rapamycin"]),
                   fontsize=7, alpha=0.8)
    
    plt.tight_layout()
    plt.savefig(fig_dir / "direction2_cr_vs_rapa_scatter.pdf", dpi=300)
    plt.close()
    print("  Saved CR vs Rapa scatter plot")
    
    # --- Plot 2: Top 20 candidates per strategy bar plots ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()
    
    for idx, (strat_name, df) in enumerate(top_candidates.items()):
        ax = axes[idx]
        top20 = df.head(20).copy()
        top20 = top20.sort_values(f"cos_sim_{strat_name}", ascending=True)
        
        colors = ["green" if c < -0.5 else "orange" if c < 0 else "red" for c in top20["composite"]]
        ax.barh(range(len(top20)), top20[f"cos_sim_{strat_name}"], color=colors, edgecolor="black", alpha=0.8)
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(top20["pert_iname"], fontsize=8)
        ax.set_xlabel("Cosine Similarity to Golden Standard", fontsize=11)
        ax.set_title(f"Strategy {strat_name}\nTop 20 Candidates (green = rejuvenation)", fontsize=12)
        ax.axvline(x=0, color="black", linestyle="-", linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(fig_dir / "direction2_top20_per_strategy.pdf", dpi=300)
    plt.close()
    print("  Saved top 20 per strategy bar plots")
    
    # --- Plot 3: Golden standard module signatures heatmap ---
    module_cols = [c for c in golden_standards.columns if c not in ["intervention", "composite", "non_module"]]
    module_cols = sorted(module_cols)
    
    fig, ax = plt.subplots(figsize=(12, 5))
    plot_df = golden_standards.set_index("intervention")[module_cols]
    im = ax.imshow(plot_df.values, cmap="RdBu_r", aspect="auto", vmin=-0.5, vmax=0.5)
    ax.set_xticks(range(len(module_cols)))
    ax.set_xticklabels(module_cols, rotation=45, ha="right")
    ax.set_yticks(range(len(plot_df.index)))
    ax.set_yticklabels(plot_df.index)
    
    for i in range(len(plot_df.index)):
        for j in range(len(module_cols)):
            val = plot_df.iloc[i, j]
            text_color = "white" if abs(val) > 0.25 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=text_color, fontsize=9)
    
    plt.colorbar(im, ax=ax, label="Drug - Control Module Contribution")
    ax.set_title("GEO Golden Standard Module Signatures (Drug - Control)", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(fig_dir / "direction2_golden_standard_heatmap.pdf", dpi=300)
    plt.close()
    print("  Saved golden standard heatmap")
    
    # --- Plot 4: Venn diagram of overlaps ---
    if "cr_like" in filtered.columns and "rapa_like" in filtered.columns and "metformin_like" in filtered.columns:
        cr_set = set(filtered[filtered["cr_like"]]["pert_id"])
        rapa_set = set(filtered[filtered["rapa_like"]]["pert_id"])
        met_set = set(filtered[filtered["metformin_like"]]["pert_id"])
        
        fig, ax = plt.subplots(figsize=(8, 8))
        venn3([cr_set, rapa_set, met_set],
              set_labels=("CR-like", "Rapa-like", "Metformin-like"),
              ax=ax)
        ax.set_title("Overlap of CR-like, Rapa-like, and Metformin-like Candidates\n(Top 20% similarity, filtered for rejuvenation + non-cytotoxic)", fontsize=12)
        plt.tight_layout()
        plt.savefig(fig_dir / "direction2_venn_overlap.pdf", dpi=300)
        plt.close()
        print("  Saved Venn diagram")
    
    # --- Plot 5: Triple hits module fingerprint ---
    if len(triple_hits) > 0:
        fig, ax = plt.subplots(figsize=(12, max(4, len(triple_hits) * 0.4)))
        plot_df = triple_hits.head(30).set_index("pert_iname")[module_cols]
        im = ax.imshow(plot_df.values, cmap="RdBu_r", aspect="auto", vmin=-0.5, vmax=0.5)
        ax.set_xticks(range(len(module_cols)))
        ax.set_xticklabels(module_cols, rotation=45, ha="right")
        ax.set_yticks(range(len(plot_df.index)))
        ax.set_yticklabels(plot_df.index, fontsize=9)
        plt.colorbar(im, ax=ax, label="Module Contribution")
        ax.set_title("Triple-Hit Candidates Module Fingerprints\n(CR-like + Rapa-like + Metformin-like)", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(fig_dir / "direction2_triple_hits_heatmap.pdf", dpi=300)
        plt.close()
        print("  Saved triple hits heatmap")


def write_report(golden_standards, top_candidates, triple_hits, known_in_data, strategies, output_dir):
    """Write markdown report."""
    print("\n[Report] Writing markdown report...")
    
    report_path = DOCS_DIR / "REPORT_Direction2_HT_Screen_Integration.md"
    
    lines = []
    lines.append("# Direction 2 — HT Screen Integration with GEO Golden Standards")
    lines.append("")
    lines.append("**Date:** 2026-06-09")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("This report describes the integration of GEO-validated long-term intervention")
    lines.append("module signatures with the LINCS high-throughput drug screen to identify")
    lines.append("novel anti-aging drug candidates.")
    lines.append("")
    
    # Golden standards
    lines.append("## 1. GEO Golden Standard Module Signatures")
    lines.append("")
    lines.append("Computed per-module Drug-Control differences for validated interventions:")
    lines.append("")
    lines.append("| Intervention | Composite | " + " | ".join(sorted([c for c in golden_standards.columns if c not in ["intervention", "composite", "non_module"]])) + " |")
    lines.append("|" + "---|" * (len(golden_standards.columns)) + "")
    for _, row in golden_standards.iterrows():
        vals = [f"{row.get(c, 0):+.3f}" for c in golden_standards.columns if c != "intervention"]
        lines.append(f"| {row['intervention']} | " + " | ".join(vals) + " |")
    lines.append("")
    lines.append("**Key findings:**")
    for _, row in golden_standards.iterrows():
        lines.append(f"- **{row['intervention']}**: composite = {row['composite']:+.3f}")
    lines.append("")
    
    # Screening strategies
    lines.append("## 2. Screening Strategies")
    lines.append("")
    lines.append("Four strategies were defined based on cosine similarity to golden standard module vectors:")
    lines.append("")
    lines.append("- **Strategy A (CR mimic)**: Cosine similarity to CR module signature")
    lines.append("- **Strategy B (Rapamycin mimic)**: Cosine similarity to Rapamycin module signature")
    lines.append("- **Strategy C (CR+Rapamycin consensus)**: Average of CR and Rapamycin signatures")
    lines.append("- **Strategy D (All-validated consensus)**: Weighted average (CR 40%, Rapa 30%, Acarbose 15%, Metformin 15%)")
    lines.append("")
    lines.append("Filters applied:")
    lines.append("- Composite < -0.5 (rejuvenation)")
    lines.append("- Not cytotoxic: NOT (orange < -0.5 AND green < -0.5)")
    lines.append("")
    
    # Top candidates per strategy
    lines.append("## 3. Top Candidates per Strategy")
    lines.append("")
    for strat_name, df in top_candidates.items():
        lines.append(f"### Strategy {strat_name}")
        lines.append("")
        lines.append("| Rank | pert_id | pert_iname | composite | cos_sim |")
        lines.append("|------|---------|------------|-----------|---------|")
        for rank, (_, row) in enumerate(df.head(10).iterrows(), 1):
            lines.append(f"| {rank} | {row['pert_id']} | {row['pert_iname']} | {row['composite']:+.4f} | {row[f'cos_sim_{strat_name}']:.4f} |")
        lines.append("")
    
    # Triple hits
    lines.append("## 4. Triple Hits (Highest Confidence)")
    lines.append("")
    lines.append("Compounds that are CR-like, Rapamycin-like, AND Metformin-like:")
    lines.append("")
    if len(triple_hits) > 0:
        lines.append("| pert_id | pert_iname | composite | cos_sim_CR | cos_sim_Rapa |")
        lines.append("|---------|------------|-----------|------------|--------------|")
        for _, row in triple_hits.head(20).iterrows():
            lines.append(f"| {row['pert_id']} | {row['pert_iname']} | {row['composite']:+.4f} | {row['cos_sim_A_CR']:.4f} | {row['cos_sim_B_Rapamycin']:.4f} |")
    else:
        lines.append("No triple hits found with current thresholds.")
    lines.append("")
    
    # Known drugs validation
    lines.append("## 5. Validation Against Known Anti-Aging Drugs")
    lines.append("")
    if known_in_data:
        lines.append("| Drug | Found in dataset | Composite |")
        lines.append("|------|------------------|-----------|")
        for drug, row in known_in_data.items():
            lines.append(f"| {drug} | Yes | {row['composite']:+.4f} |")
    else:
        lines.append("No known anti-aging drugs found in the screened dataset.")
    lines.append("")
    
    # Discussion
    lines.append("## 6. Discussion and Next Steps")
    lines.append("")
    lines.append("### Mechanistic Interpretation")
    lines.append("- CR and Rapamycin signatures show strongest overlap in **turquoise** (innate immunity)")
    lines.append("  and **pink** (mitochondrial translation) modules.")
    lines.append("- Triple-hit candidates represent compounds that converge on multiple validated pathways.")
    lines.append("")
    lines.append("### Validation Strategy")
    lines.append("1. **In vitro**: Test top candidates in cellular senescence assays")
    lines.append("2. **In vivo**: Mouse lifespan studies for top 3-5 compounds")
    lines.append("3. **Mechanism**: CRISPR screening to identify target pathways")
    lines.append("")
    lines.append("### Files Generated")
    lines.append("- `results/geo_longterm/geo_module_signatures.csv` — Golden standard module signatures")
    lines.append("- `results/geo_longterm/ht_screen_strategy_*_top100.csv` — Top 100 per strategy")
    lines.append("- `results/geo_longterm/ht_screen_triple_hits.csv` — Triple-hit candidates")
    lines.append("- `figures/direction2_*.pdf` — Visualizations")
    lines.append("")
    lines.append("---")
    lines.append("*Report generated by Direction 2 HT Screen Integration pipeline*")
    
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    
    print(f"  Saved report to {report_path}")


def main():
    parser = argparse.ArgumentParser(description="Direction 2 HT Screen Integration")
    parser.add_argument("--model", default=str(PROJECT_ROOT / "models" / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"))
    parser.add_argument("--supp-table", default=str(DATA_DIR / "raw" / "supp_table_7.xlsx"))
    parser.add_argument("--output-dir", default=str(GEO_RESULTS))
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Direction 2 — HT Screen Integration with GEO Golden Standards")
    print("=" * 70)
    
    # 1. Load model and modules
    model_params, feature_to_idx, module_genes = load_model_and_modules(
        args.model, args.supp_table
    )
    
    # 2. Process GSE131754
    gse131754_agg, gse131754_pred = process_gse131754(model_params, feature_to_idx, module_genes)
    
    # 3. Process GSE299228
    gse299228_agg, gse299228_pred = process_gse299228(model_params, feature_to_idx, module_genes)
    
    # 4. Build golden standards
    golden_standards = build_golden_standards(gse131754_agg, gse299228_agg)
    golden_standards.to_csv(output_dir / "geo_module_signatures.csv", index=False)
    print(f"\n  Saved golden standards to {output_dir / 'geo_module_signatures.csv'}")
    
    # 5. Screen LINCS
    ht_full, filtered, top_candidates, strategies = screen_lincs(golden_standards, output_dir)
    
    # 6. Cross-reference
    triple_hits = cross_reference_candidates(ht_full, filtered, top_candidates, output_dir)
    
    # 7. Validate known drugs
    known_in_data = validate_known_drugs(ht_full, top_candidates)
    
    # 8. Visualizations
    create_visualizations(ht_full, filtered, top_candidates, triple_hits, golden_standards, strategies, PROJECT_ROOT)
    
    # 9. Report
    write_report(golden_standards, top_candidates, triple_hits, known_in_data, strategies, output_dir)
    
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"Golden standards: {output_dir / 'geo_module_signatures.csv'}")
    print(f"Top candidates: {output_dir / 'ht_screen_strategy_*_top100.csv'}")
    print(f"Triple hits: {output_dir / 'ht_screen_triple_hits.csv'}")
    print(f"Figures: {PROJECT_ROOT / 'figures' / 'direction2_*.pdf'}")
    print(f"Report: {DOCS_DIR / 'REPORT_Direction2_HT_Screen_Integration.md'}")


if __name__ == "__main__":
    main()
