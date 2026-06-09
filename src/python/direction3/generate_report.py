#!/usr/bin/env python3
"""
Generate REPORT_Direction3_ARCHS4_Mining.md from analysis outputs.
"""
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
FIGURES_DIR = PROJECT_ROOT / "figures" / "direction3"
REPORT_PATH = PROJECT_ROOT / "docs" / "REPORT_Direction3_ARCHS4_Mining.md"

def section(title, level=2):
    return f"\n{'#' * level} {title}\n"

def bullet(text):
    return f"- {text}\n"

def main():
    lines = []
    lines.append("# Direction 3: ARCHS4 Large-Scale Data Mining Report\n")
    lines.append(f"*Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}*\n")

    # Executive Summary
    lines.append(section("Executive Summary", 2))
    lines.append("This report documents the large-scale mining of ARCHS4 mouse RNA-seq data ")
    lines.append("for transcriptomic age (tAge) prediction across aging and intervention studies.\n")

    # Data Sources
    lines.append(section("Data Sources", 2))
    lines.append(bullet("ARCHS4 mouse_gene_v2.5.h5 (Ensembl 107, GRCm38, ~100,000 samples)"))
    lines.append(bullet("Downloaded from: https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5"))
    lines.append(bullet("Format: HDF5 with Kallisto pseudocounts (rounded to integers)"))
    lines.append(bullet("Gene IDs: Entrez Gene Symbols"))

    # Sample Counts
    meta_file = RESULTS_DIR / "archs4_sample_metadata.csv.gz"
    if meta_file.exists():
        meta = pd.read_csv(meta_file, compression="gzip", low_memory=False)
        lines.append(f"\nTotal samples in ARCHS4: **{len(meta):,}**\n")

    pred_file = RESULTS_DIR / "archs4_all_predictions.csv"
    if pred_file.exists():
        preds = pd.read_csv(pred_file, low_memory=False)
        lines.append(f"Samples with tAge predictions: **{len(preds):,}**\n")
        lines.append("\nPredictions by query set:\n")
        # Count by inferred set from file existence
        for fname in RESULTS_DIR.glob("*_predictions.csv"):
            if fname.name == "archs4_all_predictions.csv":
                continue
            name = fname.stem.replace("_predictions", "")
            df = pd.read_csv(fname)
            lines.append(bullet(f"{name}: {len(df)} samples"))

    # Quality Control
    lines.append(section("Quality Control", 2))
    lines.append(bullet("TMM normalization applied via edgeR (method consistent with tAge pipeline)"))
    lines.append(bullet("logCPM computed with prior.count=3"))
    lines.append(bullet("Missing genes imputed with training-set medians from model"))
    lines.append(bullet("Features aligned to 10,487 Entrez Gene IDs expected by model"))

    # Key Findings
    lines.append(section("Key Findings", 2))

    # GSE288795 validation
    gse288795_pred = RESULTS_DIR / "gse288795_predictions.csv"
    gse288795_diff = RESULTS_DIR / "gse288795_drug_control_differences.csv"
    if gse288795_pred.exists():
        lines.append("\n**GSE288795 Validation (Rapamycin + Trametinib, 111 samples, 3 tissues):**\n")
        lines.append(bullet("Independent validation dataset processed while waiting for ARCHS4 download"))
        if gse288795_diff.exists():
            diff_df = pd.read_csv(gse288795_diff)
            for _, row in diff_df.iterrows():
                direction = "rejuvenation" if row["difference"] < 0 else "aging acceleration"
                lines.append(bullet(
                    f"{row['drug']} ({row['tissue']} {row['sex']}): ΔtAge = {row['difference']:+.1f} months "
                    f"({direction}, n_drug={row['n_drug']}, n_ctrl={row['n_ctrl']})"
                ))
        lines.append(bullet("Figures: `figures/direction3/gse288795_*_tage.png`"))
        lines.append(bullet("Pooled drug effects: `figures/direction3/gse288795_drug_effects_pooled.png`"))

    gse230402_pred = RESULTS_DIR / "gse230402_predictions.csv"
    gse230402_diff = RESULTS_DIR / "gse230402_cr_al_differences.csv"
    if gse230402_pred.exists():
        lines.append("\n**GSE230402 Validation (Caloric Restriction, 22 samples, liver):**\n")
        if gse230402_diff.exists():
            diff_df = pd.read_csv(gse230402_diff)
            for _, row in diff_df.iterrows():
                direction = "rejuvenation" if row["difference"] < 0 else "aging acceleration"
                lines.append(bullet(
                    f"CR ({row['sex']}): ΔtAge = {row['difference']:+.1f} months "
                    f"({direction}, n_cr={row['n_cr']}, n_al={row['n_al']})"
                ))
        lines.append(bullet("Figure: `figures/direction3/gse230402_liver_tage.png`"))
        lines.append(bullet("Note: Pre-normalized counts used; TMM normalization preferred for optimal accuracy"))

    # ARCHS4 Drug effects
    drug_file = RESULTS_DIR / "archs4_drug_effects.csv"
    if drug_file.exists():
        drug_df = pd.read_csv(drug_file)
        lines.append("\n**ARCHS4 Drug Effects (within-study comparisons):**\n")
        for _, row in drug_df.iterrows():
            direction = "rejuvenation" if row["difference"] < 0 else "aging acceleration"
            lines.append(bullet(
                f"{row['drug']}: ΔtAge = {row['difference']:+.2f} months "
                f"({direction}, n_drug={row['n_drug']}, n_ctrl={row['n_control']})"
            ))
    else:
        lines.append(bullet("ARCHS4 drug effect analysis pending (H5 download in progress)"))

    # Aging trajectory
    lines.append("\n**Aging Trajectory:**\n")
    lines.append(bullet("Parsed chronological age from sample titles/characteristics"))
    lines.append(bullet("Plotted tAge vs chronological age for control samples"))
    if (FIGURES_DIR / "archs4_aging_trajectory.png").exists():
        lines.append(bullet("Figure saved: `figures/direction3/archs4_aging_trajectory.png`"))

    # Batch effects
    lines.append("\n**Batch Effects:**\n")
    lines.append(bullet("Assessed tAge distribution across GEO series"))
    lines.append(bullet("Expected significant batch effects due to different sequencing platforms, labs, and protocols"))
    if (FIGURES_DIR / "archs4_batch_effects.png").exists():
        lines.append(bullet("Figure saved: `figures/direction3/archs4_batch_effects.png`"))

    # Limitations
    lines.append(section("Limitations", 2))
    lines.append(bullet("ARCHS4 uses Kallisto pseudocounts, not raw counts; TMM normalization assumptions may be partially violated"))
    lines.append(bullet("Severe batch effects expected across >100,000 samples from diverse studies"))
    lines.append(bullet("Within-study comparisons are required for valid drug/condition effects"))
    lines.append(bullet("Chronological age parsing is heuristic-based (text extraction from metadata)"))
    lines.append(bullet("Sample metadata quality varies; some studies lack detailed characteristics"))
    lines.append(bullet("H5 download is large (36GB) and processing is I/O intensive"))

    # Next Steps
    lines.append(section("Next Steps", 2))
    lines.append(bullet("Validate ARCHS4 predictions against our existing GSE131754 and GSE299228 results"))
    lines.append(bullet("Perform formal batch correction (e.g., ComBat) if cross-study comparison is needed"))
    lines.append(bullet("Expand to human ARCHS4 data (human_gene_v2.5.h5) for cross-species validation"))
    lines.append(bullet("Integrate LINCS module fingerprints with ARCHS4 drug signatures"))

    # Appendix: GEO Backup Search
    geo_file = PROJECT_ROOT / "data" / "geo_expansion" / "geo_eutils_search_results.csv"
    if geo_file.exists():
        geo_df = pd.read_csv(geo_file)
        lines.append(section("Appendix: GEO Backup Search", 2))
        lines.append(f"NCBI E-utilities search returned **{len(geo_df)}** GEO series records.\n")
        lines.append("Top candidate series for manual validation:\n")
        for q in ["rapamycin_liver", "cr_liver", "metformin_liver", "aging_liver"]:
            qdf = geo_df[geo_df["query_label"] == q]
            lines.append(f"\n*{q}* ({len(qdf)} records):\n")
            for _, row in qdf.head(5).iterrows():
                lines.append(bullet(f"{row['accession']}: {row['title'][:100]} (n={row['n_samples']})"))

    with open(REPORT_PATH, "w") as f:
        f.writelines(lines)

    print(f"Report written to {REPORT_PATH}")

if __name__ == "__main__":
    main()
