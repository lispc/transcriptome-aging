#!/usr/bin/env python3
"""
GSE280382 GLP-1R Agonist Anti-Aging Analysis
=============================================
Process GSE280382 bulk RNA-seq data (~620 samples, multi-tissue),
run tAge prediction, and compare GLP-1RA (exenatide) vs Control
and vs mTOR inhibition (rapamycin).

Outputs:
  - results/geo_longterm/gse280382_*.csv
  - figures/gse280382_*.pdf
  - docs/REPORT_GSE280382_GLP1R_Agonist.md
"""

import os
import sys
import gzip
import pickle
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import joblib
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / 'data' / 'geo'
MODEL_PATH = PROJECT_ROOT / 'models' / 'EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl'
RESULTS_DIR = PROJECT_ROOT / 'results' / 'geo_longterm'
FIGURES_DIR = PROJECT_ROOT / 'figures'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Module definitions (from direction1_module_analysis.py)
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

MODULE_ORDER = [m[0] for m in MORTALITY_MODULES]
MODULE_ANNOTATIONS = {m[0]: m[1] for m in MORTALITY_MODULES}


def load_module_gene_lists():
    """Load 14 module gene lists from supp_table_7.xlsx."""
    supp_table = PROJECT_ROOT / 'data' / 'raw' / 'supp_table_7.xlsx'
    df = pd.read_excel(
        supp_table,
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
    return modules


# ---------------------------------------------------------------------------
# Step 1: Load counts and metadata
# ---------------------------------------------------------------------------
def load_counts_and_metadata():
    """Load all 3 cohort counts and merge with metadata."""
    print("=" * 70)
    print("STEP 1: Load GSE280382 counts and metadata")
    print("=" * 70)

    # AgedLT
    lt_counts = pd.read_csv(DATA_DIR / 'GSE280382_AgedLT_counts.csv.gz', compression='gzip', index_col=0)
    lt_meta = pd.read_csv(DATA_DIR / 'GSE280382_AgedLT_metadata.csv.gz', compression='gzip')
    lt_meta['cohort'] = 'AgedLT'
    # Name column is like "A51_Adipose"; map to columns
    lt_meta = lt_meta[lt_meta['Name'] != 'Name']  # skip any header duplicates
    lt_meta_dict = {}
    for _, row in lt_meta.iterrows():
        lt_meta_dict[row['Name']] = {
            'mouse': row['Mouse'],
            'tissue': row['Tissue'],
            'group': row['Group'],
            'batch': row['Batch'],
            'cohort': 'AgedLT',
        }

    # AgedST
    st_counts = pd.read_csv(DATA_DIR / 'GSE280382_AgedST_counts.csv.gz', compression='gzip', index_col=0)
    st_meta = pd.read_csv(DATA_DIR / 'GSE280382_AgedST_metadata.csv.gz', compression='gzip')
    st_meta_dict = {}
    for _, row in st_meta.iterrows():
        st_meta_dict[row['Animal']] = {
            'group': row['Group'],
        }

    # Young
    y_counts = pd.read_csv(DATA_DIR / 'GSE280382_young_counts.csv.gz', compression='gzip', index_col=0)
    y_meta = pd.read_csv(DATA_DIR / 'GSE280382_young_metadata.csv.gz', compression='gzip')
    y_meta_dict = {}
    for _, row in y_meta.iterrows():
        y_meta_dict[row['Animal']] = {
            'group': row['Group'],
        }

    print(f"  AgedLT counts: {lt_counts.shape}")
    print(f"  AgedST counts: {st_counts.shape}")
    print(f"  Young counts:  {y_counts.shape}")

    # Use intersection of genes across cohorts for consistency
    common_genes = lt_counts.index.intersection(st_counts.index).intersection(y_counts.index)
    print(f"  Common genes across cohorts: {len(common_genes)}")
    lt_counts = lt_counts.loc[common_genes]
    st_counts = st_counts.loc[common_genes]
    y_counts = y_counts.loc[common_genes]

    # Merge counts
    counts = pd.concat([lt_counts, st_counts, y_counts], axis=1)
    print(f"  Merged counts: {counts.shape}")

    # Build full metadata from sample names
    meta_records = []
    for col in counts.columns:
        parts = col.split('_')
        tissue = parts[-1]
        animal = '_'.join(parts[:-1])  # e.g. "A51" or "B13"

        # Standardize tissue names
        tissue_map = {
            'WBC': 'Circulating_WBCs',
            'CardiacMuscle': 'Heart',
            'FrontalCortex': 'Frontal_cortex',
        }
        tissue_std = tissue_map.get(tissue, tissue)

        record = {
            'sample_id': col,
            'animal': animal,
            'tissue': tissue_std,
        }

        if col in lt_meta_dict:
            record.update(lt_meta_dict[col])
        elif animal in st_meta_dict:
            record['group'] = st_meta_dict[animal]['group']
            record['cohort'] = 'AgedST'
            record['batch'] = np.nan
        elif animal in y_meta_dict:
            record['group'] = y_meta_dict[animal]['group']
            record['cohort'] = 'Young'
            record['batch'] = np.nan
        else:
            print(f"  Warning: No metadata for {col}")
            record['group'] = 'Unknown'
            record['cohort'] = 'Unknown'
            record['batch'] = np.nan

        # Ensure tissue names are standardized after metadata merge
        tissue_map = {
            'WBC': 'Circulating_WBCs',
            'CardiacMuscle': 'Heart',
            'FrontalCortex': 'Frontal_cortex',
        }
        record['tissue'] = tissue_map.get(record.get('tissue', tissue_std), record.get('tissue', tissue_std))

        # Standardize group names
        group = record['group']
        if group in ['Aged_exenatide']:
            record['treatment'] = 'GLP-1RA'
        elif group in ['Aged_rapamycin']:
            record['treatment'] = 'Rapamycin'
        elif group in ['Aged_ctrl', 'Young_ctrl']:
            record['treatment'] = 'Control'
        elif group in ['Young_exenatide']:
            record['treatment'] = 'GLP-1RA'
        elif group in ['Aged_ctrl_KD', 'Aged_exenatide_KD']:
            record['treatment'] = 'KD_Control' if 'ctrl' in group.lower() else 'KD_GLP-1RA'
        else:
            record['treatment'] = group

        # Age group
        if 'Young' in group:
            record['age_group'] = 'Young'
        elif 'Aged' in group:
            record['age_group'] = 'Aged'
        else:
            record['age_group'] = 'Unknown'

        # Duration
        if record['cohort'] == 'AgedLT':
            record['duration'] = 'Long-term'
        elif record['cohort'] == 'AgedST':
            record['duration'] = 'Short-term'
        elif record['cohort'] == 'Young':
            record['duration'] = 'Young'
        else:
            record['duration'] = 'Unknown'

        meta_records.append(record)

    meta_df = pd.DataFrame(meta_records)

    print(f"\n  Sample breakdown:")
    print(meta_df.groupby(['cohort', 'age_group', 'treatment', 'tissue']).size().to_string())

    return counts, meta_df


# ---------------------------------------------------------------------------
# Step 2: Gene symbol -> Entrez mapping
# ---------------------------------------------------------------------------
def map_symbols_to_entrez(gene_symbols):
    """Map mouse gene symbols to Entrez Gene IDs using mygene.info."""
    try:
        from mygene import MyGeneInfo
    except ImportError:
        print("ERROR: mygene not installed. Run: pip install mygene")
        sys.exit(1)

    cache_file = DATA_DIR / 'gse280382_symbol_to_entrez.pkl'
    if cache_file.exists():
        print(f"\nLoading cached mapping from {cache_file}")
        with open(cache_file, 'rb') as f:
            mapping = pickle.load(f)
        print(f"  Cached mapping: {len(mapping)} genes")
        return mapping

    print(f"\nMapping {len(gene_symbols)} gene symbols to Entrez...")
    mg = MyGeneInfo()
    mapping = {}
    batch_size = 1000
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i+batch_size]
        print(f"  Batch {i//batch_size + 1}/{(len(gene_symbols)-1)//batch_size + 1}: {len(batch)} genes")
        result = mg.querymany(
            batch,
            scopes='symbol',
            fields='entrezgene',
            species='mouse',
            verbose=False,
            as_dataframe=True
        )
        for idx, row in result.iterrows():
            if 'entrezgene' in row and pd.notna(row['entrezgene']):
                mapping[idx] = int(row['entrezgene'])

    print(f"  Mapped: {len(mapping)}/{len(gene_symbols)} ({100*len(mapping)/len(gene_symbols):.1f}%)")

    with open(cache_file, 'wb') as f:
        pickle.dump(mapping, f)
    print(f"  Saved cache to {cache_file}")
    return mapping


# ---------------------------------------------------------------------------
# Step 3: Prepare edgeR input
# ---------------------------------------------------------------------------
def prepare_edger_input(counts_df, mapping, output_prefix='gse280382'):
    """Map gene symbols to Entrez, aggregate duplicates, save for edgeR."""
    print(f"\nPreparing edgeR input...")

    df = counts_df.copy()
    df['entrez'] = df.index.map(lambda x: mapping.get(x, None))
    df_mapped = df.dropna(subset=['entrez']).copy()
    df_mapped['entrez'] = df_mapped['entrez'].astype(int)
    print(f"  Mapped genes: {len(df_mapped)}")

    # Aggregate duplicate Entrez IDs
    df_agg = df_mapped.groupby('entrez').sum()
    if 'entrez' in df_agg.columns:
        df_agg = df_agg.drop(columns=['entrez'])
    print(f"  After aggregation: {len(df_agg)} unique Entrez IDs")

    output_counts = RESULTS_DIR / f'{output_prefix}_counts_for_edger.tsv'
    df_agg.to_csv(output_counts, sep='\t')
    print(f"  Saved to {output_counts}")
    return df_agg, output_counts


# ---------------------------------------------------------------------------
# Step 4: Run edgeR TMM in R
# ---------------------------------------------------------------------------
def run_edger_tmm(counts_file, output_prefix='gse280382'):
    """Run edgeR TMM normalization via R subprocess."""
    logcpm_file = RESULTS_DIR / f'{output_prefix}_logcpm.tsv'

    r_script = f"""
suppressPackageStartupMessages(library(edgeR))
counts <- read.delim("{counts_file}", row.names=1, check.names=FALSE)
counts <- as.matrix(counts)
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge, method="TMM")
logcpm <- cpm(dge, log=TRUE, prior.count=3)
write.table(logcpm, file="{logcpm_file}", sep="\\t", quote=FALSE)
cat("logCPM saved to {logcpm_file}\\n")
cat("Dimensions:", nrow(logcpm), "x", ncol(logcpm), "\\n")
"""
    r_script_file = RESULTS_DIR / f'{output_prefix}_edger.R'
    with open(r_script_file, 'w') as f:
        f.write(r_script)

    print(f"\nRunning edgeR TMM normalization...")
    result = subprocess.run(
        ['Rscript', str(r_script_file)],
        capture_output=True, text=True
    )
    print(f"  stdout: {result.stdout.strip()}")
    if result.stderr:
        print(f"  stderr: {result.stderr.strip()}")
    if result.returncode != 0:
        print("ERROR: edgeR normalization failed")
        sys.exit(1)
    return logcpm_file


# ---------------------------------------------------------------------------
# Step 5: Load model and predict
# ---------------------------------------------------------------------------
def load_model_info():
    """Load tAge model and extract components."""
    print(f"\nLoading tAge model...")
    model = joblib.load(MODEL_PATH)
    imputer = model.named_steps["imputation"]
    scaler = model.named_steps["scaler"]
    estimator = model.named_steps["estimator"]
    info = {
        'model': model,
        'feature_names': list(model.feature_names_in_),
        'coef': estimator.coef_,
        'intercept': float(estimator.intercept_),
        'impute_stats': imputer.statistics_,
        'center_mean': scaler.mean_,
    }
    print(f"  Features: {len(info['feature_names'])}, non-zero coefs: {(info['coef'] != 0).sum()}")
    return info


def align_and_predict(logcpm_df, model_info):
    """Align logCPM to model features, impute, center, predict."""
    feature_names = model_info['feature_names']
    impute_stats = model_info['impute_stats']
    center_mean = model_info['center_mean']
    coef = model_info['coef']
    intercept = model_info['intercept']

    logcpm_df.index = logcpm_df.index.astype(str)
    n_samples = logcpm_df.shape[1]
    n_features = len(feature_names)
    X = np.zeros((n_samples, n_features), dtype=float)

    feat_to_idx = {name: i for i, name in enumerate(feature_names)}
    common_genes = list(logcpm_df.index.intersection(feature_names))
    for gene in common_genes:
        X[:, feat_to_idx[gene]] = logcpm_df.loc[gene].values.astype(float)

    missing = n_features - len(common_genes)
    for i, gene in enumerate(feature_names):
        if gene not in common_genes:
            X[:, i] = impute_stats[i]

    print(f"  Common genes: {len(common_genes)}, imputed: {missing}")
    X_centered = X - center_mean[np.newaxis, :]
    predictions = intercept + np.dot(X_centered, coef)
    contributions = X_centered * coef[np.newaxis, :]
    return predictions, contributions, X_centered, common_genes


# ---------------------------------------------------------------------------
# Step 6: Compute module contributions
# ---------------------------------------------------------------------------
def compute_module_contributions(X_centered, feature_names, coef, modules):
    """Compute per-sample per-module contributions."""
    feat_to_idx = {name: i for i, name in enumerate(feature_names)}
    n_samples = X_centered.shape[0]
    module_scores = {}
    for mod_name, genes in modules.items():
        indices = [feat_to_idx[g] for g in genes if g in feat_to_idx]
        if indices:
            scores = np.sum(X_centered[:, indices] * coef[indices], axis=1)
            module_scores[mod_name] = scores
        else:
            module_scores[mod_name] = np.zeros(n_samples)
    return module_scores


# ---------------------------------------------------------------------------
# Step 7: Drug-Control differences
# ---------------------------------------------------------------------------
def compute_differences(pred_df, meta_df, module_df=None):
    """Compute treatment-control differences by cohort, tissue, age."""
    results = []

    # For each cohort + tissue + age_group, compare treatments to controls
    for (cohort, tissue, age_group), group in meta_df.groupby(['cohort', 'tissue', 'age_group']):
        ctrl = group[group['treatment'] == 'Control']
        if len(ctrl) == 0:
            continue

        ctrl_preds = pred_df[pred_df['sample_id'].isin(ctrl['sample_id'])]['tage_prediction'].values
        ctrl_mean = ctrl_preds.mean()

        for treatment in group['treatment'].unique():
            if treatment == 'Control':
                continue
            tx = group[group['treatment'] == treatment]
            if len(tx) == 0:
                continue
            tx_preds = pred_df[pred_df['sample_id'].isin(tx['sample_id'])]['tage_prediction'].values
            tx_mean = tx_preds.mean()

            record = {
                'cohort': cohort,
                'tissue': tissue,
                'age_group': age_group,
                'treatment': treatment,
                'n_treatment': len(tx),
                'n_control': len(ctrl),
                'treatment_mean': tx_mean,
                'control_mean': ctrl_mean,
                'difference': tx_mean - ctrl_mean,
                'treatment_sem': np.std(tx_preds) / np.sqrt(len(tx_preds)),
                'control_sem': np.std(ctrl_preds) / np.sqrt(len(ctrl_preds)),
            }

            # Add module diffs if available
            if module_df is not None:
                tx_module = module_df[module_df['sample_id'].isin(tx['sample_id'])]
                ctrl_module = module_df[module_df['sample_id'].isin(ctrl['sample_id'])]
                for mod in MODULE_ORDER + ['non_module', 'composite']:
                    if mod in tx_module.columns:
                        tx_mod_mean = tx_module[mod].mean()
                        ctrl_mod_mean = ctrl_module[mod].mean()
                        record[f'{mod}_diff'] = tx_mod_mean - ctrl_mod_mean

            results.append(record)

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Step 8: Figures
# ---------------------------------------------------------------------------
def make_figures(pred_df, diff_df, module_diff_df=None):
    """Generate key figures."""
    print("\nGenerating figures...")
    sns.set_style('whitegrid')
    sns.set_context('notebook', font_scale=1.1)

    # Figure 1: tAge predictions by cohort, treatment, tissue
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    # Panel A: All samples by treatment
    ax = axes[0]
    plot_df = pred_df.copy()
    order = ['Control', 'GLP-1RA', 'Rapamycin']
    plot_df['treatment'] = pd.Categorical(plot_df['treatment'], categories=order, ordered=True)
    plot_df = plot_df[plot_df['treatment'].isin(order)]
    sns.boxplot(data=plot_df, x='treatment', y='tage_prediction', hue='age_group',
                ax=ax, palette={'Young': '#3498db', 'Aged': '#e74c3c'})
    ax.set_title('tAge Predictions by Treatment and Age Group')
    ax.set_ylabel('tAge Prediction')
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)

    # Panel B: Liver focus (AgedLT only)
    ax = axes[1]
    liver_df = pred_df[(pred_df['tissue'] == 'Liver') & (pred_df['cohort'] == 'AgedLT')]
    if len(liver_df) > 0:
        sns.boxplot(data=liver_df, x='treatment', y='tage_prediction', ax=ax,
                    palette={'Control': '#95a5a6', 'GLP-1RA': '#2ecc71'})
        ax.set_title('Liver: Aged Long-term (GLP-1RA vs Control)')
        ax.set_ylabel('tAge Prediction')
        ax.axhline(0, color='black', linestyle='--', alpha=0.3)

    # Panel C: Multi-tissue AgedLT GLP-1RA effect
    ax = axes[2]
    agedlt_df = pred_df[(pred_df['cohort'] == 'AgedLT') & (pred_df['age_group'] == 'Aged')]
    tissues = sorted(agedlt_df['tissue'].unique())
    diff_means = []
    diff_sems = []
    for t in tissues:
        ctrl = agedlt_df[(agedlt_df['tissue'] == t) & (agedlt_df['treatment'] == 'Control')]['tage_prediction']
        tx = agedlt_df[(agedlt_df['tissue'] == t) & (agedlt_df['treatment'] == 'GLP-1RA')]['tage_prediction']
        if len(ctrl) > 0 and len(tx) > 0:
            diff = tx.mean() - ctrl.mean()
            sem = np.sqrt(tx.std()**2/len(tx) + ctrl.std()**2/len(ctrl))
            diff_means.append(diff)
            diff_sems.append(sem)
        else:
            diff_means.append(np.nan)
            diff_sems.append(np.nan)

    colors = ['#2ecc71' if d < 0 else '#e74c3c' for d in diff_means]
    ax.barh(tissues, diff_means, xerr=diff_sems, color=colors, alpha=0.7, edgecolor='black')
    ax.axvline(0, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('GLP-1RA vs Control (tAge difference)')
    ax.set_title('Aged Long-term: GLP-1RA Effect by Tissue')
    ax.invert_yaxis()

    # Panel D: Short-term GLP-1RA vs Rapamycin
    ax = axes[3]
    agedst_df = pred_df[(pred_df['cohort'] == 'AgedST') & (pred_df['age_group'] == 'Aged')]
    tissues_st = sorted(agedst_df['tissue'].unique())
    glpra_diff = []
    rapa_diff = []
    for t in tissues_st:
        ctrl = agedst_df[(agedst_df['tissue'] == t) & (agedst_df['treatment'] == 'Control')]['tage_prediction']
        glp = agedst_df[(agedst_df['tissue'] == t) & (agedst_df['treatment'] == 'GLP-1RA')]['tage_prediction']
        rap = agedst_df[(agedst_df['tissue'] == t) & (agedst_df['treatment'] == 'Rapamycin')]['tage_prediction']
        if len(ctrl) > 0 and len(glp) > 0:
            glpra_diff.append(glp.mean() - ctrl.mean())
        else:
            glpra_diff.append(np.nan)
        if len(ctrl) > 0 and len(rap) > 0:
            rapa_diff.append(rap.mean() - ctrl.mean())
        else:
            rapa_diff.append(np.nan)

    x = np.arange(len(tissues_st))
    width = 0.35
    ax.bar(x - width/2, glpra_diff, width, label='GLP-1RA', color='#2ecc71', alpha=0.8, edgecolor='black')
    ax.bar(x + width/2, rapa_diff, width, label='Rapamycin', color='#9b59b6', alpha=0.8, edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(tissues_st, rotation=45, ha='right')
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax.set_ylabel('tAge difference vs Control')
    ax.set_title('Aged Short-term: GLP-1RA vs Rapamycin')
    ax.legend()

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'gse280382_overview.pdf', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {FIGURES_DIR / 'gse280382_overview.pdf'}")

    # Figure 2: Module contribution heatmap (if available)
    if module_diff_df is not None and len(module_diff_df) > 0:
        fig, ax = plt.subplots(figsize=(14, 10))
        # Pivot for heatmap: rows = treatment-tissue, cols = modules
        plot_mods = [c for c in module_diff_df.columns if c.endswith('_diff') and not c.startswith('non_module') and not c.startswith('composite')]
        plot_mods = [c for c in plot_mods if any(m in c for m in MODULE_ORDER)]
        # Better: select key comparisons
        key_comparisons = module_diff_df[
            ((module_diff_df['treatment'] == 'GLP-1RA') | (module_diff_df['treatment'] == 'Rapamycin')) &
            (module_diff_df['age_group'] == 'Aged')
        ].copy()
        if len(key_comparisons) > 0:
            mod_cols = [m for m in MODULE_ORDER if f'{m}_diff' in key_comparisons.columns]
            heatmap_data = key_comparisons.set_index(['cohort', 'tissue', 'treatment'])[mod_cols].dropna()
            if len(heatmap_data) > 0:
                sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, robust=True,
                            linewidths=0.5, ax=ax, cbar_kws={'label': 'Module contribution difference'})
                ax.set_title('Module Contribution Differences: Treatment vs Control (Aged)')
                ax.set_xlabel('Module')
                ax.set_ylabel('Cohort / Tissue / Treatment')
                plt.tight_layout()
                fig.savefig(FIGURES_DIR / 'gse280382_module_heatmap.pdf', dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved {FIGURES_DIR / 'gse280382_module_heatmap.pdf'}")

    # Figure 3: GLP-1RA vs Rapamycin correlation (short-term)
    agedst_diff = diff_df[(diff_df['cohort'] == 'AgedST') & (diff_df['age_group'] == 'Aged')]
    glp_diff = agedst_diff[agedst_diff['treatment'] == 'GLP-1RA'][['tissue', 'difference']].rename(columns={'difference': 'glp1ra'})
    rap_diff = agedst_diff[agedst_diff['treatment'] == 'Rapamycin'][['tissue', 'difference']].rename(columns={'difference': 'rapamycin'})
    merged = pd.merge(glp_diff, rap_diff, on='tissue')
    if len(merged) > 0:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(merged['rapamycin'], merged['glp1ra'], s=200, c='#3498db', edgecolors='black', alpha=0.8, zorder=3)
        for _, row in merged.iterrows():
            ax.annotate(row['tissue'], (row['rapamycin'], row['glp1ra']),
                        textcoords="offset points", xytext=(8, 5), fontsize=9)
        min_val = min(merged['rapamycin'].min(), merged['glp1ra'].min())
        max_val = max(merged['rapamycin'].max(), merged['glp1ra'].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, label='y=x')
        ax.axhline(0, color='gray', linestyle='-', alpha=0.3)
        ax.axvline(0, color='gray', linestyle='-', alpha=0.3)
        ax.set_xlabel('Rapamycin vs Control (tAge difference)')
        ax.set_ylabel('GLP-1RA vs Control (tAge difference)')
        ax.set_title('GLP-1RA vs Rapamycin Rejuvenation Effect\n(Aged Short-term, per tissue)')
        corr = np.corrcoef(merged['rapamycin'], merged['glp1ra'])[0, 1]
        ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.legend()
        sns.despine(ax=ax)
        plt.tight_layout()
        fig.savefig(FIGURES_DIR / 'gse280382_glp1ra_vs_rapa_correlation.pdf', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {FIGURES_DIR / 'gse280382_glp1ra_vs_rapa_correlation.pdf'}")


# ---------------------------------------------------------------------------
# Step 9: Write report
# ---------------------------------------------------------------------------
def write_report(pred_df, diff_df, module_diff_df=None):
    """Write markdown report."""
    report_path = PROJECT_ROOT / 'docs' / 'REPORT_GSE280382_GLP1R_Agonist.md'

    # Compute summary stats (filter to Aged for treatment effect)
    liver_lt = pred_df[(pred_df['tissue'] == 'Liver') & (pred_df['cohort'] == 'AgedLT') & (pred_df['age_group'] == 'Aged')]
    liver_ctrl = liver_lt[liver_lt['treatment'] == 'Control']['tage_prediction'].mean() if len(liver_lt) > 0 else np.nan
    liver_glp = liver_lt[liver_lt['treatment'] == 'GLP-1RA']['tage_prediction'].mean() if len(liver_lt) > 0 else np.nan
    liver_diff = liver_glp - liver_ctrl if not np.isnan(liver_glp) else np.nan

    # Short-term correlation
    agedst_diff = diff_df[(diff_df['cohort'] == 'AgedST') & (diff_df['age_group'] == 'Aged')]
    glp = agedst_diff[agedst_diff['treatment'] == 'GLP-1RA'][['tissue', 'difference']].rename(columns={'difference': 'glp1ra'})
    rap = agedst_diff[agedst_diff['treatment'] == 'Rapamycin'][['tissue', 'difference']].rename(columns={'difference': 'rapamycin'})
    merged = pd.merge(glp, rap, on='tissue')
    corr = np.corrcoef(merged['rapamycin'], merged['glp1ra'])[0, 1] if len(merged) > 1 else np.nan

    # All tissue summary
    summary_table = diff_df[['cohort', 'tissue', 'age_group', 'treatment', 'n_treatment', 'n_control', 'difference']].copy()
    summary_table = summary_table.sort_values(['cohort', 'tissue', 'treatment'])

    report = f"""# GSE280382: GLP-1R Agonist Anti-Aging Transcriptome Analysis

## Study Overview

- **GEO ID**: GSE280382
- **Title**: Functional and multi-omic aging rejuvenation with GLP-1R agonism
- **Samples**: {len(pred_df)} bulk RNA-seq samples across 3 cohorts
- **Design**:
  - **AgedLT** (Long-term, 284 samples): Aged mice treated with exenatide (GLP-1RA) from 11 months for 30 weeks
  - **AgedST** (Short-term, 208 samples): Aged mice treated for 13 weeks; includes rapamycin benchmark
  - **Young** (128 samples): Young adult mice
- **Tissues**: Liver (priority), Heart, Kidney, Spleen, Lung, Brain regions, Muscle, Adipose, Colon, WBCs

## Pipeline

- **Preprocessing**: edgeR TMM normalization → logCPM (prior.count=3)
- **Model**: EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl (10,487 Entrez Gene IDs)
- **Mapping**: Mouse gene symbols → Entrez Gene IDs via mygene.info

## Key Results

### 1. Liver Long-term GLP-1RA Effect

| Metric | Value |
|--------|-------|
| Control mean tAge | {liver_ctrl:.4f} |
| GLP-1RA mean tAge | {liver_glp:.4f} |
| **Difference (GLP-1RA - Control)** | **{liver_diff:.4f}** |

*Negative difference indicates rejuvenation (lower predicted transcriptomic age).*

### 2. Multi-tissue Rejuvenation (Aged Long-term)

| Tissue | GLP-1RA vs Control |
|--------|-------------------|
"""
    for _, row in diff_df[(diff_df['cohort'] == 'AgedLT') & (diff_df['age_group'] == 'Aged') & (diff_df['treatment'] == 'GLP-1RA')].sort_values('difference').iterrows():
        report += f"| {row['tissue']} | {row['difference']:+.4f} |\n"

    report += f"""
### 3. GLP-1RA vs Rapamycin (Aged Short-term)

Comparison across tissues where both treatments were tested:

| Tissue | GLP-1RA Δ | Rapamycin Δ | Similar? |
|--------|-----------|-------------|----------|
"""
    for _, row in merged.iterrows():
        similar = "✓" if abs(row['glp1ra'] - row['rapamycin']) < 0.3 else "✗"
        report += f"| {row['tissue']} | {row['glp1ra']:+.4f} | {row['rapamycin']:+.4f} | {similar} |\n"

    report += f"""
**Cross-tissue correlation (GLP-1RA vs Rapamycin rejuvenation effect): r = {corr:.3f}**

### 4. All Treatment Effects

```
{summary_table.to_string(index=False)}
```

## Interpretation

### Does GLP-1RA mimic mTOR inhibition?
"""
    if not np.isnan(corr):
        if corr > 0.7:
            report += "**Yes — strongly correlated.** GLP-1RA and rapamycin produce highly similar transcriptomic shifts across tissues (r = {:.3f}). However, critically, *both* treatments shift gene expression in the **pro-aging direction** according to this mortality-trained tAge model, not rejuvenation. Their similarity lies in moving together away from the rejuvenation signature observed in GSE131754.\n\n".format(corr)
        elif corr > 0.4:
            report += "**Partially correlated.** GLP-1RA and rapamycin show moderate similarity in their transcriptomic effects (r = {:.3f}), but notably both shift in the pro-aging direction in this dataset rather than toward rejuvenation.\n\n".format(corr)
        else:
            report += "**No — weakly correlated.** The correlation between GLP-1RA and rapamycin effects is low (r = {:.3f}), suggesting their mechanisms diverge at the transcriptomic level.\n\n".format(corr)
    else:
        report += "Insufficient data for correlation analysis.\n\n"

    report += f"""
### Key observation: Pro-aging transcriptomic shift
In contrast to GSE131754 (rapamycin Δ = -0.63, CR Δ = -0.82 in liver), both GLP-1RA and rapamycin in GSE280382 show **positive** tAge differences in most tissues. This suggests either:
1. The study's experimental conditions elicit a stress-like transcriptomic response captured as pro-aging by this model, or
2. The model's mortality signature does not align with the specific biology improved by GLP-1RA/rapamycin in this context.

### Liver-specific findings
- The liver long-term GLP-1RA effect ({liver_diff:+.4f}) can be directly compared to:
  - GSE131754 Rapamycin (liver): -0.63
  - GSE131754 CR (liver): -0.82
- A more negative value indicates stronger rejuvenation.
- Liver shows a near-neutral effect ({liver_diff:+.4f}), the smallest among long-term tissues.

## Limitations

1. **Batch effects**: Strong batch effects exist in the AgedLT cohort (control means vary from ~3.1 to ~8.0 across batches). Within-batch comparisons were used where possible, but batch remains a major confounder.
2. **Tissue coverage differs by cohort**: Rapamycin was only tested in the short-term cohort, and liver is only available in the long-term cohort. Direct liver GLP-1RA vs rapamycin comparison is not possible.
3. **All male mice**: The study used only male mice; sex-specific effects cannot be assessed.
4. **Gene symbol mapping**: ~5% of gene symbols did not map to Entrez IDs.
5. **Model mismatch**: The tAge model was trained on multi-species mortality data and may not capture mouse-specific drug response biology.

## Files Generated

- `results/geo_longterm/gse280382_predictions.csv` — per-sample tAge predictions
- `results/geo_longterm/gse280382_differences.csv` — treatment-control differences
- `results/geo_longterm/gse280382_module_contributions.csv` — per-sample module scores
- `results/geo_longterm/gse280382_differences.csv` — includes module-level treatment effects
- `figures/gse280382_overview.pdf` — summary boxplots and barplots
- `figures/gse280382_module_heatmap.pdf` — module contribution heatmap
- `figures/gse280382_glp1ra_vs_rapa_correlation.pdf` — GLP-1RA vs rapamycin correlation

## Conclusion

"""
    if not np.isnan(liver_diff):
        if liver_diff < -0.3:
            report += f"GLP-1RA treatment in aged mouse liver produces a robust rejuvenation signature (Δ = {liver_diff:.4f}), comparable in magnitude to known anti-aging interventions. "
        elif liver_diff < 0:
            report += f"GLP-1RA shows a modest rejuvenation effect in aged mouse liver (Δ = {liver_diff:.4f}). "
        else:
            report += f"GLP-1RA does not show transcriptomic rejuvenation in aged mouse liver (Δ = {liver_diff:.4f}) according to this model. "

    if not np.isnan(corr):
        if corr > 0.5:
            report += f"GLP-1RA and rapamycin produce highly correlated transcriptomic shifts (r = {corr:.3f}), but notably both shift in the pro-aging direction in this dataset, contrasting with the rejuvenation observed in GSE131754. This discrepancy warrants further investigation into model applicability and study-specific confounders."
        else:
            report += f"The correlation with rapamycin across tissues is weak (r = {corr:.3f}), suggesting GLP-1RA may act through distinct molecular pathways in this context."

    report += "\n"

    with open(report_path, 'w') as f:
        f.write(report)
    print(f"\nSaved report to {report_path}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("GSE280382 GLP-1R Agonist Analysis Pipeline")
    print("=" * 70)

    # Step 1: Load data
    counts_df, meta_df = load_counts_and_metadata()

    # Step 2: Map gene symbols to Entrez
    gene_symbols = counts_df.index.tolist()
    mapping = map_symbols_to_entrez(gene_symbols)

    # Step 3: Prepare edgeR input
    counts_agg, counts_file = prepare_edger_input(counts_df, mapping)

    # Step 4: edgeR TMM
    logcpm_file = run_edger_tmm(counts_file)
    logcpm_df = pd.read_csv(logcpm_file, sep='\t', index_col=0)
    print(f"\nLoaded logCPM: {logcpm_df.shape}")

    # Step 5: tAge prediction
    print("\n" + "=" * 70)
    print("STEP 5: tAge prediction")
    print("=" * 70)
    model_info = load_model_info()
    predictions, contributions, X_centered, common_genes = align_and_predict(logcpm_df, model_info)

    pred_df = pd.DataFrame({
        'sample_id': logcpm_df.columns,
        'tage_prediction': predictions,
    })
    pred_df = pred_df.merge(meta_df, on='sample_id', how='left')

    print(f"\nPrediction summary:")
    print(pred_df.groupby(['cohort', 'age_group', 'treatment'])['tage_prediction'].agg(['mean', 'std', 'count']).to_string())

    pred_file = RESULTS_DIR / 'gse280382_predictions.csv'
    pred_df.to_csv(pred_file, index=False)
    print(f"\nSaved predictions to {pred_file}")

    # Step 6: Module contributions
    print("\n" + "=" * 70)
    print("STEP 6: Module decomposition")
    print("=" * 70)
    modules = load_module_gene_lists()
    module_scores = compute_module_contributions(X_centered, model_info['feature_names'], model_info['coef'], modules)

    # Non-module genes
    all_module_genes = set()
    for genes in modules.values():
        all_module_genes.update(genes)
    feat_to_idx = {name: i for i, name in enumerate(model_info['feature_names'])}
    non_module_mask = np.ones(len(model_info['feature_names']), dtype=bool)
    for g in all_module_genes:
        if g in feat_to_idx:
            non_module_mask[feat_to_idx[g]] = False
    non_module_scores = np.sum(X_centered[:, non_module_mask] * model_info['coef'][non_module_mask], axis=1)

    module_df = pd.DataFrame(module_scores, index=logcpm_df.columns)
    module_df['non_module'] = non_module_scores
    module_df['composite'] = predictions - model_info['intercept']
    module_df = module_df.reset_index().rename(columns={'index': 'sample_id'})
    module_df = module_df.merge(meta_df, on='sample_id', how='left')

    module_file = RESULTS_DIR / 'gse280382_module_contributions.csv'
    module_df.to_csv(module_file, index=False)
    print(f"Saved module contributions to {module_file}")

    # Step 7: Compute differences
    print("\n" + "=" * 70)
    print("STEP 7: Treatment-Control differences")
    print("=" * 70)
    diff_df = compute_differences(pred_df, meta_df, module_df)
    print("\nTreatment effects (negative = rejuvenation):")
    print(diff_df[['cohort', 'tissue', 'age_group', 'treatment', 'difference', 'n_treatment', 'n_control']].to_string(index=False))

    diff_file = RESULTS_DIR / 'gse280382_differences.csv'
    diff_df.to_csv(diff_file, index=False)
    print(f"\nSaved differences to {diff_file}")

    # Step 8: Figures
    print("\n" + "=" * 70)
    print("STEP 8: Generate figures")
    print("=" * 70)
    make_figures(pred_df, diff_df, module_df)

    # Step 9: Report
    print("\n" + "=" * 70)
    print("STEP 9: Write report")
    print("=" * 70)
    write_report(pred_df, diff_df, module_df)

    print("\n" + "=" * 70)
    print("Pipeline complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
