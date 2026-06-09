#!/usr/bin/env python3
"""
GSE288795 Analysis: Trametinib + Rapamycin Combination
========================================================
Process GSE288795 RNA-seq count data, run tAge prediction,
and test additivity of trametinib + rapamycin combination.

Study: "The geroprotectors trametinib and rapamycin combine additively
to extend mouse healthspan and lifespan" (Gkioni et al., 2025)

Tissues: Muscle, Kidney, Spleen (NOT liver — collected at 24 months)
Samples: 111 (bulk RNA-seq)
Design: Control, Rapamycin, Trametinib, Combo; Male & Female

Author: Automated pipeline
Date: 2026-06-09
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
SUPP_TABLE = PROJECT_ROOT / 'data' / 'raw' / 'supp_table_7.xlsx'

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


# ============================================================================
# STEP 1: Load GSE288795 counts
# ============================================================================

def load_gse288795_counts():
    """Load GSE288795 count matrix from Excel."""
    filepath = DATA_DIR / 'GSE288795_count_data.xlsx'
    print(f"Loading GSE288795 from {filepath}")
    
    # Row 0 = sample titles, Row 1 = library names, Row 2+ = data
    df = pd.read_excel(filepath, sheet_name='RNA Seq count data',
                       header=0, index_col=0, skiprows=[1])
    
    # Remove the first row if it's still there (shouldn't be with skiprows)
    # Ensure index is gene IDs
    df.index.name = 'gene_id'
    
    print(f"  Shape: {df.shape}")
    print(f"  Genes: {df.shape[0]:,}")
    print(f"  Samples: {df.shape[1]}")
    print(f"  Sample groups: {df.columns.tolist()[:5]}...")
    
    return df


def parse_sample_metadata_gse288795(columns):
    """Parse GSE288795 column names into metadata."""
    meta = []
    for col in columns:
        parts = col.split()
        # Format: "Sex Tissue Treatment Biol Rep N"
        # e.g., "Male Muscle Rapamycin Biol Rep 1"
        # e.g., "Female Muscle Rapamycin/Trametinib Biol Rep 1"
        sex = parts[0]
        tissue = parts[1]
        
        # Treatment may be "Rapamycin", "Trametinib", "Rapamycin/Trametinib", or "Control"
        if parts[2] == 'Rapamycin/Trametinib':
            treatment = 'Combo'
        elif parts[2] == 'Rapamycin':
            treatment = 'Rapamycin'
        elif parts[2] == 'Trametinib':
            treatment = 'Trametinib'
        elif parts[2] == 'Control':
            treatment = 'Control'
        else:
            treatment = parts[2]
        
        # Replicate number is last element
        rep = parts[-1]
        
        meta.append({
            'sample_id': col,
            'sex': sex,
            'tissue': tissue,
            'treatment': treatment,
            'replicate': rep,
            'age': '24m',
        })
    
    return pd.DataFrame(meta)


# ============================================================================
# STEP 2: Ensembl ID → Entrez Gene ID mapping (reuse cached)
# ============================================================================

def load_cached_mapping():
    """Load cached Ensembl→Entrez mapping from GSE131754."""
    cache_file = DATA_DIR / 'ensembl_to_entrez_mapping_gse131754.pkl'
    if os.path.exists(cache_file):
        print(f"Loading cached mapping from {cache_file}")
        with open(cache_file, 'rb') as f:
            mapping = pickle.load(f)
        print(f"  Cached mapping: {len(mapping)} genes")
        return mapping
    else:
        print("ERROR: No cached mapping found.")
        sys.exit(1)


# ============================================================================
# STEP 3: Prepare data for edgeR TMM normalization
# ============================================================================

def prepare_edger_input(counts_df, mapping, output_prefix='gse288795'):
    """
    Prepare counts matrix for edgeR TMM normalization in R.
    """
    print(f"\nPreparing edgeR input...")
    
    df = counts_df.copy()
    df['entrez'] = df.index.map(lambda x: mapping.get(x.split('.')[0], None))
    
    df_mapped = df.dropna(subset=['entrez']).copy()
    df_mapped['entrez'] = df_mapped['entrez'].astype(int)
    
    print(f"  Mapped genes: {len(df_mapped)}")
    
    df_agg = df_mapped.groupby('entrez').sum()
    if 'entrez' in df_agg.columns:
        df_agg = df_agg.drop(columns=['entrez'])
    
    print(f"  After aggregation: {len(df_agg)} unique Entrez IDs")
    
    output_counts = RESULTS_DIR / f'{output_prefix}_counts_for_edger.tsv'
    df_agg.to_csv(output_counts, sep='\t')
    print(f"  Saved to {output_counts}")
    
    return df_agg, output_counts


def run_edger_tmm(counts_file, output_prefix='gse288795'):
    """Run edgeR TMM normalization in R via subprocess."""
    
    logcpm_file = RESULTS_DIR / f'{output_prefix}_logcpm.tsv'
    
    r_script = f"""
suppressPackageStartupMessages(library(edgeR))

# Read counts
counts <- read.delim("{counts_file}", row.names=1, check.names=FALSE)
counts <- as.matrix(counts)

# Create DGEList
dge <- DGEList(counts=counts)

# TMM normalization
dge <- calcNormFactors(dge, method="TMM")

# logCPM (prior.count=3 as in tAge pipeline)
logcpm <- cpm(dge, log=TRUE, prior.count=3)

# Save
write.table(logcpm, file="{logcpm_file}", sep="\\t", quote=FALSE)
cat("logCPM saved to {logcpm_file}\\n")
cat("Dimensions:", nrow(logcpm), "x", ncol(logcpm), "\\n")
"""
    
    r_script_file = RESULTS_DIR / f'{output_prefix}_edger.R'
    with open(r_script_file, 'w') as f:
        f.write(r_script)
    
    print(f"\nRunning edgeR TMM normalization...")
    print(f"  R script: {r_script_file}")
    
    result = subprocess.run(
        ['Rscript', str(r_script_file)],
        capture_output=True,
        text=True
    )
    
    print(f"  stdout: {result.stdout}")
    if result.stderr:
        print(f"  stderr: {result.stderr}")
    
    if result.returncode != 0:
        print("ERROR: edgeR normalization failed")
        sys.exit(1)
    
    return logcpm_file


# ============================================================================
# STEP 4: tAge prediction
# ============================================================================

def load_tage_model(model_path):
    """Load the pre-trained tAge model."""
    print(f"\nLoading tAge model from {model_path}")
    
    model = joblib.load(model_path)
    
    imputer = model.named_steps["imputation"]
    scaler = model.named_steps["scaler"]
    estimator = model.named_steps["estimator"]
    
    feature_names = list(model.feature_names_in_)
    coef = estimator.coef_
    intercept = float(estimator.intercept_)
    impute_stats = imputer.statistics_
    center_mean = scaler.mean_
    
    print(f"  Features: {len(feature_names)}")
    print(f"  Non-zero coefs: {(coef != 0).sum()}")
    print(f"  Intercept: {intercept:.4f}")
    
    return {
        'model': model,
        'coef': coef,
        'intercept': intercept,
        'feature_names': feature_names,
        'impute_stats': impute_stats,
        'center_mean': center_mean,
    }


def align_and_impute(logcpm_df, model_info):
    """Align logCPM matrix with model features."""
    feature_names = model_info['feature_names']
    impute_stats = model_info['impute_stats']
    center_mean = model_info['center_mean']
    
    logcpm_df.index = logcpm_df.index.astype(str)
    
    common_genes = list(logcpm_df.index.intersection(feature_names))
    print(f"  Common genes: {len(common_genes)}")
    
    n_samples = logcpm_df.shape[1]
    n_features = len(feature_names)
    X = np.zeros((n_samples, n_features), dtype=float)
    
    feat_to_idx = {name: i for i, name in enumerate(feature_names)}
    
    for gene in common_genes:
        col_idx = feat_to_idx[gene]
        X[:, col_idx] = logcpm_df.loc[gene].values.astype(float)
    
    missing_count = 0
    for i, gene in enumerate(feature_names):
        if gene not in common_genes:
            X[:, i] = impute_stats[i]
            missing_count += 1
    
    print(f"  Missing genes imputed: {missing_count}")
    
    X_centered = X - center_mean[np.newaxis, :]
    
    return X_centered, feature_names


def tage_predict(logcpm_df, model_info):
    """Run tAge prediction on logCPM matrix."""
    coef = model_info['coef']
    intercept = model_info['intercept']
    
    X_centered, feature_names = align_and_impute(logcpm_df, model_info)
    
    predictions = intercept + np.dot(X_centered, coef)
    contributions = X_centered * coef[np.newaxis, :]
    
    return predictions, contributions, feature_names, X_centered


# ============================================================================
# STEP 5: Module decomposition
# ============================================================================

def load_module_gene_lists():
    """Load 14 module gene lists from supp_table_7.xlsx."""
    print(f"\nLoading module gene lists from {SUPP_TABLE}")
    df = pd.read_excel(
        SUPP_TABLE,
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


# ============================================================================
# STEP 6: Compute Drug - Control differences
# ============================================================================

def compute_drug_control_differences(predictions, metadata_df, module_df=None):
    """Compute Drug - Control differences per treatment/tissue/sex."""
    results = []
    module_diff_records = []
    
    control_samples = metadata_df[
        metadata_df['treatment'] == 'Control'
    ]['sample_id'].tolist()
    
    control_preds = {s: predictions[s] for s in control_samples if s in predictions}
    
    for tissue in metadata_df['tissue'].unique():
        for sex in metadata_df['sex'].unique():
            for treatment in metadata_df['treatment'].unique():
                if treatment == 'Control':
                    continue
                
                drug_samples = metadata_df[
                    (metadata_df['treatment'] == treatment) &
                    (metadata_df['tissue'] == tissue) &
                    (metadata_df['sex'] == sex)
                ]['sample_id'].tolist()
                
                control_samples_match = metadata_df[
                    (metadata_df['treatment'] == 'Control') &
                    (metadata_df['tissue'] == tissue) &
                    (metadata_df['sex'] == sex)
                ]['sample_id'].tolist()
                
                drug_preds = [predictions[s] for s in drug_samples if s in predictions]
                control_preds_match = [predictions[s] for s in control_samples_match if s in predictions]
                
                if len(drug_preds) == 0 or len(control_preds_match) == 0:
                    continue
                
                drug_mean = np.mean(drug_preds)
                control_mean = np.mean(control_preds_match)
                diff = drug_mean - control_mean
                
                results.append({
                    'treatment': treatment,
                    'tissue': tissue,
                    'sex': sex,
                    'age': '24m',
                    'n_drug': len(drug_preds),
                    'n_control': len(control_preds_match),
                    'drug_mean': drug_mean,
                    'control_mean': control_mean,
                    'difference': diff,
                    'drug_sem': np.std(drug_preds) / np.sqrt(len(drug_preds)),
                    'control_sem': np.std(control_preds_match) / np.sqrt(len(control_preds_match)),
                })
                
                # Module-level differences
                if module_df is not None:
                    drug_module = module_df[module_df['sample_id'].isin(drug_samples)][MODULE_ORDER + ['non_module', 'composite']].mean()
                    ctrl_module = module_df[module_df['sample_id'].isin(control_samples_match)][MODULE_ORDER + ['non_module', 'composite']].mean()
                    mod_diff = drug_module - ctrl_module
                    
                    rec = {
                        'treatment': treatment,
                        'tissue': tissue,
                        'sex': sex,
                        'age': '24m',
                        'n_drug': len(drug_samples),
                        'n_control': len(control_preds_match),
                    }
                    for col in MODULE_ORDER + ['non_module', 'composite']:
                        rec[f'{col}_drug'] = drug_module[col]
                        rec[f'{col}_control'] = ctrl_module[col]
                        rec[f'{col}_diff'] = mod_diff[col]
                    module_diff_records.append(rec)
    
    return pd.DataFrame(results), pd.DataFrame(module_diff_records)


# ============================================================================
# STEP 7: Visualizations
# ============================================================================

def plot_predictions(pred_df, output_file='gse288795_predictions.pdf'):
    """Plot tAge predictions by treatment, tissue, and sex."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    
    tissues = ['Muscle', 'Kidney', 'Spleen']
    colors = {'Control': '#7f7f7f', 'Rapamycin': '#2ca02c', 
              'Trametinib': '#1f77b4', 'Combo': '#d62728'}
    
    for ax, tissue in zip(axes, tissues):
        tissue_df = pred_df[pred_df['tissue'] == tissue]
        
        positions = []
        labels = []
        pos = 0
        
        for sex in ['Male', 'Female']:
            sex_df = tissue_df[tissue_df['sex'] == sex]
            for treatment in ['Control', 'Rapamycin', 'Trametinib', 'Combo']:
                treat_df = sex_df[sex_df['treatment'] == treatment]
                if len(treat_df) > 0:
                    ax.scatter([pos] * len(treat_df), treat_df['tage_prediction'],
                              color=colors[treatment], alpha=0.6, s=50)
                    ax.boxplot(treat_df['tage_prediction'].values, positions=[pos],
                              widths=0.5, patch_artist=True,
                              boxprops=dict(facecolor=colors[treatment], alpha=0.3),
                              medianprops=dict(color='black'))
                positions.append(pos)
                labels.append(f"{sex}\n{treatment}")
                pos += 1
            pos += 0.5  # gap between sexes
        
        ax.set_title(tissue, fontsize=12, fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.set_ylabel('tAge prediction')
    
    fig.suptitle('GSE288795: tAge Predictions by Treatment, Tissue, and Sex', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved figure: {FIGURES_DIR / output_file}")
    plt.close()


def plot_drug_differences(diff_df, output_file='gse288795_drug_differences.pdf'):
    """Plot Drug - Control differences."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    tissues = ['Muscle', 'Kidney', 'Spleen']
    treatments = ['Rapamycin', 'Trametinib', 'Combo']
    colors = {'Rapamycin': '#2ca02c', 'Trametinib': '#1f77b4', 'Combo': '#d62728'}
    hatches = {'Male': '', 'Female': '///'}
    
    x = np.arange(len(tissues))
    width = 0.25
    
    for i, treatment in enumerate(treatments):
        for j, sex in enumerate(['Male', 'Female']):
            vals = []
            errs = []
            for tissue in tissues:
                row = diff_df[
                    (diff_df['treatment'] == treatment) &
                    (diff_df['tissue'] == tissue) &
                    (diff_df['sex'] == sex)
                ]
                if len(row) > 0:
                    vals.append(row['difference'].values[0])
                    errs.append(row['drug_sem'].values[0] + row['control_sem'].values[0])
                else:
                    vals.append(np.nan)
                    errs.append(0)
            
            offset = (i - 1) * width + (j - 0.5) * width * 0.4
            ax.bar(x + offset, vals, width * 0.4, yerr=errs, 
                   label=f'{treatment} {sex}',
                   color=colors[treatment], hatch=hatches[sex],
                   edgecolor='black', alpha=0.8, capsize=3)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.set_ylabel('Drug - Control tAge Difference\n(negative = rejuvenation)', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(tissues, fontsize=11)
    ax.set_title('GSE288795: Rejuvenation Effects by Treatment, Tissue, and Sex', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=8)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved figure: {FIGURES_DIR / output_file}")
    plt.close()


def plot_module_heatmap(module_diff_df, output_file='gse288795_module_heatmap.pdf'):
    """Plot module-level differences as heatmap."""
    if len(module_diff_df) == 0:
        print("  No module diff data for heatmap")
        return
    
    # Create a pivot for each tissue
    tissues = module_diff_df['tissue'].unique()
    fig, axes = plt.subplots(1, len(tissues), figsize=(18, 6), sharey=True)
    if len(tissues) == 1:
        axes = [axes]
    
    for ax, tissue in zip(axes, tissues):
        tissue_df = module_diff_df[module_diff_df['tissue'] == tissue].copy()
        
        # Build row labels
        tissue_df['label'] = tissue_df['sex'] + ' ' + tissue_df['treatment']
        
        diff_cols = [c for c in tissue_df.columns if c.endswith('_diff') and not c.startswith('non_module') and not c.startswith('composite')]
        diff_cols = [c for c in diff_cols if c.replace('_diff', '') in MODULE_ORDER]
        
        # Reorder columns by MODULE_ORDER
        ordered_cols = [f"{m}_diff" for m in MODULE_ORDER if f"{m}_diff" in diff_cols]
        
        pivot = tissue_df.set_index('label')[ordered_cols]
        
        # Shorten column names
        pivot.columns = [c.replace('_diff', '') for c in pivot.columns]
        
        sns.heatmap(pivot, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                   linewidths=0.5, ax=ax, cbar_kws={'label': 'Drug - Control'})
        ax.set_title(f'{tissue}', fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('')
    
    fig.suptitle('GSE288795: Module-Level Rejuvenation Effects', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved figure: {FIGURES_DIR / output_file}")
    plt.close()


def plot_additivity_test(diff_df, output_file='gse288795_additivity_test.pdf'):
    """Test additivity: is Combo ≈ Rapa + Tram?"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    tissues = ['Muscle', 'Kidney', 'Spleen']
    
    for ax, tissue in zip(axes, tissues):
        tissue_diff = diff_df[diff_df['tissue'] == tissue]
        
        for sex in ['Male', 'Female']:
            sex_diff = tissue_diff[tissue_diff['sex'] == sex]
            
            rapa = sex_diff[sex_diff['treatment'] == 'Rapamycin']['difference'].values
            tram = sex_diff[sex_diff['treatment'] == 'Trametinib']['difference'].values
            combo = sex_diff[sex_diff['treatment'] == 'Combo']['difference'].values
            
            if len(rapa) == 0 or len(tram) == 0 or len(combo) == 0:
                continue
            
            expected = rapa[0] + tram[0]
            observed = combo[0]
            
            color = '#1f77b4' if sex == 'Male' else '#ff7f0e'
            ax.scatter(expected, observed, color=color, s=150, label=sex, edgecolors='black', zorder=5)
            
            # Error bars
            rapa_err = sex_diff[sex_diff['treatment'] == 'Rapamycin']['drug_sem'].values[0] + sex_diff[sex_diff['treatment'] == 'Rapamycin']['control_sem'].values[0]
            tram_err = sex_diff[sex_diff['treatment'] == 'Trametinib']['drug_sem'].values[0] + sex_diff[sex_diff['treatment'] == 'Trametinib']['control_sem'].values[0]
            combo_err = sex_diff[sex_diff['treatment'] == 'Combo']['drug_sem'].values[0] + sex_diff[sex_diff['treatment'] == 'Combo']['control_sem'].values[0]
            
            expected_err = np.sqrt(rapa_err**2 + tram_err**2)
            ax.errorbar(expected, observed, xerr=expected_err, yerr=combo_err,
                       fmt='none', color=color, alpha=0.5, capsize=4)
        
        # Diagonal (perfect additivity)
        lims = ax.get_xlim() + ax.get_ylim()
        min_val = min(lims)
        max_val = max(lims)
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=1, label='Perfect additivity')
        
        ax.set_xlabel('Expected (Rapa + Tram)', fontsize=10)
        ax.set_ylabel('Observed (Combo)', fontsize=10)
        ax.set_title(tissue, fontsize=11, fontweight='bold')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    
    fig.suptitle('GSE288795: Additivity Test — Is Combo ≈ Rapa + Trametinib?', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved figure: {FIGURES_DIR / output_file}")
    plt.close()


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("GSE288795: Trametinib + Rapamycin Combination Analysis")
    print("=" * 70)
    
    # --- Step 1: Load data ---
    print("\n" + "=" * 70)
    print("STEP 1: Load GSE288795 counts")
    print("=" * 70)
    counts_df = load_gse288795_counts()
    
    metadata = parse_sample_metadata_gse288795(counts_df.columns)
    print(f"\nSample metadata summary:")
    print(metadata.groupby(['tissue', 'sex', 'treatment']).size().to_string())
    
    # --- Step 2: ID mapping ---
    print("\n" + "=" * 70)
    print("STEP 2: Load Ensembl → Entrez mapping")
    print("=" * 70)
    mapping = load_cached_mapping()
    
    # --- Step 3: Prepare for edgeR ---
    print("\n" + "=" * 70)
    print("STEP 3: Prepare edgeR input")
    print("=" * 70)
    counts_agg, counts_file = prepare_edger_input(counts_df, mapping)
    
    # --- Step 4: TMM normalization ---
    print("\n" + "=" * 70)
    print("STEP 4: edgeR TMM normalization")
    print("=" * 70)
    logcpm_file = run_edger_tmm(counts_file)
    
    print(f"\nLoading logCPM from {logcpm_file}")
    logcpm_df = pd.read_csv(logcpm_file, sep='\t', index_col=0)
    print(f"  Shape: {logcpm_df.shape}")
    
    # --- Step 5: tAge prediction ---
    print("\n" + "=" * 70)
    print("STEP 5: tAge prediction")
    print("=" * 70)
    model_info = load_tage_model(MODEL_PATH)
    
    predictions, gene_contributions, feature_names, X_centered = tage_predict(logcpm_df, model_info)
    
    pred_df = pd.DataFrame({
        'sample_id': logcpm_df.columns,
        'tage_prediction': predictions,
    })
    pred_df = pred_df.merge(metadata, on='sample_id', how='left')
    
    print(f"\nPredictions summary:")
    print(pred_df.groupby(['tissue', 'sex', 'treatment'])['tage_prediction'].agg(['mean', 'std', 'count']).to_string())
    
    pred_file = RESULTS_DIR / 'gse288795_tage_predictions.csv'
    pred_df.to_csv(pred_file, index=False)
    print(f"\nSaved predictions to {pred_file}")
    
    # --- Step 6: Module decomposition ---
    print("\n" + "=" * 70)
    print("STEP 6: Module decomposition")
    print("=" * 70)
    modules = load_module_gene_lists()
    
    module_scores = compute_module_contributions(X_centered, feature_names, model_info['coef'], modules)
    module_df = pd.DataFrame(module_scores, index=logcpm_df.columns)
    module_df.index.name = 'sample_id'
    module_df = module_df.reset_index()
    module_df = module_df.merge(metadata, on='sample_id', how='left')
    
    # Add non-module and composite
    non_module_genes = set()
    for genes in modules.values():
        non_module_genes.update(genes)
    feat_to_idx = {name: i for i, name in enumerate(feature_names)}
    non_module_mask = np.ones(len(feature_names), dtype=bool)
    for g in non_module_genes:
        if g in feat_to_idx:
            non_module_mask[feat_to_idx[g]] = False
    non_module_scores = np.sum(X_centered[:, non_module_mask] * model_info['coef'][non_module_mask], axis=1)
    module_df['non_module'] = non_module_scores
    module_df['composite'] = predictions - model_info['intercept']
    
    module_file = RESULTS_DIR / 'gse288795_per_sample_module_contributions.csv'
    module_df.to_csv(module_file, index=False)
    print(f"  Saved per-sample module contributions to {module_file}")
    
    # --- Step 7: Drug - Control differences ---
    print("\n" + "=" * 70)
    print("STEP 7: Drug - Control differences")
    print("=" * 70)
    pred_dict = dict(zip(logcpm_df.columns, predictions))
    diff_df, module_diff_df = compute_drug_control_differences(pred_dict, metadata, module_df)
    
    print(f"\nDrug - Control differences (negative = rejuvenation):")
    print(diff_df[['treatment', 'tissue', 'sex', 'difference', 'n_drug', 'n_control']].to_string(index=False))
    
    diff_file = RESULTS_DIR / 'gse288795_drug_control_differences.csv'
    diff_df.to_csv(diff_file, index=False)
    print(f"\nSaved differences to {diff_file}")
    
    if len(module_diff_df) > 0:
        module_diff_file = RESULTS_DIR / 'gse288795_module_drug_control_differences.csv'
        module_diff_df.to_csv(module_diff_file, index=False)
        print(f"Saved module differences to {module_diff_file}")
    
    # --- Step 8: Additivity analysis ---
    print("\n" + "=" * 70)
    print("STEP 8: Additivity analysis")
    print("=" * 70)
    
    additivity_results = []
    for tissue in diff_df['tissue'].unique():
        for sex in diff_df['sex'].unique():
            tissue_sex = diff_df[(diff_df['tissue'] == tissue) & (diff_df['sex'] == sex)]
            
            rapa = tissue_sex[tissue_sex['treatment'] == 'Rapamycin']
            tram = tissue_sex[tissue_sex['treatment'] == 'Trametinib']
            combo = tissue_sex[tissue_sex['treatment'] == 'Combo']
            
            if len(rapa) == 0 or len(tram) == 0 or len(combo) == 0:
                continue
            
            expected = rapa['difference'].values[0] + tram['difference'].values[0]
            observed = combo['difference'].values[0]
            deviation = observed - expected
            
            additivity_results.append({
                'tissue': tissue,
                'sex': sex,
                'rapa_diff': rapa['difference'].values[0],
                'tram_diff': tram['difference'].values[0],
                'combo_diff': observed,
                'expected_additive': expected,
                'deviation': deviation,
                'is_synergistic': deviation < -0.1,  # stronger than additive
                'is_antagonistic': deviation > 0.1,   # weaker than additive
            })
    
    add_df = pd.DataFrame(additivity_results)
    print(f"\nAdditivity test:")
    print(add_df[['tissue', 'sex', 'rapa_diff', 'tram_diff', 'combo_diff', 'expected_additive', 'deviation']].to_string(index=False))
    
    add_file = RESULTS_DIR / 'gse288795_additivity_test.csv'
    add_df.to_csv(add_file, index=False)
    print(f"\nSaved additivity test to {add_file}")
    
    # --- Step 9: Figures ---
    print("\n" + "=" * 70)
    print("STEP 9: Generate figures")
    print("=" * 70)
    
    plot_predictions(pred_df)
    plot_drug_differences(diff_df)
    plot_module_heatmap(module_diff_df)
    plot_additivity_test(diff_df)
    
    # --- Step 10: Summary ---
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    print("\nOverall rejuvenation effects (Drug - Control, negative = younger):")
    summary = diff_df.groupby('treatment')['difference'].mean().sort_values()
    print(summary.to_string())
    
    print("\nBy tissue:")
    tissue_summary = diff_df.groupby(['tissue', 'treatment'])['difference'].mean().unstack()
    print(tissue_summary.to_string())
    
    print("\nBy sex:")
    sex_summary = diff_df.groupby(['sex', 'treatment'])['difference'].mean().unstack()
    print(sex_summary.to_string())
    
    print("\nAdditivity conclusion:")
    for _, row in add_df.iterrows():
        direction = "synergistic" if row['is_synergistic'] else ("antagonistic" if row['is_antagonistic'] else "additive")
        print(f"  {row['tissue']} {row['sex']}: Combo={row['combo_diff']:.3f}, Expected={row['expected_additive']:.3f}, Deviation={row['deviation']:.3f} → {direction}")
    
    print("\n" + "=" * 70)
    print("Pipeline complete!")
    print("=" * 70)
    
    return pred_df, diff_df, add_df, module_diff_df


if __name__ == '__main__':
    main()
