#!/usr/bin/env python3
"""
GEO Long-term Drug Analysis Pipeline
=====================================
Process GSE131754 and GSE299228 long-term drug treatment RNA-seq data,
run tAge prediction, and compare with LINCS short-term signatures.

Steps:
1. Load GEO counts matrices
2. Map Ensembl IDs to Entrez Gene IDs (using mygene.info)
3. TMM normalization + logCPM (via edgeR in R)
4. tAge prediction for each sample
5. Compute Drug - Control differences
6. Compare with LINCS module fingerprints

Author: Analysis pipeline
Date: 2026-06-09
"""

import os
import sys
import gzip
import pickle
import subprocess
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings('ignore')

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / 'data' / 'geo'
MODEL_PATH = PROJECT_ROOT / 'models' / 'EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl'
RESULTS_DIR = PROJECT_ROOT / 'results' / 'geo_longterm'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# STEP 1: Load and prepare GSE131754 counts
# ============================================================================

def load_gse131754_counts():
    """Load GSE131754 featureCounts matrix."""
    filepath = DATA_DIR / 'GSE131754_Interventions_assigned_reads.txt.gz'
    print(f"Loading GSE131754 from {filepath}")
    
    # Read tab-delimited, first column is gene ID
    df = pd.read_csv(filepath, sep='\t', compression='gzip', 
                     header=0, index_col=0, low_memory=False)
    
    # Remove quotes from column names (GEO format artifact)
    df.columns = [c.strip('"') for c in df.columns]
    df.index.name = 'gene_id'
    
    print(f"  Shape: {df.shape}")
    print(f"  Genes: {df.shape[0]:,}")
    print(f"  Samples: {df.shape[1]}")
    print(f"  Sample groups: {df.columns.tolist()[:5]}...")
    
    return df


def parse_sample_metadata_gse131754(columns):
    """Parse GSE131754 column names into metadata."""
    meta = []
    for col in columns:
        parts = col.split('_')
        # Format: INTERVENTION_AGE_SEX_REP
        # e.g., ACA_12m_F_1, CON_6m_M_2, GHRCON_5m_M_1
        if col.startswith('CON_'):
            intervention = 'Control'
        elif col.startswith('GHRCON_'):
            intervention = 'GHRKO_Control'
        elif col.startswith('MRCON_'):
            intervention = 'MR_Control'
        elif col.startswith('SNELLCON_'):
            intervention = 'Snell_Control'
        else:
            intervention = parts[0]
        
        age = parts[1] if len(parts) > 1 else 'unknown'
        sex = parts[2] if len(parts) > 2 else 'unknown'
        rep = parts[3] if len(parts) > 3 else '1'
        
        meta.append({
            'sample_id': col,
            'intervention': intervention,
            'age': age,
            'sex': sex,
            'replicate': rep
        })
    
    return pd.DataFrame(meta)


# ============================================================================
# STEP 2: Ensembl ID → Entrez Gene ID mapping
# ============================================================================

def map_ensembl_to_entrez(ensembl_ids, species='mouse', batch_size=1000):
    """Map Ensembl gene IDs to Entrez Gene IDs using mygene.info."""
    try:
        import mygene
    except ImportError:
        print("ERROR: mygene package not installed. Run: pip install mygene")
        sys.exit(1)
    
    print(f"\nMapping {len(ensembl_ids)} Ensembl IDs to Entrez...")
    mg = mygene.MyGeneInfo()
    
    # Remove version suffixes (e.g., ENSMUSG00000102693.1 → ENSMUSG00000102693)
    clean_ids = [gid.split('.')[0] for gid in ensembl_ids]
    
    mapping = {}
    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i+batch_size]
        print(f"  Batch {i//batch_size + 1}/{(len(clean_ids)-1)//batch_size + 1}: {len(batch)} IDs")
        
        result = mg.querymany(
            batch,
            scopes='ensembl.gene',
            fields='entrezgene',
            species=species,
            verbose=False,
            as_dataframe=True
        )
        
        for idx, row in result.iterrows():
            ensembl_id = idx
            if 'entrezgene' in row and pd.notna(row['entrezgene']):
                mapping[ensembl_id] = int(row['entrezgene'])
    
    print(f"  Mapped: {len(mapping)}/{len(ensembl_ids)} ({100*len(mapping)/len(ensembl_ids):.1f}%)")
    
    return mapping


def load_or_create_mapping(counts_df, cache_file=None):
    """Load mapping from cache or create new one."""
    if cache_file is None:
        cache_file = DATA_DIR / 'ensembl_to_entrez_mapping_gse131754.pkl'
    
    if os.path.exists(cache_file):
        print(f"Loading cached mapping from {cache_file}")
        with open(cache_file, 'rb') as f:
            mapping = pickle.load(f)
        print(f"  Cached mapping: {len(mapping)} genes")
        return mapping
    
    mapping = map_ensembl_to_entrez(counts_df.index.tolist())
    
    with open(cache_file, 'wb') as f:
        pickle.dump(mapping, f)
    print(f"Saved mapping to {cache_file}")
    
    return mapping


# ============================================================================
# STEP 3: Prepare data for edgeR TMM normalization (via R)
# ============================================================================

def prepare_edger_input(counts_df, mapping, output_prefix='gse131754'):
    """
    Prepare counts matrix for edgeR TMM normalization in R.
    
    - Map Ensembl IDs to Entrez IDs
    - Aggregate duplicate Entrez IDs (sum counts)
    - Save as TSV for R processing
    """
    print(f"\nPreparing edgeR input...")
    
    # Create copy with mapped IDs
    df = counts_df.copy()
    df['entrez'] = df.index.map(lambda x: mapping.get(x.split('.')[0], None))
    
    # Remove unmapped genes
    df_mapped = df.dropna(subset=['entrez']).copy()
    df_mapped['entrez'] = df_mapped['entrez'].astype(int)
    
    print(f"  Mapped genes: {len(df_mapped)}")
    
    # Aggregate duplicate Entrez IDs (sum counts)
    df_agg = df_mapped.groupby('entrez').sum()
    
    # Drop the 'entrez' column if it exists after groupby
    if 'entrez' in df_agg.columns:
        df_agg = df_agg.drop(columns=['entrez'])
    
    print(f"  After aggregation: {len(df_agg)} unique Entrez IDs")
    
    # Save for R
    output_counts = RESULTS_DIR / f'{output_prefix}_counts_for_edger.tsv'
    df_agg.to_csv(output_counts, sep='\t')
    print(f"  Saved to {output_counts}")
    
    return df_agg, output_counts


def run_edger_tmm(counts_file, output_prefix='gse131754'):
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
    """Load the pre-trained tAge model using joblib (sklearn pipeline)."""
    print(f"\nLoading tAge model from {model_path}")
    
    import joblib
    model = joblib.load(model_path)
    
    print(f"  Model type: {type(model)}")
    print(f"  Pipeline steps: {list(model.named_steps.keys())}")
    
    # Extract components
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
    """
    Align logCPM matrix with model features.
    
    Steps:
    1. Match model features (Entrez IDs as strings)
    2. Impute missing genes with training set medians (from imputer)
    3. Center using scaler mean (x - mean)
    """
    feature_names = model_info['feature_names']
    impute_stats = model_info['impute_stats']
    center_mean = model_info['center_mean']
    
    print(f"\nAligning genes with model features...")
    print(f"  Model expects: {len(feature_names)} features")
    print(f"  Data has: {len(logcpm_df)} genes")
    
    # Ensure logcpm index is string
    logcpm_df.index = logcpm_df.index.astype(str)
    
    # Find intersection
    common_genes = list(logcpm_df.index.intersection(feature_names))
    print(f"  Common genes: {len(common_genes)}")
    
    # Create aligned matrix (samples x features)
    n_samples = logcpm_df.shape[1]
    n_features = len(feature_names)
    X = np.zeros((n_samples, n_features), dtype=float)
    
    # Build feature name to index mapping
    feat_to_idx = {name: i for i, name in enumerate(feature_names)}
    
    # Fill with available data
    for gene in common_genes:
        col_idx = feat_to_idx[gene]
        X[:, col_idx] = logcpm_df.loc[gene].values.astype(float)
    
    # Impute missing genes with training medians
    missing_count = 0
    for i, gene in enumerate(feature_names):
        if gene not in common_genes:
            X[:, i] = impute_stats[i]
            missing_count += 1
    
    print(f"  Missing genes imputed: {missing_count}")
    
    # Center (scaler step: x - mean)
    X_centered = X - center_mean[np.newaxis, :]
    
    return X_centered, feature_names


def tage_predict(logcpm_df, model_info):
    """
    Run tAge prediction on logCPM matrix.
    
    Preprocessing: impute → center → predict
    Returns per-sample predictions.
    """
    coef = model_info['coef']
    intercept = model_info['intercept']
    
    # Align and preprocess
    X_centered, feature_names = align_and_impute(logcpm_df, model_info)
    
    # Prediction: y = intercept + dot(coef, x_centered)
    predictions = intercept + np.dot(X_centered, coef)
    
    # Per-gene contributions for interpretability
    contributions = X_centered * coef[np.newaxis, :]
    
    return predictions, contributions, feature_names


# ============================================================================
# STEP 5: Compute Drug - Control differences
# ============================================================================

def compute_drug_control_differences(predictions, metadata_df):
    """Compute Drug - Control differences for each intervention."""
    results = []
    
    # Get control samples
    control_samples = metadata_df[
        metadata_df['intervention'].str.contains('Control', case=False)
    ]['sample_id'].tolist()
    
    control_preds = {s: predictions[s] for s in control_samples if s in predictions}
    
    # Group by intervention
    for intervention in metadata_df['intervention'].unique():
        if 'Control' in intervention:
            continue
        
        drug_samples = metadata_df[
            metadata_df['intervention'] == intervention
        ]['sample_id'].tolist()
        
        drug_preds = [predictions[s] for s in drug_samples if s in predictions]
        
        if len(drug_preds) == 0:
            continue
        
        # Find matching controls
        # Strategy: match by age and sex
        drug_meta = metadata_df[metadata_df['intervention'] == intervention].iloc[0]
        
        # Find controls with same age/sex pattern
        if intervention in ['MR', 'GHRKO', 'Snell']:
            # Genetic models: match by age and sex
            matching_controls = metadata_df[
                (metadata_df['intervention'].str.contains('Control')) &
                (metadata_df['age'] == drug_meta['age']) &
                (metadata_df['sex'] == drug_meta['sex'])
            ]
        else:
            # Drug interventions: CON with same age and sex
            matching_controls = metadata_df[
                (metadata_df['intervention'] == 'Control') &
                (metadata_df['age'] == drug_meta['age']) &
                (metadata_df['sex'] == drug_meta['sex'])
            ]
        
        control_ids = matching_controls['sample_id'].tolist()
        matching_control_preds = [predictions[s] for s in control_ids if s in predictions]
        
        if len(matching_control_preds) == 0:
            print(f"  Warning: No matching controls for {intervention}")
            continue
        
        drug_mean = np.mean(drug_preds)
        control_mean = np.mean(matching_control_preds)
        diff = drug_mean - control_mean
        
        results.append({
            'intervention': intervention,
            'age': drug_meta['age'],
            'sex': drug_meta['sex'],
            'n_drug': len(drug_preds),
            'n_control': len(matching_control_preds),
            'drug_mean': drug_mean,
            'control_mean': control_mean,
            'difference': diff,
            'drug_sem': np.std(drug_preds) / np.sqrt(len(drug_preds)),
            'control_sem': np.std(matching_control_preds) / np.sqrt(len(matching_control_preds))
        })
    
    return pd.DataFrame(results)


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("GEO Long-term Drug Analysis Pipeline")
    print("=" * 70)
    
    # --- Step 1: Load data ---
    print("\n" + "=" * 70)
    print("STEP 1: Load GSE131754 counts")
    print("=" * 70)
    counts_df = load_gse131754_counts()
    
    metadata = parse_sample_metadata_gse131754(counts_df.columns)
    print(f"\nSample metadata:")
    print(metadata.groupby(['intervention', 'age', 'sex']).size().to_string())
    
    # --- Step 2: ID mapping ---
    print("\n" + "=" * 70)
    print("STEP 2: Map Ensembl IDs to Entrez IDs")
    print("=" * 70)
    mapping = load_or_create_mapping(counts_df)
    
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
    
    # Load logCPM
    print(f"\nLoading logCPM from {logcpm_file}")
    logcpm_df = pd.read_csv(logcpm_file, sep='\t', index_col=0)
    print(f"  Shape: {logcpm_df.shape}")
    
    # --- Step 5: tAge prediction ---
    print("\n" + "=" * 70)
    print("STEP 5: tAge prediction")
    print("=" * 70)
    model_info = load_tage_model(MODEL_PATH)
    
    predictions, contributions, aligned = tage_predict(logcpm_df, model_info)
    
    # Create prediction dataframe
    pred_df = pd.DataFrame({
        'sample_id': logcpm_df.columns,
        'tage_prediction': predictions
    })
    pred_df = pred_df.merge(metadata, on='sample_id', how='left')
    
    print(f"\nPredictions summary:")
    print(pred_df.groupby(['intervention', 'age', 'sex'])['tage_prediction'].mean().to_string())
    
    # Save predictions
    pred_file = RESULTS_DIR / 'gse131754_tage_predictions.csv'
    pred_df.to_csv(pred_file, index=False)
    print(f"\nSaved predictions to {pred_file}")
    
    # --- Step 6: Drug - Control differences ---
    print("\n" + "=" * 70)
    print("STEP 6: Drug - Control differences")
    print("=" * 70)
    diff_df = compute_drug_control_differences(
        dict(zip(logcpm_df.columns, predictions)),
        metadata
    )
    
    print(f"\nDrug - Control differences (negative = rejuvenation):")
    print(diff_df[['intervention', 'age', 'sex', 'difference', 'n_drug', 'n_control']].to_string(index=False))
    
    diff_file = RESULTS_DIR / 'gse131754_drug_control_differences.csv'
    diff_df.to_csv(diff_file, index=False)
    print(f"\nSaved differences to {diff_file}")
    
    # --- Step 7: Module decomposition (if module info available) ---
    # This will be done in a separate step
    
    print("\n" + "=" * 70)
    print("Pipeline complete!")
    print("=" * 70)
    
    return pred_df, diff_df


if __name__ == '__main__':
    main()
