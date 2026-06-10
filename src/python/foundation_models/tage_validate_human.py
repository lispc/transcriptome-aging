"""
tAge validation on human AIDA PBMC data.
Test: does tAge predict chronological age in human PBMCs?
"""
import json
import re
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import joblib
from pathlib import Path

DATA_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data")
OUTPUT_DIR = DATA_DIR / "geneformer/tokenized"


def load_mapping():
    """Load ENSG → Entrez mapping."""
    with open(DATA_DIR / "gene_id_mappings/ensembl_to_entrez_geneformer.json") as f:
        return json.load(f)


def parse_age(stage_str):
    """Parse age from AIDA development_stage string."""
    match = re.search(r'(\d+)-year-old', str(stage_str).lower())
    if match:
        return float(match.group(1))
    return np.nan


def compute_pseudobulk_per_donor(adata_path, max_cells_per_donor=1000, n_donors=None):
    """
    Compute pseudobulk per donor by subsampling cells.
    Uses backed mode to handle large files.
    """
    print("Loading AIDA metadata...")
    meta = pd.read_csv(OUTPUT_DIR / "aida_v1_metadata.csv", index_col=0)
    
    # Get donor ages
    donor_ages = meta.groupby('donor_id')['age_numeric'].first()
    print(f"Total donors: {len(donor_ages)}")
    print(f"Age range: {donor_ages.min():.0f} - {donor_ages.max():.0f}")
    
    if n_donors:
        donor_ages = donor_ages.head(n_donors)
        print(f"Using first {n_donors} donors for validation")
    
    donors = donor_ages.index.tolist()
    
    # Open backed
    adata = sc.read_h5ad(adata_path, backed='r')
    gene_names = list(adata.var_names)
    
    pseudobulk_list = []
    donor_age_list = []
    
    for i, donor in enumerate(donors):
        if i % 50 == 0:
            print(f"  Processing donor {i+1}/{len(donors)}: {donor}")
        
        # Get cell indices for this donor
        donor_mask = meta['donor_id'] == donor
        donor_indices = meta[donor_mask].index.tolist()
        
        # Subsample if too many cells
        if len(donor_indices) > max_cells_per_donor:
            np.random.seed(42)
            donor_indices = np.random.choice(donor_indices, max_cells_per_donor, replace=False).tolist()
        
        # Get expression (backed mode)
        # Convert index names to integer positions
        all_obs_names = list(adata.obs_names)
        pos_map = {name: idx for idx, name in enumerate(all_obs_names)}
        positions = [pos_map[idx] for idx in donor_indices if idx in pos_map]
        
        if len(positions) == 0:
            continue
        
        expr = adata.X[positions, :]
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        
        # Pseudobulk: mean across cells
        bulk = np.mean(expr, axis=0)
        
        pseudobulk_list.append(bulk)
        donor_age_list.append(donor_ages[donor])
    
    adata.file.close()
    
    # Create DataFrame
    pseudobulk_df = pd.DataFrame(pseudobulk_list, index=donors, columns=gene_names)
    age_series = pd.Series(donor_age_list, index=donors, name='chronological_age')
    
    return pseudobulk_df, age_series


def map_genes_and_predict_tage(pseudobulk_df, mapping_dict):
    """Map ENSG to Entrez and run tAge prediction."""
    print(f"\nPseudobulk shape: {pseudobulk_df.shape}")
    
    # Map column names
    mapped_cols = {}
    for col in pseudobulk_df.columns:
        if col in mapping_dict:
            mapped_cols[col] = mapping_dict[col]
    
    print(f"Genes mapped to Entrez: {len(mapped_cols)} / {len(pseudobulk_df.columns)}")
    
    # Rename columns
    pseudobulk_mapped = pseudobulk_df.rename(columns=mapped_cols)
    
    # Aggregate duplicate Entrez IDs (take mean)
    pseudobulk_mapped = pseudobulk_mapped.T.groupby(level=0).mean().T
    print(f"After dedup: {pseudobulk_mapped.shape}")
    
    # Load tAge model
    print("\nLoading tAge model...")
    model = joblib.load(DATA_DIR / "models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    
    # tAge model expects specific features
    # We need to reindex to match model features
    # Get feature names from model
    imputer = model.named_steps['imputation']
    # Features are stored in the pipeline's feature names
    # Since we don't have feature_names_in_, we'll use the imputer's statistics length
    n_features = len(imputer.statistics_)
    print(f"tAge model expects {n_features} features")
    
    # For now, we need to know what the feature names are
    # The model was trained on Entrez IDs, so we need the actual list
    # Let's check if there's a way to get them
    
    # Alternative: use the model's SelectKBest to get feature names
    selector = model.named_steps['kbest']
    print(f"SelectKBest selected {selector.get_support().sum()} features")
    
    # We need the original feature names. Since they're not stored in the model,
    # we'll try to infer from the pipeline or use all available mapped genes
    # For validation, we'll use the intersection of available genes and model features
    
    # Load features from a separate file if available
    features_path = DATA_DIR / "models/tage_features.txt"
    if features_path.exists():
        with open(features_path) as f:
            tage_features = [line.strip() for line in f]
    else:
        # Use available mapped genes as features
        # This is an approximation - the model will impute missing features
        tage_features = list(pseudobulk_mapped.columns)
        print(f"No feature list found, using {len(tage_features)} available genes")
    
    # Reindex
    X = pseudobulk_mapped.reindex(columns=tage_features)
    print(f"Input shape after reindex: {X.shape}")
    print(f"Non-null columns: {X.notna().all().sum()}")
    
    # Fill missing with mean
    X = X.fillna(X.mean())
    
    # Predict
    print("Predicting tAge...")
    y_pred = model.predict(X)
    
    # tAge formula: (y - 5.5) * 48 = months
    tAge_months = (y_pred - 5.5) * 48
    
    return tAge_months


def main():
    print("=" * 60)
    print("tAge Validation on Human AIDA PBMC Data")
    print("=" * 60)
    
    # Load mapping
    mapping_dict = load_mapping()
    print(f"Loaded mapping: {len(mapping_dict)} ENSG → Entrez")
    
    # Compute pseudobulk
    pseudobulk_df, true_ages = compute_pseudobulk_per_donor(
        DATA_DIR / "scrna_aging/human/aida_v1.h5ad",
        max_cells_per_donor=500,  # Subsample for speed
        n_donors=100,  # Start with 100 donors for quick validation
    )
    
    # Predict tAge
    tAge_pred = map_genes_and_predict_tage(pseudobulk_df, mapping_dict)
    
    # Evaluate
    from scipy.stats import pearsonr, spearmanr
    
    results = pd.DataFrame({
        'donor_id': pseudobulk_df.index,
        'chronological_age': true_ages.values,
        'tAge_months': tAge_pred,
    })
    
    # Remove any NaN
    results = results.dropna()
    
    print(f"\n{'='*60}")
    print("Results")
    print(f"{'='*60}")
    print(f"Valid donors: {len(results)}")
    
    if len(results) > 10:
        r_pearson, p_pearson = pearsonr(results['chronological_age'], results['tAge_months'])
        r_spearman, p_spearman = spearmanr(results['chronological_age'], results['tAge_months'])
        mae = np.mean(np.abs(results['chronological_age'] - results['tAge_months'] / 12))
        
        print(f"Pearson r: {r_pearson:.3f} (p={p_pearson:.2e})")
        print(f"Spearman ρ: {r_spearman:.3f} (p={p_spearman:.2e})")
        print(f"MAE (years): {mae:.1f}")
        
        # Save results
        results.to_csv(OUTPUT_DIR / "tage_validation_aida.csv", index=False)
        print(f"\nResults saved to {OUTPUT_DIR / 'tage_validation_aida.csv'}")
    else:
        print("Too few valid results for correlation analysis")
    
    print(results.head(10).to_string())


if __name__ == "__main__":
    main()
