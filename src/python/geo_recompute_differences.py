#!/usr/bin/env python3
"""Recompute Drug-Control differences with proper matching by age/sex."""

import pandas as pd
import numpy as np
from pathlib import Path

RESULTS_DIR = Path('/home/scroll/zzhang/transcriptome-aging/results/geo_longterm')

# Load predictions
pred = pd.read_csv(RESULTS_DIR / 'gse131754_tage_predictions.csv')

# Define control mapping for each intervention
control_map = {
    'ACA': 'Control',
    'CR': 'Control',
    'EST': 'Control',
    'PROT': 'Control',
    'RAP': 'Control',
    'GHRKO': 'GHRKO_Control',
    'MR': 'MR_Control',
    'SNELL': 'Snell_Control',
}

results = []

for intervention, control_name in control_map.items():
    drug_data = pred[pred['intervention'] == intervention]
    control_data = pred[pred['intervention'] == control_name]
    
    # Group by age and sex
    for (age, sex), drug_group in drug_data.groupby(['age', 'sex']):
        ctrl_group = control_data[(control_data['age'] == age) & (control_data['sex'] == sex)]
        
        if len(ctrl_group) == 0:
            print(f"Warning: No control for {intervention} {age} {sex}")
            continue
        
        drug_mean = drug_group['tage_prediction'].mean()
        ctrl_mean = ctrl_group['tage_prediction'].mean()
        diff = drug_mean - ctrl_mean
        
        results.append({
            'intervention': intervention,
            'age': age,
            'sex': sex,
            'n_drug': len(drug_group),
            'n_control': len(ctrl_group),
            'drug_mean': drug_mean,
            'control_mean': ctrl_mean,
            'difference': diff,
            'drug_sem': drug_group['tage_prediction'].sem(),
            'control_sem': ctrl_group['tage_prediction'].sem(),
        })

diff_df = pd.DataFrame(results)

print("=" * 70)
print("GSE131754 Drug - Control Differences (negative = rejuvenation)")
print("=" * 70)
print(diff_df[['intervention', 'age', 'sex', 'difference', 'drug_mean', 'control_mean', 'n_drug', 'n_control']].to_string(index=False))

# Also compute pooled (all ages/sexes combined) for each intervention
print("\n" + "=" * 70)
print("Pooled differences (all ages/sexes combined)")
print("=" * 70)

pooled = []
for intervention, control_name in control_map.items():
    drug_data = pred[pred['intervention'] == intervention]
    control_data = pred[pred['intervention'] == control_name]
    
    # For pooled, match controls by age/sex for each drug sample
    diffs = []
    for _, drug_row in drug_data.iterrows():
        ctrl_match = control_data[
            (control_data['age'] == drug_row['age']) & 
            (control_data['sex'] == drug_row['sex'])
        ]
        if len(ctrl_match) > 0:
            diffs.append(drug_row['tage_prediction'] - ctrl_match['tage_prediction'].mean())
    
    if diffs:
        pooled.append({
            'intervention': intervention,
            'n_samples': len(diffs),
            'mean_diff': np.mean(diffs),
            'sem': np.std(diffs) / np.sqrt(len(diffs)),
            'std': np.std(diffs),
        })

pooled_df = pd.DataFrame(pooled).sort_values('mean_diff')
print(pooled_df.to_string(index=False))

# Save
pooled_df.to_csv(RESULTS_DIR / 'gse131754_pooled_differences.csv', index=False)
diff_df.to_csv(RESULTS_DIR / 'gse131754_drug_control_differences_v2.csv', index=False)
print(f"\nSaved to {RESULTS_DIR}")
