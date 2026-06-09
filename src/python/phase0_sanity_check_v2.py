"""
Phase 0 Sanity Check v2: Matched control analysis.

For each drug sample, find the closest-matched control (same tissue, species, strain, sex, similar age)
and compute control-subtracted expression before prediction.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from tage_predict_official import predict_tAge


def find_matched_controls(drug_ann, control_ann, max_age_diff=30):
    """For each drug sample, find matched controls with same tissue/strain/sex and similar age."""
    matched = []
    for _, drug_row in drug_ann.iterrows():
        mask = (
            (control_ann["Tissue"] == drug_row["Tissue"]) &
            (control_ann["Species"] == drug_row["Species"]) &
            (control_ann["Strain"] == drug_row["Strain"]) &
            (control_ann["Sex"] == drug_row["Sex"])
        )
        candidates = control_ann[mask].copy()
        if len(candidates) == 0:
            # Relax strain match
            mask = (
                (control_ann["Tissue"] == drug_row["Tissue"]) &
                (control_ann["Species"] == drug_row["Species"]) &
                (control_ann["Sex"] == drug_row["Sex"])
            )
            candidates = control_ann[mask].copy()
        
        if len(candidates) == 0:
            matched.append(None)
            continue
        
        candidates["age_diff"] = np.abs(candidates["Age.days"] - drug_row["Age.days"])
        best = candidates.sort_values("age_diff").iloc[0]
        matched.append(best["Sample"])
    
    return matched


def control_subtract(expr_df, control_samples):
    """Subtract control group median from expression matrix."""
    control_expr = expr_df[control_samples]
    control_median = control_expr.median(axis=1)
    return expr_df.sub(control_median, axis=0)


def run_sanity_check(drug_name, keyword, expr, ann):
    print(f"\n=== {drug_name} ===")
    drug_mask = ann["Intervention"].str.contains(keyword, case=False, na=False)
    control_mask = ann["Intervention"] == "Control"
    
    drug_ann = ann[drug_mask].copy()
    control_ann_all = ann[control_mask].copy()
    
    if len(drug_ann) == 0:
        print(f"  No samples found")
        return None
    
    # Find matched controls
    matched_control_samples = find_matched_controls(drug_ann, control_ann_all)
    valid_idx = [i for i, s in enumerate(matched_control_samples) if s is not None]
    
    if len(valid_idx) == 0:
        print(f"  No matched controls found")
        return None
    
    drug_ann = drug_ann.iloc[valid_idx].reset_index(drop=True)
    matched_controls = pd.DataFrame([control_ann_all[control_ann_all["Sample"] == s].iloc[0] 
                                      for s in [matched_control_samples[i] for i in valid_idx]])
    
    drug_samples = drug_ann["Sample"].tolist()
    control_samples = matched_controls["Sample"].tolist()
    
    # Combine for control subtraction
    combined_samples = list(set(drug_samples + control_samples))
    combined_expr = expr[combined_samples]
    
    # Control subtraction
    control_median = expr[control_samples].median(axis=1)
    diff_expr = combined_expr.sub(control_median, axis=0)
    
    # Predict
    sub_expr = diff_expr.T
    sub_expr.columns = sub_expr.columns.map(str)
    meta = pd.DataFrame(index=sub_expr.index)
    meta["group"] = ["drug" if s in drug_samples else "control" for s in sub_expr.index]
    
    model_path = "models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
    preds = predict_tAge(
        model_path=model_path,
        exprs_data_df=sub_expr,
        annotation_data_df=meta,
        species="mouse",
        return_std=False,
        prefix="EN_"
    )
    
    drug_pred = preds[preds["group"] == "drug"]["EN_tAge"]
    control_pred = preds[preds["group"] == "control"]["EN_tAge"]
    
    result = {
        "drug": drug_name,
        "n_drug": len(drug_pred),
        "n_control": len(control_pred),
        "drug_mean": drug_pred.mean(),
        "drug_median": drug_pred.median(),
        "drug_std": drug_pred.std(),
        "control_mean": control_pred.mean(),
        "control_median": control_pred.median(),
        "control_std": control_pred.std(),
        "delta_mean": drug_pred.mean() - control_pred.mean(),
        "delta_median": drug_pred.median() - control_pred.median(),
    }
    
    print(f"  Matched drug samples: {result['n_drug']}")
    print(f"  Matched control samples: {result['n_control']}")
    print(f"  Drug mean tAge: {result['drug_mean']:.3f} ± {result['drug_std']:.3f}")
    print(f"  Control mean tAge: {result['control_mean']:.3f} ± {result['control_std']:.3f}")
    print(f"  Delta (Drug - Control): {result['delta_mean']:.3f}")
    
    if result['delta_mean'] < 0:
        print(f"  ✅ Direction: Rejuvenation (lower tAge)")
    else:
        print(f"  ❌ Direction: Accelerated ageing (higher tAge)")
    
    return result


def main():
    print("Loading rodent data...")
    expr = pd.read_csv("data/raw/Expression_data_relative_rodents_Scaled.csv", index_col=0)
    ann = pd.read_excel("data/raw/Data_annotation_relative_rodents.xlsx")
    print(f"Expression: {expr.shape}, Annotation: {ann.shape}")
    
    drugs = {
        "Rapamycin": "Rapamycin",
        "Acarbose": "Acarbose",
        "Canagliflozin": "Canagliflozin",
    }
    
    results = []
    for name, keyword in drugs.items():
        res = run_sanity_check(name, keyword, expr, ann)
        if res:
            results.append(res)
    
    df = pd.DataFrame(results)
    df = df[["drug", "n_drug", "n_control", "drug_mean", "control_mean", 
             "delta_mean", "delta_median"]]
    df.to_csv("results/phase0_sanity_check_v2.csv", index=False)
    print(f"\n=== Summary ===")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
