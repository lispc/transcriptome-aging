"""
Phase 0 Sanity Check: Validate tAge pipeline on known drugs using rodent ITP data.

Uses the pre-processed rodent Scaled expression matrix from Zenodo.
For each drug, compares drug-treated samples vs. matched controls.
Expected directions:
- Rapamycin, Acarbose, Canagliflozin -> decreased mortality tAge (rejuvenation)
- Doxorubicin, Dexamethasone -> increased mortality tAge (pro-ageing)
"""

import argparse
import joblib
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from tage_predict_official import predict_tAge


def load_rodent_data():
    """Load pre-processed rodent expression matrix and annotation."""
    expr_path = Path("data/raw/Expression_data_relative_rodents_Scaled.csv")
    ann_path = Path("data/raw/Data_annotation_relative_rodents.xlsx")
    
    print(f"Loading expression matrix from {expr_path} ...")
    expr = pd.read_csv(expr_path, index_col=0)
    print(f"  Expression matrix shape: {expr.shape}")
    
    print(f"Loading annotation from {ann_path} ...")
    ann = pd.read_excel(ann_path)
    print(f"  Annotation shape: {ann.shape}")
    
    return expr, ann


def run_predictions(expr, sample_ids, model_path):
    """Run tAge predictions on a subset of samples."""
    # Ensure sample IDs exist in expression matrix
    missing = set(sample_ids) - set(expr.columns)
    if missing:
        print(f"  Warning: {len(missing)} samples missing from expression matrix")
        sample_ids = [s for s in sample_ids if s in expr.columns]
    
    sub_expr = expr[sample_ids].T  # samples x genes
    sub_expr.columns = sub_expr.columns.map(str)
    
    meta = pd.DataFrame(index=sub_expr.index)
    
    results = predict_tAge(
        model_path=model_path,
        exprs_data_df=sub_expr,
        annotation_data_df=meta,
        species="mouse",
        return_std=False,
        prefix="EN_"
    )
    return results["EN_tAge"]


def compare_drug_vs_control(expr, ann, drug_keyword, control_keyword="Control"):
    """Find drug and control samples, run predictions, and compare."""
    drug_mask = ann["Intervention"].str.contains(drug_keyword, case=False, na=False)
    control_mask = ann["Intervention"] == control_keyword
    
    # Match by tissue and species to ensure fair comparison
    drug_ann = ann[drug_mask].copy()
    control_ann = ann[control_mask].copy()
    
    if len(drug_ann) == 0:
        print(f"  No samples found for '{drug_keyword}'")
        return None
    
    # Use all matched controls (same tissue/species if possible)
    matched_controls = control_ann.copy()
    
    drug_samples = drug_ann["Sample"].tolist()
    control_samples = matched_controls["Sample"].tolist()
    
    print(f"  Drug samples: {len(drug_samples)}")
    print(f"  Control samples: {len(control_samples)}")
    
    model_path = "models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
    
    drug_pred = run_predictions(expr, drug_samples, model_path)
    control_pred = run_predictions(expr, control_samples, model_path)
    
    return {
        "drug_mean": drug_pred.mean(),
        "drug_std": drug_pred.std(),
        "control_mean": control_pred.mean(),
        "control_std": control_pred.std(),
        "delta_mean": drug_pred.mean() - control_pred.mean(),
        "delta_median": drug_pred.median() - control_pred.median(),
        "n_drug": len(drug_pred),
        "n_control": len(control_pred)
    }


def main():
    parser = argparse.ArgumentParser(description="Phase 0 Sanity Check")
    parser.add_argument("--output", default="results/phase0_sanity_check.csv")
    args = parser.parse_args()
    
    expr, ann = load_rodent_data()
    
    # Define drugs to test
    drugs = {
        "Rapamycin": "Rapamycin",
        "Acarbose": "Acarbose",
        "Canagliflozin": "Canagliflozin",
        # Pro-ageing controls (not in rodent ITP data, but check anyway)
        # "Doxorubicin": "Doxorubicin",
        # "Dexamethasone": "Dexamethasone",
    }
    
    results = []
    for drug_name, keyword in drugs.items():
        print(f"\n=== {drug_name} ===")
        res = compare_drug_vs_control(expr, ann, keyword)
        if res:
            res["drug"] = drug_name
            results.append(res)
            print(f"  Drug mean tAge: {res['drug_mean']:.3f} ± {res['drug_std']:.3f}")
            print(f"  Control mean tAge: {res['control_mean']:.3f} ± {res['control_std']:.3f}")
            print(f"  Delta (Drug - Control): {res['delta_mean']:.3f}")
            if res['delta_mean'] < 0:
                print(f"  ✅ Direction: Rejuvenation (lower tAge)")
            else:
                print(f"  ❌ Direction: Accelerated ageing (higher tAge)")
    
    if results:
        df = pd.DataFrame(results)
        df = df[["drug", "n_drug", "n_control", "drug_mean", "control_mean", "delta_mean", "delta_median"]]
        df.to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")
        print(df.to_string(index=False))
    else:
        print("No results generated.")


if __name__ == "__main__":
    main()
