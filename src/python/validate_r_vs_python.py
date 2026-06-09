"""
Validation script: compare R tAge predictions vs direct Python predictions.

Steps:
1. Load the scaled_diff matrix produced by R preprocessing.
2. Run prediction using the official tAge Python inference logic.
3. Load R predictions.
4. Compare numerically (should match to machine precision).
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Import the official tAge prediction logic copied from R package
sys.path.insert(0, str(Path(__file__).parent))
from tage_predict_official import predict_tAge


def main():
    # Paths
    scaled_diff_path = Path("data/test/r_scaled_diff_matrix.csv")
    metadata_path = Path("data/test/r_metadata.csv")
    r_predictions_path = Path("data/test/r_predictions.csv")
    model_path = Path("models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
    
    if not scaled_diff_path.exists():
        print(f"ERROR: {scaled_diff_path} not found. Run src/r/run_validation.R first.")
        sys.exit(1)
    
    # Load R-preprocessed matrix (genes x samples)
    scaled_diff = pd.read_csv(scaled_diff_path, index_col=0)
    print(f"Loaded scaled_diff matrix: {scaled_diff.shape}")
    
    # Load metadata
    meta = pd.read_csv(metadata_path, index_col=0)
    
    # Ensure orientation matches predict_tAge expectations
    # The official script expects exprs_data_df with samples as rows, genes as columns
    # If features are in index, it auto-transposes
    exprs_df = scaled_diff.T
    exprs_df.columns = exprs_df.columns.map(str)
    print(f"Expression matrix for prediction (samples x genes): {exprs_df.shape}")
    
    # Run Python prediction using official logic
    py_results = predict_tAge(
        model_path=model_path,
        exprs_data_df=exprs_df,
        annotation_data_df=meta,
        species="mouse",
        return_std=False,
        prefix="EN_"
    )
    
    # Save Python predictions
    py_results.to_csv("data/test/python_predictions.csv")
    print("\nPython predictions saved to data/test/python_predictions.csv")
    print(py_results[["EN_tAge"]].head())
    
    # Compare with R if available
    if r_predictions_path.exists():
        r_results = pd.read_csv(r_predictions_path, index_col=0)
        
        # Align by sample ID (rownames)
        common_samples = py_results.index.intersection(r_results.index)
        if len(common_samples) == 0:
            # Try matching by first column if index doesn't align
            common_samples = py_results.index.intersection(r_results.iloc[:, 0])
        
        if len(common_samples) == 0:
            print("\nWARNING: Could not align R and Python results by sample ID.")
            print("Python index:", list(py_results.index[:5]))
            print("R index:", list(r_results.index[:5]))
            return
        
        py_vals = py_results.loc[common_samples, "EN_tAge"].values
        r_col = [c for c in r_results.columns if "scaled_diff_EN_tAge" in c or "EN_tAge" in c]
        if not r_col:
            print("\nWARNING: No tAge column found in R results.")
            print("R columns:", list(r_results.columns))
            return
        r_vals = r_results.loc[common_samples, r_col[0]].values
        
        diff = py_vals - r_vals
        print(f"\n=== R vs Python Comparison ({len(common_samples)} samples) ===")
        print(f"Max absolute difference: {np.max(np.abs(diff)):.6e}")
        print(f"Mean absolute difference: {np.mean(np.abs(diff)):.6e}")
        print(f"Correlation: {np.corrcoef(py_vals, r_vals)[0,1]:.10f}")
        
        if np.max(np.abs(diff)) < 1e-10:
            print("\n✅ PASS: R and Python predictions are numerically identical.")
        elif np.max(np.abs(diff)) < 1e-6:
            print("\n✅ PASS: R and Python predictions match within floating-point tolerance.")
        else:
            print("\n❌ FAIL: Differences exceed expected tolerance. Investigate.")
            print("Sample-by-sample diff:")
            for i, s in enumerate(common_samples):
                print(f"  {s}: Py={py_vals[i]:.6f}, R={r_vals[i]:.6f}, diff={diff[i]:.6e}")
    else:
        print(f"\n{r_predictions_path} not found yet. Run R validation first.")


if __name__ == "__main__":
    main()
