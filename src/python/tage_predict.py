"""
Python inference script for tAge models.
Loads a scikit-learn Pipeline (.pkl) and runs prediction on a gene expression matrix.

Expected input: pandas DataFrame with columns = model feature names (Entrez IDs).
The model pipeline handles imputation, centering, and prediction internally.
"""

import argparse
import joblib
import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=UserWarning)


def load_model(model_path: str):
    """Load a scikit-learn Pipeline from a pickle file."""
    model = joblib.load(model_path)
    return model


def predict(model, expr_df: pd.DataFrame) -> pd.DataFrame:
    """
    Run prediction.
    
    Parameters
    ----------
    model : sklearn Pipeline
    expr_df : pd.DataFrame
        Rows = samples, columns = gene features (must match model.feature_names_in_)
    
    Returns
    -------
    pd.DataFrame with columns ['sample_id', 'prediction']
    """
    # Ensure column order matches model training order
    features = list(model.feature_names_in_)
    missing = set(features) - set(expr_df.columns)
    if missing:
        raise ValueError(f"Missing {len(missing)} features from input: {list(missing)[:10]}...")
    
    X = expr_df[features].values.astype(float)
    y_pred = model.predict(X)
    
    return pd.DataFrame({
        "sample_id": expr_df.index,
        "prediction": y_pred
    })


def construct_test_data(features: list, n_samples: int = 5, seed: int = 42) -> pd.DataFrame:
    """Construct a synthetic expression matrix for testing."""
    rng = np.random.default_rng(seed)
    # Use small random values centered around 0 (relative expression scale)
    data = rng.normal(loc=0.0, scale=0.5, size=(n_samples, len(features)))
    df = pd.DataFrame(data, columns=features)
    df.index = [f"sample_{i+1}" for i in range(n_samples)]
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="tAge Python inference")
    parser.add_argument("--model", required=True, help="Path to .pkl model")
    parser.add_argument("--input", help="Path to CSV expression matrix (samples x genes)")
    parser.add_argument("--output", default="results/python_predictions.csv", help="Output CSV path")
    parser.add_argument("--test", action="store_true", help="Run on synthetic test data")
    args = parser.parse_args()
    
    model = load_model(args.model)
    print(f"Loaded model: {args.model}")
    print(f"Expected features: {len(model.feature_names_in_)}")
    
    if args.test:
        expr_df = construct_test_data(list(model.feature_names_in_))
        expr_df.to_csv("data/test/synthetic_test_matrix.csv")
        print(f"Saved synthetic test data to data/test/synthetic_test_matrix.csv")
    elif args.input:
        expr_df = pd.read_csv(args.input, index_col=0)
    else:
        raise ValueError("Provide --input or use --test")
    
    results = predict(model, expr_df)
    results.to_csv(args.output, index=False)
    print(f"Saved predictions to {args.output}")
    print(results)
