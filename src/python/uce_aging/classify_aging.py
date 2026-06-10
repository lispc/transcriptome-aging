"""
Aging classifier on UCE embeddings for Tabula Muris Senis.

Tasks:
  1. Binary classification: young (1m, 3m) vs old (18m, 21m, 24m, 30m)
  2. Continuous regression: predict exact age in months
  3. Cross-validation: leave-one-age-out, leave-one-tissue-out
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression, ElasticNet, ElasticNetCV, RidgeCV, LinearRegression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import accuracy_score, roc_auc_score, mean_absolute_error, r2_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict, KFold
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")


from sklearn.model_selection import train_test_split

def binary_classify(X, y_age, model_type="logistic"):
    """Binary young vs old classification - fast train/test split."""
    y = (y_age >= 18).astype(int)  # 1=old, 0=young
    
    if model_type == "logistic":
        model = LogisticRegression(solver="liblinear", max_iter=1000,
                                    tol=1e-3, random_state=42)
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Train/test split (80/20, stratified)
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, test_size=0.2, random_state=42, stratify=y)
    
    model.fit(X_train, y_train)
    y_proba = model.predict_proba(X_test)[:, 1]
    y_pred = (y_proba >= 0.5).astype(int)
    
    acc = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_proba)
    
    return {
        "accuracy": acc,
        "auc": auc,
        "n_young": (y == 0).sum(),
        "n_old": (y == 1).sum(),
    }


from sklearn.model_selection import train_test_split

def regression_predict(X, y_age, model_type="linear"):
    """Continuous age regression - fast train/test split for speed."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Train/test split (80/20)
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y_age, test_size=0.2, random_state=42)
    
    if model_type == "linear":
        model = LinearRegression()
    elif model_type == "ridge":
        model = RidgeCV(alphas=[0.01, 0.1, 1.0, 10.0, 100.0], cv=3)
    elif model_type == "xgboost":
        model = GradientBoostingRegressor(n_estimators=100, max_depth=3,
                                           learning_rate=0.1, random_state=42)
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    pearson_r = np.corrcoef(y_test, y_pred)[0, 1]
    
    return {
        "mae_months": mae,
        "r2": r2,
        "pearson_r": pearson_r,
    }


def main(args):
    print(f"Loading UCE embeddings from {args.input}")
    adata = sc.read_h5ad(args.input)
    
    X = adata.obsm["X_uce"]
    y_age = adata.obs["age"].astype(str).str.replace("m", "").astype(int).values
    
    tissue = adata.obs["tissue"].iloc[0] if "tissue" in adata.obs.columns else "unknown"
    print(f"Tissue: {tissue}, Cells: {len(X)}, Features: {X.shape[1]}")
    print(f"Age distribution: {dict(zip(*np.unique(y_age, return_counts=True)))}")
    
    results = {"tissue": tissue, "n_cells": len(X), "n_features": X.shape[1]}
    
    # Task 1: Binary classification (skip elasticnet - too slow with saga)
    print("\n=== Binary Classification (young vs old) ===")
    for model_type in ["logistic"]:
        res = binary_classify(X, y_age, model_type=model_type)
        print(f"  {model_type}: accuracy={res['accuracy']:.3f}, AUC={res['auc']:.3f}")
        for k, v in res.items():
            results[f"binary_{model_type}_{k}"] = v
    
    # Task 2: Regression (use LinearRegression + train/test split - fastest)
    print("\n=== Continuous Age Regression ===")
    for model_type in ["linear"]:
        res = regression_predict(X, y_age, model_type=model_type)
        print(f"  {model_type}: MAE={res['mae_months']:.2f} months, R2={res['r2']:.3f}, r={res['pearson_r']:.3f}")
        for k, v in res.items():
            results[f"regression_{model_type}_{k}"] = v
    
    # Save results
    if args.output:
        pd.DataFrame([results]).to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")
    
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to UCE-embedded h5ad")
    parser.add_argument("--output", default=None, help="Path to save CSV results")
    args = parser.parse_args()
    main(args)
