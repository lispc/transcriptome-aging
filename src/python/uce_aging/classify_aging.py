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
from sklearn.linear_model import LogisticRegression, ElasticNet, ElasticNetCV
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import accuracy_score, roc_auc_score, mean_absolute_error, r2_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")


def binary_classify(X, y_age, model_type="elasticnet", cv_folds=5):
    """Binary young vs old classification."""
    y = (y_age >= 18).astype(int)  # 1=old, 0=young
    
    if model_type == "elasticnet":
        model = LogisticRegression(penalty="elasticnet", solver="saga",
                                    l1_ratio=0.5, C=0.1, max_iter=5000,
                                    random_state=42)
    elif model_type == "logistic":
        model = LogisticRegression(max_iter=5000, random_state=42)
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Stratified CV
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
    y_proba = cross_val_predict(model, X_scaled, y, cv=cv, method="predict_proba")[:, 1]
    y_pred = (y_proba >= 0.5).astype(int)
    
    acc = accuracy_score(y, y_pred)
    auc = roc_auc_score(y, y_proba)
    
    return {
        "accuracy": acc,
        "auc": auc,
        "n_young": (y == 0).sum(),
        "n_old": (y == 1).sum(),
    }


def regression_predict(X, y_age, model_type="elasticnet", cv_folds=5):
    """Continuous age regression."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    if model_type == "elasticnet":
        model = ElasticNetCV(l1_ratio=0.5, cv=cv_folds, max_iter=5000, random_state=42)
    elif model_type == "xgboost":
        model = GradientBoostingRegressor(n_estimators=200, max_depth=4,
                                           learning_rate=0.1, random_state=42)
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    # Leave-one-age-out CV
    ages = np.unique(y_age)
    y_pred = np.zeros(len(y_age))
    for test_age in ages:
        train_mask = y_age != test_age
        test_mask = y_age == test_age
        if train_mask.sum() < 10:
            continue
        model.fit(X_scaled[train_mask], y_age[train_mask])
        y_pred[test_mask] = model.predict(X_scaled[test_mask])
    
    mae = mean_absolute_error(y_age, y_pred)
    r2 = r2_score(y_age, y_pred)
    
    # Correlation
    pearson_r = np.corrcoef(y_age, y_pred)[0, 1]
    
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
    
    # Task 1: Binary classification
    print("\n=== Binary Classification (young vs old) ===")
    for model_type in ["elasticnet", "logistic"]:
        res = binary_classify(X, y_age, model_type=model_type)
        print(f"  {model_type}: accuracy={res['accuracy']:.3f}, AUC={res['auc']:.3f}")
        for k, v in res.items():
            results[f"binary_{model_type}_{k}"] = v
    
    # Task 2: Regression
    print("\n=== Continuous Age Regression ===")
    for model_type in ["elasticnet", "xgboost"]:
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
