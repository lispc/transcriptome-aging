"""
Cross-tissue aging validation: train on one tissue, test on another.
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score


def cross_tissue_validate(train_path, test_path, model_type="elasticnet"):
    adata_train = sc.read_h5ad(train_path)
    adata_test = sc.read_h5ad(test_path)
    
    X_train = adata_train.obsm["X_uce"]
    y_train = (adata_train.obs["age"].astype(str).str.replace("m", "").astype(int) >= 18).astype(int).values
    
    X_test = adata_test.obsm["X_uce"]
    y_test = (adata_test.obs["age"].astype(str).str.replace("m", "").astype(int) >= 18).astype(int).values
    
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    
    if model_type == "elasticnet":
        model = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.1, max_iter=5000)
    else:
        model = LogisticRegression(max_iter=5000)
    
    model.fit(X_train_s, y_train)
    y_proba = model.predict_proba(X_test_s)[:, 1]
    y_pred = (y_proba >= 0.5).astype(int)
    
    return {
        "train_tissue": adata_train.obs["tissue"].iloc[0],
        "test_tissue": adata_test.obs["tissue"].iloc[0],
        "accuracy": accuracy_score(y_test, y_pred),
        "auc": roc_auc_score(y_test, y_proba),
        "n_train": len(X_train),
        "n_test": len(X_test),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    
    res = cross_tissue_validate(args.train, args.test)
    print(res)
    if args.output:
        pd.DataFrame([res]).to_csv(args.output, index=False)
