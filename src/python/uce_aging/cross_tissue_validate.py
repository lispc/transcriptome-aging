"""
Cross-tissue validation: train classifier on Spleen, test on Lung (and vice versa).
"""
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler


def load_data(path):
    adata = sc.read_h5ad(path)
    X = adata.obsm["X_uce"]
    y_age = adata.obs["age"].astype(str).str.replace("m", "").astype(int).values
    y = (y_age >= 18).astype(int)
    return X, y, y_age


def main():
    print("Loading Spleen...")
    X_spleen, y_spleen, _ = load_data("results/uce_spleen/tms_spleen_uce_adata.h5ad")
    print(f"  Spleen: {len(X_spleen)} cells")
    
    print("Loading Lung...")
    X_lung, y_lung, _ = load_data("results/uce_lung/tms_lung_uce_adata.h5ad")
    print(f"  Lung: {len(X_lung)} cells")
    
    scaler_s = StandardScaler().fit(X_spleen)
    scaler_l = StandardScaler().fit(X_lung)
    
    # Train on Spleen, test on Lung
    print("\n=== Train on Spleen, Test on Lung ===")
    model = LogisticRegression(max_iter=5000, random_state=42)
    model.fit(scaler_s.transform(X_spleen), y_spleen)
    proba = model.predict_proba(scaler_l.transform(X_lung))[:, 1]
    acc = accuracy_score(y_lung, proba >= 0.5)
    auc = roc_auc_score(y_lung, proba)
    print(f"  Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    
    # Train on Lung, test on Spleen
    print("\n=== Train on Lung, Test on Spleen ===")
    model = LogisticRegression(max_iter=5000, random_state=42)
    model.fit(scaler_l.transform(X_lung), y_lung)
    proba = model.predict_proba(scaler_s.transform(X_spleen))[:, 1]
    acc = accuracy_score(y_spleen, proba >= 0.5)
    auc = roc_auc_score(y_spleen, proba)
    print(f"  Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    
    # Same-tissue baselines
    print("\n=== Same-tissue baselines (for reference) ===")
    from sklearn.model_selection import cross_val_predict, StratifiedKFold
    cv = StratifiedKFold(5, shuffle=True, random_state=42)
    
    proba = cross_val_predict(LogisticRegression(max_iter=5000, random_state=42),
                               scaler_s.transform(X_spleen), y_spleen, cv=cv, method="predict_proba")[:, 1]
    print(f"  Spleen CV: AUC={roc_auc_score(y_spleen, proba):.3f}")
    
    proba = cross_val_predict(LogisticRegression(max_iter=5000, random_state=42),
                               scaler_l.transform(X_lung), y_lung, cv=cv, method="predict_proba")[:, 1]
    print(f"  Lung CV: AUC={roc_auc_score(y_lung, proba):.3f}")


if __name__ == "__main__":
    main()
