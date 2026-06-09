#!/usr/bin/env python3
"""
Unified GEO Pipeline for Direction 3: ARCHS4 Backup via GEO Expansion.
Processes multiple mouse liver GEO datasets through the tAge prediction model.
"""
import warnings
import pickle
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import joblib
import h5py
from sklearn.exceptions import InconsistentVersionWarning

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings("ignore", category=UserWarning)

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
MODEL_PATH = PROJECT_ROOT / "models" / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
COUNTS_DIR = PROJECT_ROOT / "data" / "geo_expansion" / "counts"
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------
def load_model_info():
    model = joblib.load(MODEL_PATH)
    imputer = model.named_steps["imputation"]
    scaler = model.named_steps["scaler"]
    estimator = model.named_steps["estimator"]
    feature_names = list(model.feature_names_in_)
    return {
        "model": model,
        "feature_names": feature_names,
        "coef": estimator.coef_,
        "intercept": float(estimator.intercept_),
        "impute_stats": imputer.statistics_,
        "center_mean": scaler.mean_,
    }


# ---------------------------------------------------------------------------
# Gene mapping helpers
# ---------------------------------------------------------------------------
def get_mygene():
    try:
        import mygene
    except ImportError:
        raise ImportError("mygene is required for gene mapping")
    return mygene.MyGeneInfo()


def map_ensembl_to_entrez(ensembl_ids: List[str], species="mouse", batch_size=1000) -> Dict[str, int]:
    mg = get_mygene()
    clean = [g.split(".")[0] for g in ensembl_ids]
    mapping = {}
    for i in range(0, len(clean), batch_size):
        batch = clean[i:i + batch_size]
        result = mg.querymany(batch, scopes="ensembl.gene", fields="entrezgene",
                              species=species, verbose=False, as_dataframe=True)
        for idx, row in result.iterrows():
            if "entrezgene" in row and pd.notna(row["entrezgene"]):
                mapping[idx] = int(row["entrezgene"])
    return mapping


def map_symbols_to_entrez(symbols: List[str], species="mouse", batch_size=1000) -> Dict[str, int]:
    mg = get_mygene()
    mapping = {}
    for i in range(0, len(symbols), batch_size):
        batch = symbols[i:i + batch_size]
        result = mg.querymany(batch, scopes="symbol,alias", fields="entrezgene",
                              species=species, verbose=False, as_dataframe=True)
        for idx, row in result.iterrows():
            if "entrezgene" in row and pd.notna(row["entrezgene"]):
                mapping[idx] = int(row["entrezgene"])
    return mapping


def load_or_build_mapping(cache_path: Path, ids: List[str], map_type: str) -> Dict[str, int]:
    if cache_path.exists():
        with open(cache_path, "rb") as f:
            mapping = pickle.load(f)
        print(f"  Loaded cached mapping: {len(mapping)} / {len(ids)}")
        return mapping
    if map_type == "ensembl":
        mapping = map_ensembl_to_entrez(ids)
    else:
        mapping = map_symbols_to_entrez(ids)
    with open(cache_path, "wb") as f:
        pickle.dump(mapping, f)
    print(f"  Built mapping: {len(mapping)} / {len(ids)}")
    return mapping


# ---------------------------------------------------------------------------
# Prediction helper
# ---------------------------------------------------------------------------
def predict_from_logexpr(log_expr: pd.DataFrame, model_info: dict) -> np.ndarray:
    """log_expr: genes (Entrez str) x samples. Returns raw predictions (before species adj)."""
    feature_names = model_info["feature_names"]
    n_samples = log_expr.shape[1]
    n_features = len(feature_names)
    X = np.zeros((n_samples, n_features), dtype=float)
    f2i = {name: i for i, name in enumerate(feature_names)}

    common = 0
    for gene in log_expr.index:
        if gene in f2i:
            X[:, f2i[gene]] = log_expr.loc[gene].values
            common += 1
    print(f"  Common genes: {common}/{n_features}")

    for i in range(n_features):
        if np.all(X[:, i] == 0):
            X[:, i] = model_info["impute_stats"][i]

    X_centered = X - model_info["center_mean"][np.newaxis, :]
    preds = model_info["intercept"] + np.dot(X_centered, model_info["coef"])
    return preds


# ---------------------------------------------------------------------------
# Dataset processors
# ---------------------------------------------------------------------------

def process_gse282210(model_info: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Young/Old/Very Old liver partial hepatectomy."""
    print("\n=== GSE282210 ===")
    import gzip
    path = COUNTS_DIR / "GSE282210_processed_data.csv.gz"
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, index_col=0)
    print(f"  Raw shape: {df.shape}")

    # Map Ensembl -> Entrez
    genes = df.index.astype(str).tolist()
    mapping = load_or_build_mapping(RESULTS_DIR / "gse282210_ensembl_to_entrez.pkl", genes, "ensembl")
    entrez_list = []
    expr_list = []
    for g in genes:
        cg = g.split(".")[0]
        if cg in mapping:
            entrez_list.append(str(mapping[cg]))
            expr_list.append(df.loc[g].values)
    expr = pd.DataFrame(np.array(expr_list), columns=df.columns, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    log_expr = np.log2(expr + 1)
    preds = predict_from_logexpr(log_expr, model_info)
    preds_adj = preds * 48.0

    # Metadata from sample names
    meta_rows = []
    for col in df.columns:
        parts = col.split("_")
        age_group = parts[0]  # Y, O, V
        sex = "Female" if "F" in age_group else "Male"
        age_group_clean = age_group.replace("F", "").replace("M", "")
        condition = parts[1]  # Ctrl, 52h, etc.
        age_months = {"Y": 3, "O": 14, "V": 23}.get(age_group_clean, np.nan)
        # Handle VOF/VOM where age_group is already V but replace removed nothing
        if age_group.startswith("V"):
            age_group_clean = "V"
            age_months = 23
        meta_rows.append({
            "sample_id": col,
            "dataset": "GSE282210",
            "age_group": age_group_clean,
            "age_months": age_months,
            "sex": sex,
            "condition": condition,
            "tissue": "liver",
            "treatment": "PH" if condition != "Ctrl" else "Control",
        })
    meta = pd.DataFrame(meta_rows)
    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    return meta, expr


def process_gse283201(model_info: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Young vs Aged liver, control vs ethanol."""
    print("\n=== GSE283201 ===")
    path = COUNTS_DIR / "GSE283201_Read_counts.xlsx"
    df = pd.read_excel(path)
    print(f"  Raw shape: {df.shape}")

    # Columns: EMSM_ID, mgi_symbol, then samples
    genes = df["mgi_symbol"].astype(str).tolist()
    sample_cols = [c for c in df.columns if c not in ["EMSM_ID", "mgi_symbol"]]
    mapping = load_or_build_mapping(RESULTS_DIR / "gse283201_symbol_to_entrez.pkl", genes, "symbol")
    entrez_list = []
    expr_list = []
    for i, g in enumerate(genes):
        if g in mapping:
            entrez_list.append(str(mapping[g]))
            expr_list.append(df.iloc[i][sample_cols].values.astype(float))
    expr = pd.DataFrame(np.array(expr_list), columns=sample_cols, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    log_expr = np.log2(expr + 1)
    preds = predict_from_logexpr(log_expr, model_info)
    preds_adj = preds * 48.0

    meta_rows = []
    for col in sample_cols:
        parts = col.split("-")
        tissue = {"H": "hippocampus", "L": "liver"}.get(parts[0], parts[0])
        age_group = parts[1]  # Y or A
        diet = parts[2][:3]  # CDM or EDE
        age_months = {"Y": 3, "A": 22}.get(age_group, np.nan)
        meta_rows.append({
            "sample_id": col,
            "dataset": "GSE283201",
            "age_group": age_group,
            "age_months": age_months,
            "sex": "Unknown",
            "condition": diet,
            "tissue": tissue,
            "treatment": "EtOH" if diet == "EDE" else "Control",
        })
    meta = pd.DataFrame(meta_rows)
    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    return meta, expr


def process_gse253612(model_info: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Aging keto diet: 26-month liver only (no 12-month liver in file)."""
    print("\n=== GSE253612 ===")
    import gzip
    path = COUNTS_DIR / "GSE253612_Newman_2020_Mouse_RNAseq_counts_all.star.tsv.gz"
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, sep="\t", index_col=0)
    print(f"  Raw shape: {df.shape}")

    sample_cols = [c for c in df.columns]
    # Keep only liver samples
    liver_cols = [c for c in sample_cols if c.startswith("Lvr_")]
    df = df[liver_cols]
    print(f"  Liver samples: {len(liver_cols)}")

    genes = df.index.astype(str).tolist()
    mapping = load_or_build_mapping(RESULTS_DIR / "gse253612_ensembl_to_entrez.pkl", genes, "ensembl")
    entrez_list = []
    expr_list = []
    for g in genes:
        cg = g.split(".")[0]
        if cg in mapping:
            entrez_list.append(str(mapping[cg]))
            expr_list.append(df.loc[g].values)
    expr = pd.DataFrame(np.array(expr_list), columns=liver_cols, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    log_expr = np.log2(expr + 1)
    preds = predict_from_logexpr(log_expr, model_info)
    preds_adj = preds * 48.0

    meta_rows = []
    for col in liver_cols:
        parts = col.split("_")
        age_months = int(parts[1])  # 26
        diet = parts[2]  # CD or KD
        sex = "Unknown"
        meta_rows.append({
            "sample_id": col,
            "dataset": "GSE253612",
            "age_group": "Old" if age_months >= 20 else "Young",
            "age_months": age_months,
            "sex": sex,
            "condition": diet,
            "tissue": "liver",
            "treatment": "Keto" if diet == "KD" else "Control",
        })
    meta = pd.DataFrame(meta_rows)
    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    return meta, expr


def process_gse305103(model_info: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Corylin treatment in 50-week (young) and 130-week (old) mice."""
    print("\n=== GSE305103 ===")
    path = COUNTS_DIR / "GSE305103_Fulltable_Liver_counts_TPM.xlsx"
    df = pd.read_excel(path)
    print(f"  Raw shape: {df.shape}")

    # Use Counts_ columns for liver
    count_cols = [c for c in df.columns if c.startswith("Counts_") and "_L(" in c]
    print(f"  Liver count columns: {len(count_cols)}")

    # Already has EntrezID column! Convert float NaN -> int -> str cleanly
    entrez_ids = df["EntrezID"].apply(lambda x: str(int(x)) if pd.notna(x) else "nan").tolist()
    expr = pd.DataFrame(df[count_cols].values.astype(float), columns=count_cols, index=entrez_ids)
    expr = expr[expr.index != "nan"]
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    log_expr = np.log2(expr + 1)
    preds = predict_from_logexpr(log_expr, model_info)
    preds_adj = preds * 48.0

    meta_rows = []
    for col in count_cols:
        # Parse: Counts_OC_F_B1_L(OC_F_L)  -> OC_F_B1
        core = col.replace("Counts_", "").split("(")[0]
        parts = core.split("_")
        group = parts[0]  # OC, OSD, YSD
        sex_code = parts[1] if len(parts) > 1 else ""
        sex = "Female" if sex_code == "F" else "Male"
        # Map group to age/treatment
        if group.startswith("Y"):
            age_months = 12
            treatment = "Control"
        elif group.startswith("O") and len(group) == 2:
            age_months = 30
            treatment = "Control"
        else:
            age_months = 30
            treatment = "Corylin"
        meta_rows.append({
            "sample_id": col,
            "dataset": "GSE305103",
            "age_group": "Young" if age_months == 12 else "Old",
            "age_months": age_months,
            "sex": sex,
            "condition": treatment,
            "tissue": "liver",
            "treatment": treatment,
        })
    meta = pd.DataFrame(meta_rows)
    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    return meta, expr


def process_gse221286(model_info: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """RagC mutant vs WT, young liver."""
    print("\n=== GSE221286 ===")
    import gzip
    path = COUNTS_DIR / "GSE221286_Young_RagC_S74N_vs_Young_RagC_wt.normalizedCounts.txt.gz"
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, sep="\t", index_col=0)
    print(f"  Raw shape: {df.shape}")

    sample_cols = [c for c in df.columns]
    genes = df.index.astype(str).tolist()
    mapping = load_or_build_mapping(RESULTS_DIR / "gse221286_symbol_to_entrez.pkl", genes, "symbol")
    entrez_list = []
    expr_list = []
    for g in genes:
        if g in mapping:
            entrez_list.append(str(mapping[g]))
            expr_list.append(df.loc[g].values)
    expr = pd.DataFrame(np.array(expr_list), columns=sample_cols, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    log_expr = np.log2(expr + 1)
    preds = predict_from_logexpr(log_expr, model_info)
    preds_adj = preds * 48.0

    # Need genotype info from series matrix. Based on GEO metadata:
    # YAA00160, YAA00170, YAA00180, YAA00189 = RagC S74N
    # YAA00164, YAA00185, YAA00190, YAA00198 = WT
    mutant_samples = ["YAA00160", "YAA00170", "YAA00180", "YAA00189"]
    meta_rows = []
    for col in sample_cols:
        treatment = "RagC_mutant" if col in mutant_samples else "WT"
        meta_rows.append({
            "sample_id": col,
            "dataset": "GSE221286",
            "age_group": "Young",
            "age_months": 3,
            "sex": "Unknown",
            "condition": treatment,
            "tissue": "liver",
            "treatment": treatment,
        })
    meta = pd.DataFrame(meta_rows)
    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    return meta, expr


def process_gse288795(model_info: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Already processed script exists; re-run here for unified output."""
    print("\n=== GSE288795 ===")
    path = COUNTS_DIR / "GSE288795_count_data.xlsx"
    df_raw = pd.read_excel(path, sheet_name=0, header=None)
    descriptions = df_raw.iloc[0, 1:].tolist()
    sample_ids = [f"S{i+1}" for i in range(len(descriptions))]
    data = df_raw.iloc[2:, 1:].values.astype(float)
    genes = df_raw.iloc[2:, 0].astype(str).tolist()
    print(f"  Samples: {len(sample_ids)}, Genes: {len(genes)}")

    import re
    def parse_meta(desc):
        m = re.match(r'(Male|Female)\s+(\w+)\s+(Rapamycin/Trametinib|Rapamycin|Trametinib|Control)\s+Biol Rep\s+(\d+)', desc)
        if not m:
            return {"sex": "", "tissue": "", "treatment": "", "rep": ""}
        return {"sex": m.group(1), "tissue": m.group(2), "treatment": m.group(3), "rep": int(m.group(4))}

    meta = pd.DataFrame([parse_meta(d) for d in descriptions])
    meta["sample_id"] = sample_ids
    meta["description"] = descriptions

    mapping = load_or_build_mapping(RESULTS_DIR / "gse288795_ensembl_to_entrez.pkl",
                                    [g.split(".")[0] for g in genes], "ensembl")
    entrez_list = []
    expr_list = []
    for g in genes:
        cg = g.split(".")[0]
        if cg in mapping:
            entrez_list.append(str(mapping[cg]))
            expr_list.append(data[genes.index(g)])
    expr = pd.DataFrame(np.array(expr_list), columns=sample_ids, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    log_expr = np.log2(expr + 1)
    preds = predict_from_logexpr(log_expr, model_info)
    preds_adj = preds * 48.0
    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    meta["dataset"] = "GSE288795"
    meta["age_group"] = "Unknown"
    meta["age_months"] = np.nan
    meta["condition"] = meta["treatment"]
    meta["tissue"] = meta["tissue"]
    return meta, expr


def process_gse230402(model_info: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """CR vs AL liver."""
    print("\n=== GSE230402 ===")
    import gzip
    path = COUNTS_DIR / "GSE230402_Counts.normalized_WT.txt.gz"
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(path, sep="\t", compression="gzip")
    genes = df.iloc[:, 0].astype(str).tolist()
    samples = df.columns[1:].tolist()
    counts = df.iloc[:, 1:].values.astype(float)
    print(f"  Samples: {len(samples)}")

    mapping = load_or_build_mapping(RESULTS_DIR / "gse230402_symbol_to_entrez.pkl", genes, "symbol")
    entrez_list = []
    expr_list = []
    for i, g in enumerate(genes):
        if g in mapping:
            entrez_list.append(str(mapping[g]))
            expr_list.append(counts[i])
    expr = pd.DataFrame(np.array(expr_list), columns=samples, index=entrez_list)
    expr = expr.groupby(level=0).mean()
    print(f"  After dedup: {expr.shape}")

    log_expr = np.log2(expr + 1)
    preds = predict_from_logexpr(log_expr, model_info)
    preds_adj = preds * 48.0

    meta = pd.DataFrame({
        "sample_id": samples,
        "sex": ["Female","Male","Male","Female","Female","Female","Male","Female","Male","Female","Male","Male","Male","Female","Male","Male","Female","Female","Female","Female","Male","Male"],
        "condition": ["CR","AL","CR","AL","CR","AL","AL","CR","CR","CR","AL","AL","CR","AL","CR","AL","CR","AL","AL","CR","AL","CR"],
    })
    meta["dataset"] = "GSE230402"
    meta["age_group"] = "Unknown"
    meta["age_months"] = np.nan
    meta["tissue"] = "liver"
    meta["treatment"] = meta["condition"]
    meta["tAge"] = preds
    meta["tAge_adj"] = preds_adj
    return meta, expr


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("Unified GEO Pipeline — Direction 3")
    print("=" * 70)

    model_info = load_model_info()
    print(f"Model features: {len(model_info['feature_names'])}")

    all_meta = []
    processors = [
        process_gse282210,
        process_gse283201,
        process_gse253612,
        process_gse305103,
        process_gse221286,
        process_gse288795,
        process_gse230402,
    ]

    for proc in processors:
        try:
            meta, expr = proc(model_info)
            all_meta.append(meta)
            # Save per-dataset expression for possible future use
            gse = meta["dataset"].iloc[0]
            expr.to_csv(RESULTS_DIR / f"{gse.lower()}_expression.csv")
        except Exception as e:
            print(f"  ERROR processing {proc.__name__}: {e}")
            import traceback
            traceback.print_exc()

    combined = pd.concat(all_meta, ignore_index=True)
    combined.to_csv(RESULTS_DIR / "combined_geo_predictions.csv", index=False)
    print(f"\nCombined predictions saved: {len(combined)} samples")
    print(combined.groupby("dataset").size().to_string())


if __name__ == "__main__":
    main()
