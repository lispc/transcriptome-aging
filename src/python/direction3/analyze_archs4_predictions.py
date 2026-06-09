#!/usr/bin/env python3
"""
Analyze ARCHS4 tAge predictions:
- Aging trajectory
- Drug effects (within-study comparisons)
- Batch effect assessment
- Dose-response where available
"""
import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
FIGURES_DIR = PROJECT_ROOT / "figures" / "direction3"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 150
sns.set_style("whitegrid")

def load_all_predictions():
    f = RESULTS_DIR / "archs4_all_predictions.csv"
    if not f.exists():
        print(f"Combined predictions not found: {f}")
        sys.exit(1)
    df = pd.read_csv(f, low_memory=False)
    print(f"Loaded {len(df)} predictions")
    return df

def parse_age_from_text(text):
    """Heuristic extraction of age in months from sample text."""
    if pd.isna(text):
        return np.nan
    text = str(text).lower()
    # Look for patterns like "12 months", "12m", "12 month"
    import re
    # months
    m = re.search(r'(\d+(?:\.\d+)?)\s*(?:months?|mo|m)(?:\s*old)?', text)
    if m:
        return float(m.group(1))
    # weeks
    m = re.search(r'(\d+(?:\.\d+)?)\s*(?:weeks?|wk|w)(?:\s*old)?', text)
    if m:
        return float(m.group(1)) / 4.33
    # years
    m = re.search(r'(\d+(?:\.\d+)?)\s*(?:years?|yr|y)(?:\s*old)?', text)
    if m:
        return float(m.group(1)) * 12
    # days
    m = re.search(r'(\d+(?:\.\d+)?)\s*(?:days?|d)(?:\s*old)?', text)
    if m:
        return float(m.group(1)) / 30.44
    return np.nan

def add_parsed_age(df):
    text = df.get("Sample_title", "") + " " + df.get("Sample_characteristics_ch1", "") + " " + df.get("Sample_source_name_ch1", "")
    text = text.fillna("")
    df["parsed_age_months"] = text.apply(parse_age_from_text)
    return df

def plot_aging_trajectory(df):
    ctrl = df[
        ~df["Sample_title"].str.contains("rapamycin|metformin|dexamethasone|doxorubicin|cr|caloric|restriction|treated|drug", case=False, na=False) &
        df["parsed_age_months"].notna()
    ].copy()
    if len(ctrl) < 10:
        print("Not enough control samples with parsed age for trajectory plot")
        return

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=ctrl, x="parsed_age_months", y="tAge_adj", alpha=0.3, s=20)
    # Add lowess
    try:
        sns.regplot(data=ctrl, x="parsed_age_months", y="tAge_adj", scatter=False, lowess=True, color="red")
    except Exception as e:
        print(f"  lowess failed: {e}")
    plt.xlabel("Chronological Age (months)")
    plt.ylabel("Predicted tAge (months)")
    plt.title("Aging Trajectory — ARCHS4 Control Samples")
    plt.tight_layout()
    out = FIGURES_DIR / "archs4_aging_trajectory.png"
    plt.savefig(out)
    print(f"Saved {out}")
    plt.close()

def plot_drug_effects(df):
    drugs = {
        "rapamycin": "Rapamycin",
        "metformin": "Metformin",
        "dexamethasone": "Dexamethasone",
        "doxorubicin": "Doxorubicin",
    }
    results = []
    for key, label in drugs.items():
        drug_mask = df["Sample_title"].str.contains(key, case=False, na=False) | \
                    df["Sample_characteristics_ch1"].str.contains(key, case=False, na=False)
        drug_samples = df[drug_mask]
        if len(drug_samples) == 0:
            continue
        # Try to find controls in same series
        series_ids = drug_samples["Series_geo_accession"].dropna().unique()
        ctrl_mask = df["Series_geo_accession"].isin(series_ids) & \
                    ~df["Sample_title"].str.contains(key, case=False, na=False) & \
                    ~df["Sample_characteristics_ch1"].str.contains(key, case=False, na=False)
        ctrl_samples = df[ctrl_mask]
        if len(ctrl_samples) == 0:
            print(f"  {label}: no controls in same series")
            continue
        diff = drug_samples["tAge_adj"].mean() - ctrl_samples["tAge_adj"].mean()
        results.append({
            "drug": label,
            "n_drug": len(drug_samples),
            "n_control": len(ctrl_samples),
            "drug_mean": drug_samples["tAge_adj"].mean(),
            "ctrl_mean": ctrl_samples["tAge_adj"].mean(),
            "difference": diff,
            "drug_sem": drug_samples["tAge_adj"].sem(),
            "ctrl_sem": ctrl_samples["tAge_adj"].sem(),
        })
    if not results:
        print("No drug-control pairs found")
        return
    res_df = pd.DataFrame(results)
    print("\nDrug effects (negative = rejuvenation):")
    print(res_df.to_string(index=False))
    res_df.to_csv(RESULTS_DIR / "archs4_drug_effects.csv", index=False)

    plt.figure(figsize=(8, 5))
    x = np.arange(len(res_df))
    plt.bar(x, res_df["difference"], color=["green" if d < 0 else "red" for d in res_df["difference"]])
    plt.xticks(x, res_df["drug"], rotation=45, ha="right")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.ylabel("Drug - Control tAge (months)")
    plt.title("Drug Effects on tAge (within-study comparisons)")
    plt.tight_layout()
    out = FIGURES_DIR / "archs4_drug_effects.png"
    plt.savefig(out)
    print(f"Saved {out}")
    plt.close()

def plot_batch_effects(df):
    # Top series by sample count
    top_series = df["Series_geo_accession"].value_counts().head(15).index.tolist()
    sub = df[df["Series_geo_accession"].isin(top_series)].copy()
    if len(sub) < 10:
        print("Not enough data for batch effect plot")
        return

    plt.figure(figsize=(12, 6))
    sns.boxplot(data=sub, x="Series_geo_accession", y="tAge_adj")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Predicted tAge (months)")
    plt.title("tAge Distribution by GEO Series (batch effect check)")
    plt.tight_layout()
    out = FIGURES_DIR / "archs4_batch_effects.png"
    plt.savefig(out)
    print(f"Saved {out}")
    plt.close()

def plot_cr_dose_response(df):
    cr_mask = df["Sample_title"].str.contains("caloric restriction|calorie restriction|cr diet", case=False, na=False) | \
              df["Sample_characteristics_ch1"].str.contains("caloric restriction|calorie restriction|cr diet", case=False, na=False)
    cr = df[cr_mask].copy()
    if len(cr) == 0:
        print("No CR samples found for dose-response")
        return
    # Heuristic: search for percentages
    import re
    def extract_pct(text):
        if pd.isna(text):
            return np.nan
        text = str(text).lower()
        m = re.search(r'(\d+)%?\s*(?:caloric restriction|cr|restriction)', text)
        if m:
            return float(m.group(1))
        m = re.search(r'cr\s*(\d+)%?', text)
        if m:
            return float(m.group(1))
        return np.nan
    text = cr.get("Sample_title", "") + " " + cr.get("Sample_characteristics_ch1", "")
    text = text.fillna("")
    cr["cr_pct"] = text.apply(extract_pct)
    if cr["cr_pct"].notna().sum() < 3:
        print("Not enough CR dose info for dose-response")
        return
    # Match controls per series
    series_ids = cr["Series_geo_accession"].dropna().unique()
    ctrl_mask = df["Series_geo_accession"].isin(series_ids) & ~cr_mask
    ctrl = df[ctrl_mask]
    ctrl_mean_by_series = ctrl.groupby("Series_geo_accession")["tAge_adj"].mean().to_dict()
    cr["ctrl_mean"] = cr["Series_geo_accession"].map(ctrl_mean_by_series)
    cr["diff"] = cr["tAge_adj"] - cr["ctrl_mean"]
    cr_valid = cr[cr["cr_pct"].notna() & cr["ctrl_mean"].notna()]
    if len(cr_valid) < 3:
        print("Not enough matched CR dose data")
        return
    plt.figure(figsize=(7, 5))
    sns.scatterplot(data=cr_valid, x="cr_pct", y="diff", hue="Series_geo_accession", legend=False, s=60)
    plt.xlabel("CR Level (%)")
    plt.ylabel("CR - Control tAge (months)")
    plt.title("CR Dose-Response")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.tight_layout()
    out = FIGURES_DIR / "archs4_cr_dose_response.png"
    plt.savefig(out)
    print(f"Saved {out}")
    plt.close()

def main():
    df = load_all_predictions()
    df = add_parsed_age(df)
    plot_aging_trajectory(df)
    plot_drug_effects(df)
    plot_batch_effects(df)
    plot_cr_dose_response(df)
    print("\nAnalysis complete.")

if __name__ == "__main__":
    main()
