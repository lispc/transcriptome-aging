#!/usr/bin/env python3
"""
Combined analysis of all GEO datasets for Direction 3.
Generates figures and report.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "direction3"
FIGURES_DIR = PROJECT_ROOT / "figures" / "direction3"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 150
sns.set_style("whitegrid")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv(RESULTS_DIR / "combined_geo_predictions.csv")
print(f"Loaded {len(df)} predictions from {df['dataset'].nunique()} datasets")

# ---------------------------------------------------------------------------
# 1. Aging trajectory: tAge vs chronological age for control liver samples
# ---------------------------------------------------------------------------
print("\n--- Aging Trajectory ---")
control_mask = (
    (df["treatment"].isin(["Control", "WT", "AL"])) &
    (df["tissue"] == "liver") &
    (df["age_months"].notna())
)
ctrl_df = df[control_mask].copy()
print(f"Control samples with age: {len(ctrl_df)}")
print(ctrl_df.groupby(["dataset", "age_group"]).size().to_string())

fig, ax = plt.subplots(figsize=(8, 6))
colors = sns.color_palette("tab10", n_colors=ctrl_df["dataset"].nunique())
for (dataset, color) in zip(ctrl_df["dataset"].unique(), colors):
    sub = ctrl_df[ctrl_df["dataset"] == dataset]
    ax.scatter(sub["age_months"], sub["tAge_adj"], label=dataset, s=80, alpha=0.7, color=color)
    # Fit line per dataset if >=3 points and >1 unique age
    if len(sub) >= 3 and sub["age_months"].nunique() > 1:
        slope, intercept, r, p, se = stats.linregress(sub["age_months"], sub["tAge_adj"])
        x_line = np.array([sub["age_months"].min(), sub["age_months"].max()])
        ax.plot(x_line, slope * x_line + intercept, color=color, linestyle="--", alpha=0.5)
        print(f"  {dataset}: slope={slope:.2f} months/month, R={r:.2f}, p={p:.3f}")

ax.set_xlabel("Chronological Age (months)")
ax.set_ylabel("Predicted tAge (months)")
ax.set_title("Aging Trajectory — Control Liver Samples")
ax.legend(title="Dataset", loc="upper left", fontsize=8)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "aging_trajectory_control_liver.png")
plt.close(fig)
print(f"  Saved aging_trajectory_control_liver.png")

# ---------------------------------------------------------------------------
# 2. Drug effect comparison across studies
# ---------------------------------------------------------------------------
print("\n--- Drug Effect Comparison ---")

# Define comparisons per dataset
comparisons = []

# GSE288795 (Rapamycin, Trametinib, Combo)
for tissue in ["Muscle", "Kidney", "Spleen"]:
    for drug in ["Rapamycin", "Trametinib", "Rapamycin/Trametinib"]:
        sub = df[(df["dataset"] == "GSE288795") & (df["tissue"] == tissue)]
        ctrl = sub[sub["treatment"] == "Control"]["tAge_adj"]
        drug_vals = sub[sub["treatment"] == drug]["tAge_adj"]
        if len(ctrl) > 0 and len(drug_vals) > 0:
            comparisons.append({
                "dataset": "GSE288795",
                "comparison": f"{drug} vs Control ({tissue})",
                "drug": drug,
                "tissue": tissue,
                "delta": drug_vals.mean() - ctrl.mean(),
                "ctrl_mean": ctrl.mean(),
                "drug_mean": drug_vals.mean(),
                "n_ctrl": len(ctrl),
                "n_drug": len(drug_vals),
            })

# GSE230402 (CR vs AL)
for sex in ["Male", "Female"]:
    sub = df[(df["dataset"] == "GSE230402") & (df["sex"] == sex)]
    ctrl = sub[sub["treatment"] == "AL"]["tAge_adj"]
    drug_vals = sub[sub["treatment"] == "CR"]["tAge_adj"]
    if len(ctrl) > 0 and len(drug_vals) > 0:
        comparisons.append({
            "dataset": "GSE230402",
            "comparison": f"CR vs AL ({sex})",
            "drug": "CR",
            "tissue": "liver",
            "delta": drug_vals.mean() - ctrl.mean(),
            "ctrl_mean": ctrl.mean(),
            "drug_mean": drug_vals.mean(),
            "n_ctrl": len(ctrl),
            "n_drug": len(drug_vals),
        })

# GSE305103 (Corylin vs Old Control only; Young Control is separate)
for sex in ["Male", "Female"]:
    sub = df[(df["dataset"] == "GSE305103") & (df["sex"] == sex) & (df["age_group"] == "Old")]
    ctrl = sub[sub["treatment"] == "Control"]["tAge_adj"]
    drug_vals = sub[sub["treatment"] == "Corylin"]["tAge_adj"]
    if len(ctrl) > 0 and len(drug_vals) > 0:
        comparisons.append({
            "dataset": "GSE305103",
            "comparison": f"Corylin vs Old Control ({sex})",
            "drug": "Corylin",
            "tissue": "liver",
            "delta": drug_vals.mean() - ctrl.mean(),
            "ctrl_mean": ctrl.mean(),
            "drug_mean": drug_vals.mean(),
            "n_ctrl": len(ctrl),
            "n_drug": len(drug_vals),
        })
# Also add Young vs Old within GSE305103
for sex in ["Male", "Female"]:
    sub = df[(df["dataset"] == "GSE305103") & (df["sex"] == sex)]
    young = sub[sub["age_group"] == "Young"]["tAge_adj"]
    old = sub[sub["age_group"] == "Old"]["tAge_adj"]
    if len(young) > 0 and len(old) > 0:
        comparisons.append({
            "dataset": "GSE305103",
            "comparison": f"Old vs Young ({sex})",
            "drug": "Aging",
            "tissue": "liver",
            "delta": old.mean() - young.mean(),
            "ctrl_mean": young.mean(),
            "drug_mean": old.mean(),
            "n_ctrl": len(young),
            "n_drug": len(old),
        })

# GSE253612 (Keto vs Control at 26mo)
sub = df[df["dataset"] == "GSE253612"]
ctrl = sub[sub["treatment"] == "Control"]["tAge_adj"]
drug_vals = sub[sub["treatment"] == "Keto"]["tAge_adj"]
if len(ctrl) > 0 and len(drug_vals) > 0:
    comparisons.append({
        "dataset": "GSE253612",
        "comparison": "Keto vs Control (26mo liver)",
        "drug": "Keto",
        "tissue": "liver",
        "delta": drug_vals.mean() - ctrl.mean(),
        "ctrl_mean": ctrl.mean(),
        "drug_mean": drug_vals.mean(),
        "n_ctrl": len(ctrl),
        "n_drug": len(drug_vals),
    })

# GSE221286 (RagC mutant vs WT)
sub = df[df["dataset"] == "GSE221286"]
ctrl = sub[sub["treatment"] == "WT"]["tAge_adj"]
drug_vals = sub[sub["treatment"] == "RagC_mutant"]["tAge_adj"]
if len(ctrl) > 0 and len(drug_vals) > 0:
    comparisons.append({
        "dataset": "GSE221286",
        "comparison": "RagC mutant vs WT (3mo liver)",
        "drug": "RagC_mutant",
        "tissue": "liver",
        "delta": drug_vals.mean() - ctrl.mean(),
        "ctrl_mean": ctrl.mean(),
        "drug_mean": drug_vals.mean(),
        "n_ctrl": len(ctrl),
        "n_drug": len(drug_vals),
    })

comp_df = pd.DataFrame(comparisons)
comp_df.to_csv(RESULTS_DIR / "drug_effect_comparisons.csv", index=False)
print(f"Comparisons: {len(comp_df)}")
print(comp_df[["comparison", "delta", "n_drug", "n_ctrl"]].to_string(index=False))

# Barplot of drug effects
fig, ax = plt.subplots(figsize=(10, 6))
y_pos = np.arange(len(comp_df))
colors = ["green" if d < 0 else "red" for d in comp_df["delta"]]
ax.barh(y_pos, comp_df["delta"], color=colors, alpha=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(comp_df["comparison"], fontsize=8)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("Δ tAge (drug - control, months)")
ax.set_title("Drug/Intervention Effects Across GEO Datasets")
ax.invert_yaxis()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "drug_effects_across_datasets.png")
plt.close(fig)
print(f"  Saved drug_effects_across_datasets.png")

# ---------------------------------------------------------------------------
# 3. Batch effect check: tAge distribution by dataset
# ---------------------------------------------------------------------------
print("\n--- Batch Effects ---")
liver_df = df[df["tissue"] == "liver"].copy()
print("Mean tAge by dataset (liver only):")
print(liver_df.groupby("dataset")["tAge_adj"].agg(["mean", "std", "count"]).to_string())

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=liver_df, x="dataset", y="tAge_adj", ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_ylabel("Predicted tAge (months)")
ax.set_title("tAge Distribution by Dataset (Liver Samples)")
plt.tight_layout()
fig.savefig(FIGURES_DIR / "batch_effect_by_dataset.png")
plt.close(fig)
print(f"  Saved batch_effect_by_dataset.png")

# ANOVA across datasets for control samples
ctrl_liver = liver_df[liver_df["treatment"].isin(["Control", "WT", "AL"])]
if ctrl_liver["dataset"].nunique() > 1:
    groups = [g["tAge_adj"].values for _, g in ctrl_liver.groupby("dataset")]
    f_stat, p_val = stats.f_oneway(*groups)
    print(f"  ANOVA across control liver samples: F={f_stat:.2f}, p={p_val:.4f}")

# ---------------------------------------------------------------------------
# 4. Within-study aging trajectories
# ---------------------------------------------------------------------------
print("\n--- Within-Study Aging Trajectories ---")

# GSE282210 detailed
gse282210 = df[(df["dataset"] == "GSE282210") & (df["treatment"] == "Control")]
fig, ax = plt.subplots(figsize=(7, 5))
for sex in ["Female", "Male"]:
    sub = gse282210[gse282210["sex"] == sex]
    ax.scatter(sub["age_months"], sub["tAge_adj"], label=sex, s=80, alpha=0.7)
if len(gse282210) >= 3:
    slope, intercept, r, p, se = stats.linregress(gse282210["age_months"], gse282210["tAge_adj"])
    x_line = np.array([gse282210["age_months"].min(), gse282210["age_months"].max()])
    ax.plot(x_line, slope * x_line + intercept, "k--", alpha=0.5,
            label=f"Fit: R={r:.2f}, p={p:.3f}")
ax.set_xlabel("Chronological Age (months)")
ax.set_ylabel("Predicted tAge (months)")
ax.set_title("GSE282210 — Control Liver tAge vs Age")
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "gse282210_aging_trajectory.png")
plt.close(fig)
print(f"  Saved gse282210_aging_trajectory.png")

# GSE283201 detailed
gse283201 = df[(df["dataset"] == "GSE283201") & (df["tissue"] == "liver") & (df["treatment"] == "Control")]
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(gse283201["age_months"], gse283201["tAge_adj"], s=100, alpha=0.7, color="steelblue")
if len(gse283201) >= 3:
    slope, intercept, r, p, se = stats.linregress(gse283201["age_months"], gse283201["tAge_adj"])
    x_line = np.array([gse283201["age_months"].min(), gse283201["age_months"].max()])
    ax.plot(x_line, slope * x_line + intercept, "k--", alpha=0.5,
            label=f"Fit: R={r:.2f}, p={p:.3f}")
ax.set_xlabel("Chronological Age (months)")
ax.set_ylabel("Predicted tAge (months)")
ax.set_title("GSE283201 — Control Liver tAge vs Age")
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "gse283201_aging_trajectory.png")
plt.close(fig)
print(f"  Saved gse283201_aging_trajectory.png")

# GSE305103 detailed (Young vs Old vs Corylin)
gse305103 = df[(df["dataset"] == "GSE305103")]
fig, ax = plt.subplots(figsize=(7, 5))
sns.boxplot(data=gse305103, x="age_group", y="tAge_adj", hue="treatment", ax=ax,
            order=["Young", "Old"])
ax.set_ylabel("Predicted tAge (months)")
ax.set_title("GSE305103 — Corylin Effect on Liver tAge")
plt.tight_layout()
fig.savefig(FIGURES_DIR / "gse305103_corylin_effect.png")
plt.close(fig)
print(f"  Saved gse305103_corylin_effect.png")

# ---------------------------------------------------------------------------
# 5. Summary table for report
# ---------------------------------------------------------------------------
print("\n--- Summary Statistics ---")
summary = liver_df.groupby(["dataset", "treatment"])["tAge_adj"].agg(["mean", "std", "count"])
print(summary.to_string())
summary.to_csv(RESULTS_DIR / "liver_summary_by_dataset_treatment.csv")

print("\nAnalysis complete!")
