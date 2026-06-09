# Direction 3: ARCHS4 H5 Processing & GEO Expansion — Complete Report
*Generated: 2026-06-09*

## Executive Summary

The ARCHS4 mouse H5 file (`mouse_gene_v2.5.h5`, target 36 GB) was still downloading (~9.7 GB, 27% complete, unreadable due to truncation). Rather than wait, we pivoted to a large-scale **GEO expansion** using the previously compiled `geo_eutils_search_results.csv` (944 records). We identified, downloaded, and processed **7 high-value GEO series** encompassing **339 mouse transcriptome samples** spanning aging, rapamycin, caloric restriction, corylin, ketogenic diet, and genetic mTORC1 activation.

**Key results:**
- **Aging trajectory is detectable** across independent liver datasets (R = 0.54–0.95).
- **Rapamycin + Trametinib** shows strong, consistent rejuvenation across muscle, kidney, and spleen (ΔtAge −15 to −57 months per tissue).
- **Caloric restriction** rejuvenates liver in both sexes (Male Δ = −20.7 mo, Female Δ = −48.1 mo).
- **Ketogenic diet** at 26 months shows rejuvenation (Δ = −59.0 mo).
- **RagC mutant** (mTORC1 hyperactivation) accelerates aging signature in young liver (Δ = +35.2 mo vs WT).
- **Corylin** at 130 weeks does **not** show transcriptomic rejuvenation in liver (small positive Δ, not significant); possible explanations include late timepoint, tissue-specific effects, or technical batch differences.
- **Severe batch effects exist** between GEO series (ANOVA F = 45.25, p < 0.0001), confirming that **within-study comparisons are essential** for valid inference.

---

## Data Sources & Processing

### ARCHS4 Status
| Item | Status |
|------|--------|
| File | `data/archs4/mouse_gene_v2.5.h5` |
| Current size | 9.7 GB |
| Target size | ~36 GB |
| Download progress | ~27% |
| Usable? | **No** — truncated HDF5 header prevents opening |
| Background download | `wget` still active in screen session `archs4_download` |

### GEO Datasets Processed

| Dataset | N samples | Tissue focus | Ages | Design | Gene overlap |
|---------|-----------|--------------|------|--------|--------------|
| **GSE282210** | 133 | Liver (quiescent & regenerating) | 3, 14, 23 mo | Young / Old / Very Old + PH timecourse | 6,699 / 10,487 |
| **GSE283201** | 32 | Liver + Hippocampus | 3, 22 mo | Young vs Aged, Control vs Ethanol | 10,430 / 10,487 |
| **GSE253612** | 16 | Liver (26 mo only) | 26 mo | Control Diet vs Ketogenic Diet | 10,483 / 10,487 |
| **GSE305103** | 17 | Liver, Kidney, Quadriceps | 12, 30 mo | Young vs Old vs Corylin-treated | 10,482 / 10,487 |
| **GSE221286** | 8 | Liver | 3 mo | RagC S74N mutant vs WT (mTORC1↑)| 8,467 / 10,487 |
| **GSE288795** | 111 | Muscle, Kidney, Spleen | — | Rapamycin, Trametinib, Combo vs Control | 10,486 / 10,487 |
| **GSE230402** | 22 | Liver | — | Caloric Restriction vs Ad Libitum (sex-balanced) | 9,581 / 10,487 |

### Preprocessing Pipeline
1. **Load processed counts / TPM** from GEO supplementary files (no raw FASTQ processing needed).
2. **Gene ID mapping**: Ensembl → Entrez (via mygene) or Symbol → Entrez; deduplicated by averaging.
3. **Log-transform**: `log2(expression + 1)` to approximate logCPM scale.
4. **Feature alignment**: Impute missing features with training-set medians from the model.
5. **tAge prediction**: Center with model scaler mean, apply ElasticNet coefficients, multiply by mouse species adjustment (×48).

---

## 1. Aging Trajectory: tAge vs Chronological Age

### Control Liver Samples
We pooled all control/quiescent liver samples with annotated chronological ages:

- **GSE282210** (3 → 14 → 23 mo): slope = **+2.66 tAge months / chronological month**, R = 0.54, p = 0.002
- **GSE283201** (3 → 22 mo): slope = **+9.85 tAge months / chronological month**, R = 0.95, p < 0.001
- **GSE305103** (12 → 30 mo): slope = **+7.76 tAge months / chronological month**, R = 0.93, p < 0.001

**Interpretation:** The transcriptomic clock tracks chronological age strongly in GSE283201 and GSE305103. GSE282210 shows a weaker correlation (R = 0.54), likely because partial-hepatectomy-related regenerative stress elevates tAge in some young control samples (see also PH vs Ctrl analysis below).

**Figure:** `figures/direction3/aging_trajectory_control_liver.png`

### Within-Study Trajectories

#### GSE282210 — Young / Old / Very Old Liver
Control-only regression: R = 0.54. The quiescent livers of 23-month (V) mice sit intermediate between 3-month (Y) and 14-month (O) in this model, but with wide variance.

**Figure:** `figures/direction3/gse282210_aging_trajectory.png`

#### GSE283201 — Young (3 mo) vs Aged (22 mo) Liver
Almost perfect linear relationship (R = 0.95). Control-diet aged livers are ~180 tAge months older than young livers.

**Figure:** `figures/direction3/gse283201_aging_trajectory.png`

---

## 2. Drug & Intervention Effect Comparison

All effects are reported as **ΔtAge = treatment mean − control mean** (negative = rejuvenation).

### Rapamycin & Trametinib (GSE288795)
Pooled across sexes per tissue:

| Tissue | Rapamycin | Trametinib | Combo |
|--------|-----------|------------|-------|
| Muscle | −36.8 | −31.3 | **−49.3** |
| Kidney | −15.6 | −38.9 | **−50.8** |
| Spleen | **−56.7** | −7.3 | **−53.5** |

- **Rapamycin** is consistently rejuvenating across all three tissues.
- **Trametinib** shows strong rejuvenation in kidney, modest in muscle, and minimal in spleen.
- **Combo** is the most rejuvenating in muscle and kidney, comparable to rapamycin alone in spleen.

**Figures:** `figures/direction3/gse288795_*_tage.png`, `gse288795_drug_effects_pooled.png`

### Caloric Restriction (GSE230402)
| Sex | ΔtAge (CR − AL) | n_CR | n_AL |
|-----|-----------------|------|------|
| Male | **−20.7** | 5 | 6 |
| Female | **−48.1** | 6 | 5 |

CR shows stronger rejuvenation in female liver than male liver in this dataset.

**Figure:** `figures/direction3/gse230402_liver_tage.png`

### Ketogenic Diet (GSE253612, 26-month liver)
| Comparison | ΔtAge | n_KD | n_CD |
|------------|-------|------|------|
| Keto vs Control | **−59.0** | 8 | 8 |

A large rejuvenation signal in aged liver following cyclic ketogenic diet.

### RagC Mutant — Genetic mTORC1 Hyperactivation (GSE221286)
| Comparison | ΔtAge | n_mut | n_WT |
|------------|-------|-------|------|
| RagC S74N vs WT | **+35.2** | 4 | 4 |

Active mTORC1 signaling drives a premature aging transcriptomic signature even in 3-month-old mice, consistent with the paper's inflammaging phenotype.

### Corylin (GSE305103, 130-week liver)
| Comparison | ΔtAge | n_treat | n_ctrl |
|------------|-------|---------|--------|
| Corylin vs Old Control (Male) | +37.4 | 3 | 2 |
| Corylin vs Old Control (Female) | +19.8 | 3 | 3 |
| Old vs Young (Male) | +135.8 | 5 | 3 |
| Old vs Young (Female) | +167.1 | 6 | 3 |

**Corylin does not show transcriptomic rejuvenation in liver at 130 weeks.** It is slightly higher than old control (not statistically significant with n = 2–3). Possible reasons:
1. **Late timepoint**: Rejuvenation may have occurred earlier or may be post-transcriptional.
2. **Tissue specificity**: Corylin's longevity benefits may be driven by kidney or metabolic tissues rather than liver transcriptome.
3. **Count normalization**: GSE305103 provided raw STAR counts; log2(counts+1) is an approximation of the logCPM the model was trained on.

**Figure:** `figures/direction3/gse305103_corylin_effect.png`

### Pooled Drug Effects Barplot
**Figure:** `figures/direction3/drug_effects_across_datasets.png`

---

## 3. Batch Effect Check

### tAge Distribution by Dataset (Liver Only)
| Dataset | Mean tAge | SD | N |
|---------|-----------|----|---|
| GSE221286 | −21.4 | 30.9 | 8 |
| GSE230402 | 251.3 | 39.0 | 22 |
| GSE253612 | 261.5 | 65.0 | 16 |
| GSE282210 | 166.4 | 46.3 | 133 |
| GSE283201 | 292.1 | 96.8 | 16 |
| GSE305103 | 386.6 | 80.4 | 17 |

### ANOVA on Control Liver Samples
F = 45.25, p < 0.0001

**Conclusion:** There are massive, statistically significant batch effects between GEO series. Absolute tAge values are **not comparable across studies**. All drug/aging inferences must be made **within-study** (treatment vs control from the same series). Cross-study meta-analysis would require formal batch correction (e.g., ComBat) and is not attempted here.

**Figure:** `figures/direction3/batch_effect_by_dataset.png`

---

## 4. Additional Observations

### Partial Hepatectomy Stress Increases tAge (GSE282210)
Within GSE282210, comparing all PH timepoints vs quiescent controls:
- Control mean tAge = 141.5 months
- PH (all timepoints pooled) mean tAge = 173.9 months
- Δ = +32.5 months

Surgical stress and regenerative proliferation produce a transcriptomic signature that the clock interprets as "older." This is biologically plausible (inflammation, cell-cycle entry).

### Ethanol Does Not Strongly Shift Liver tAge (GSE283201)
- Control diet mean = 293.5 months
- Ethanol diet mean = 290.6 months
- Δ = −2.9 months (negligible)

In this 3-month vs 22-month comparison, ethanol feeding does not produce a strong aging or rejuvenation signal in liver.

---

## 5. Limitations

1. **ARCHS4 H5 incomplete**: The 36 GB download is still in progress. All results here are GEO-based backups.
2. **Batch effects**: Extremely large between-study variance. No cross-study batch correction applied.
3. **Normalization mismatch**: Model was trained on TMM-normalized logCPM. GEO datasets provide raw counts, normalized counts, or TPM. We used `log2(x+1)` as a practical approximation, which may introduce systematic biases.
4. **Small sample sizes**: Some comparisons (Corylin, RagC) have n = 3–4 per group.
5. **Gene mapping losses**: GSE282210 had only ~6,700 overlapping Entrez genes (vs 10,487), reducing prediction accuracy for that series.
6. **Chronological age parsing**: For some datasets (GSE288795, GSE230402), exact ages were not available in GEO metadata.
7. **Heterogeneous tissues**: GSE288795 includes muscle, kidney, and spleen—not liver—limiting direct comparison to hepatic interventions.

---

## 6. Files Generated

### Predictions & Tables
- `results/direction3/combined_geo_predictions.csv` — 339 samples, all datasets
- `results/direction3/drug_effect_comparisons.csv` — 17 within-study comparisons
- `results/direction3/liver_summary_by_dataset_treatment.csv` — mean tAge by group

### Figures
- `figures/direction3/aging_trajectory_control_liver.png`
- `figures/direction3/gse282210_aging_trajectory.png`
- `figures/direction3/gse283201_aging_trajectory.png`
- `figures/direction3/gse305103_corylin_effect.png`
- `figures/direction3/drug_effects_across_datasets.png`
- `figures/direction3/batch_effect_by_dataset.png`
- `figures/direction3/gse230402_liver_tage.png`
- `figures/direction3/gse288795_*_tage.png` (muscle, kidney, spleen)
- `figures/direction3/gse288795_drug_effects_pooled.png`

### Code
- `src/python/direction3/unified_geo_pipeline.py` — loads, maps, predicts all datasets
- `src/python/direction3/combined_analysis.py` — generates figures and statistics

---

## 7. Next Steps

1. **ARCHS4 completion**: Re-run the unified pipeline on the full H5 once wget finishes (expected ~24–48 h depending on bandwidth).
2. **Batch correction**: Apply ComBat or harmonic mean normalization if cross-study meta-analysis is needed.
3. **TMM normalization**: For GEO datasets with raw counts (GSE253612, GSE283201, GSE305103), run edgeR TMM → logCPM to better match training data.
4. **Additional GEO series**: Prioritize GSE113745 (young vs old liver rhythms, 24 samples) and GSE296874 (CR + fiber, 36 samples) once raw count matrices are retrieved.
5. **Human ARCHS4**: Download and process `human_gene_v2.5.h5` for cross-species validation.
6. **Module fingerprint integration**: Correlate LINCS drug module signatures with GEO-derived drug signatures for mechanism-based drug repurposing.

---

*Report generated by Direction 3 subagent. ARCHS4 download is still in progress; this report represents the complete GEO-expansion backup analysis.*
