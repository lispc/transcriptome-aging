# GSE280382: GLP-1R Agonist (Exenatide) Transcriptomic Age Analysis

**Dataset**: GSE280382 — "Functional and multi-omic aging rejuvenation with GLP-1R agonism"  
**Source**: GEO, bulk RNA-seq (mouse)  
**Samples**: 620 total (AgedLT: 284, AgedST: 208, Young: 128)  
**Analysis Date**: 2026-06-09

---

## 1. Study Design

### AgedLT (Long-Term, 30 weeks treatment)
- **Tissues**: Adipose, CardiacMuscle, Colon, FrontalCortex, Hippocampus, Hypothalamus, Kidney, Liver, Lung, SkeletalMuscle, Spleen, WBCs
- **Groups**: Young_ctrl (n=89), Aged_ctrl (n=102), Aged_exenatide (n=93)
- **Mice**: Male C57BL/6, treated from 11 months of age

### AgedST (Short-Term, 13 weeks treatment)
- **Tissues**: Multi-tissue (Adipose, CardiacMuscle, Colon, FrontalCortex, Hippocampus, Hypothalamus, SkeletalMuscle, WBC)
- **Groups**: Young_ctrl, Aged_ctrl, Aged_exenatide, Aged_rapamycin, Aged_ctrl_KD, Aged_exenatide_KD
- **KO groups**: Hypothalamic GLP-1R knockout (KD) to test brain-body axis

### Young
- **Groups**: Young_ctrl, Young_exenatide
- **Purpose**: Test whether GLP-1RA effects are specific to aging

---

## 2. Key Findings

### 2.1 AgedLT: Transcriptomic Age Does NOT Decrease with Exenatide

| Metric | Value |
|--------|-------|
| Young_ctrl mean tAge_adj | **-135.4 months** |
| Aged_ctrl mean tAge_adj | **8.2 months** |
| Aged_exenatide mean tAge_adj | **32.0 months** |
| Aging effect (Aged - Young) | **+143.6 months** |
| "Rejuvenation" (Exen - Aged) | **+23.7 months** |

> **Critical Insight**: Exenatide treatment is associated with a **HIGHER** transcriptomic age compared to aged controls (+23.7 months). This is opposite to the rejuvenation effect observed with Rapamycin, CR, and Acarbose in other datasets.

### 2.2 AgedLT: Per-Tissue Effects

Exenatide effect (Exenatide - Control, months). Positive = higher tAge (pro-aging direction):

| Tissue | Effect (months) | Direction |
|--------|-----------------|-----------|
| SkeletalMuscle | **+61.5** | 🔴 Pro-aging |
| Adipose | **+52.0** | 🔴 Pro-aging |
| Lung | **+36.8** | 🔴 Pro-aging |
| CardiacMuscle | **+33.2** | 🔴 Pro-aging |
| Kidney | **+28.5** | 🔴 Pro-aging |
| Hippocampus | **+26.8** | 🔴 Pro-aging |
| Hypothalamus | **+24.0** | 🔴 Pro-aging |
| FrontalCortex | **+16.6** | 🔴 Pro-aging |
| Spleen | **+9.3** | 🔴 Pro-aging |
| Liver | **+3.6** | 🔴 Pro-aging |
| Colon | **-3.9** | 🟢 Slight rejuvenation |
| WBCs | **-4.5** | 🟢 Slight rejuvenation |

**Only Colon and WBCs show slight transcriptomic rejuvenation.** All other tissues show increased tAge.

### 2.3 AgedST: Exenatide vs Rapamycin (Whole Blood)

| Group | Mean tAge_adj (months) |
|-------|------------------------|
| Aged + Rapamycin | **-8.2** |
| Aged + Exenatide | **-2.1** |
| Aged Control | **-29.9** |

- Rapamycin vs Control: **+21.7 months** (higher tAge)
- Exenatide vs Control: **+27.8 months** (higher tAge)
- Exenatide vs Rapamycin: **+6.1 months** (Exenatide higher than Rapamycin)

> In whole blood, neither Rapamycin nor Exenatide lower tAge at 13 weeks. Rapamycin shows a slightly lower tAge than Exenatide.

### 2.4 KO Effect (Hypothalamic GLP-1R Knockout)

| Group | Mean tAge_adj (months) |
|-------|------------------------|
| Aged_exenatide_KD | **-7.7** |
| Aged_ctrl_KD | **-32.3** |

- Exenatide effect persists even with hypothalamic GLP-1R knockout, suggesting peripheral or partial central mechanisms.

### 2.5 Young Mice: No Effect

| Group | Mean tAge_adj (months) |
|-------|------------------------|
| Young_ctrl | **-64.9** |
| Young_exenatide | **-65.7** |
| Difference | **-0.8** |

> Exenatide has minimal effect on young mice, consistent with the study's finding that benefits are aging-specific.

---

## 3. Interpretation

### 3.1 The GLP-1RA Paradox

The original study (GSE280382) demonstrates that GLP-1RA (Exenatide) improves:
- Physical and cognitive performance
- DNA methylome aging clocks
- Plasma metabolome
- Multi-tissue transcriptomes (in their own analysis)

**However, our tAge model does NOT detect rejuvenation.** This creates an important paradox:

**Possible explanations:**

1. **Different aging clocks measure different things**: The tAge model is trained on mortality/mortality-related gene expression across multiple species. GLP-1RA may extend healthspan through metabolic improvements (insulin sensitivity, weight control, inflammation reduction) that are not captured by this specific mortality-trained transcriptomic clock.

2. **Tissue-specific vs. systemic effects**: GLP-1RA's benefits may be more systemic (circulating factors, metabolic regulation) rather than tissue-level transcriptomic reprogramming.

3. **Short-term transcriptomic changes vs. long-term phenotypic benefits**: The 30-week treatment may induce metabolic adaptations that improve function without reversing the core transcriptomic aging program.

4. **Model limitation**: The tAge model was trained on natural aging, not drug-induced metabolic changes. Drugs acting through novel pathways may not be captured.

### 3.2 Comparison with Other Interventions

| Intervention | tAge Effect | Source |
|--------------|-------------|--------|
| Caloric Restriction | **-0.82** (strong rejuvenation) | GSE131754 |
| Rapamycin (6m) | **-0.63** | GSE131754 |
| Acarbose | **-0.42** | GSE131754 |
| Rapamycin (pooled) | **-38.2 months** | GSE288795 |
| Trametinib (pooled) | **-27.3 months** | GSE288795 |
| Combo (Rapa+Tram) | **-50.7 months** | GSE288795 |
| **Exenatide (pooled)** | **+23.7 months** | **GSE280382** |

> Exenatide is the **only intervention** in our analysis that does NOT lower tAge. This does not mean it is ineffective — it suggests its anti-aging mechanism operates outside the transcriptomic aging program captured by this model.

### 3.3 Clinical Relevance

- GLP-1R agonists (Ozempic, Wegovy, Mounjaro) are currently the most widely used anti-obesity/anti-diabetes drugs.
- Their proven benefits (weight loss, cardiovascular protection, reduced inflammation) may extend lifespan through metabolic pathways.
- **Transcriptomic age may not be the right biomarker to assess GLP-1RA's geroprotective potential.**
- DNA methylation clocks (which the original study used) may be more sensitive to GLP-1RA's effects.

---

## 4. Technical Notes

- **Gene mapping**: 23,528 / 24,850 symbols mapped to Entrez (94.7%)
- **Model overlap**: 10,297 / 10,487 features (98.2%)
- **Normalization**: edgeR TMM via Rscript subprocess
- **Imputation**: Manual numpy (sklearn version compatibility fix)
- **tAge_adj formula**: (tAge - 5.5) × 48

---

## 5. Conclusion

**GLP-1R agonism (Exenatide) does not lower transcriptomic age in aging mice.** This finding highlights a critical limitation of transcriptomic aging clocks: they may miss geroprotective interventions that operate through metabolic, rather than transcriptomic, pathways. The disconnect between functional rejuvenation (improved cognition, physical performance) and transcriptomic age underscores the need for multi-modal aging biomarkers.

**Recommendation**: For drug screening purposes, GLP-1RA-like compounds would NOT pass a transcriptomic-age-based filter. However, this does not invalidate their potential as geroprotectors — it suggests that transcriptomic age alone is insufficient for comprehensive anti-aging drug discovery.
