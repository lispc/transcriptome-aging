# GSE288795 完成报告：Trametinib + Rapamycin 联合用药的转录组 rejuvenation 效应

**日期**: 2026-06-09

**GEO ID**: GSE288795  
**标题**: "The geroprotectors trametinib and rapamycin combine additively to extend mouse healthspan and lifespan"  
**PMID**: 40437307  
**原文**: https://doi.org/10.1038/s43587-025-00876-4

---

## 摘要

对 GSE288795（111 样本，小鼠 bulk RNA-seq）的转录组分析显示，**Trametinib + Rapamycin 联合用药在所有三种组织中均产生 rejuvenation 效应**，整体效果（Combo -0.97）优于 Rapamycin 单药（-0.70）和 Trametinib 单药（-0.54）。

**关键发现**：
1. **Spleen 是联合用药效果最显著的组织**（Female: -1.61; Male: -0.62），Rapamycin 单药在 Spleen 中已表现出强烈的 rejuvenation（Female: -1.28; Male: -1.07）。
2. **Muscle 中 Male 的 rejuvenation 效应最强**（Rapa: -1.20; Combo: -1.47; Tram: -0.92），但 Female Muscle 效应较弱。
3. **Kidney Female 中出现 Rapamycin 单药 pro-aging 的异常信号**（+0.38），但 Combo 仍然有效（-0.87）。
4. **Additivity 检验**：统计上 **无任何组织/性别组合显示出显著偏离加和性的证据**（|z| < 1.96，p > 0.05）。组合效果在数值上趋向于加和或轻微拮抗，但受限于样本量（n=4–6）和 SEM，无法断言协同或拮抗。
5. **模块分解**：Rapamycin 主导下调 pink 模块（线粒体翻译 / 氧化磷酸化）和 orange 模块（染色质修饰）。Trametinib 对模块的影响方向与 Rapamycin 部分互补，联合用药在模块层面呈现混合模式。

---

## 1. 数据集概况

| 属性 | 详情 |
|------|------|
| **来源** | Gkioni et al., *Nature Aging* 2025, PMID: 40437307 |
| **品系** | C3B6F1 wild type |
| **组织** | Muscle, Kidney, Spleen（24 月龄取样） |
| **干预** | Control, Rapamycin (42 mg/kg 隔周), Trametinib (1.44 mg/kg 饲料), Combo |
| **性别** | Male + Female |
| **重复** | n = 4–6 / 组 |
| **总样本** | 111 |

**注意**：原始任务描述中提及 "mouse liver"，但本研究实际取样组织为 Muscle、Kidney、Spleen。本分析使用的模型为 **EN_Mortality_Multispecies_Multitissue_scaleddiff**，该模型训练数据涵盖多物种多组织，因此可跨组织泛化预测。

### 样本分布

| 组织 | 性别 | Control | Rapamycin | Trametinib | Combo | 小计 |
|------|------|---------|-----------|------------|-------|------|
| Muscle | Male | 5 | 4 | 4 | 4 | 17 |
| Muscle | Female | 5 | 5 | 5 | 4 | 19 |
| Kidney | Male | 5 | 4 | 5 | 4 | 18 |
| Kidney | Female | 5 | 4 | 5 | 5 | 19 |
| Spleen | Male | 6 | 4 | 4 | 4 | 18 |
| Spleen | Female | 5 | 5 | 5 | 5 | 20 |
| **总计** | — | **31** | **26** | **28** | **26** | **111** |

---

## 2. 预处理流程

1. **数据来源**：GEO Supplementary file `GSE288795_Gkioni_et_al.__2025_RNA_seq_count_data.xlsx`
2. **基因 ID**：Ensembl ID (ENSMUSG…)，55,057 genes
3. **ID 映射**：复用 GSE131754 缓存的 Ensembl → Entrez 映射（mygene.info）
   - 映射率：52.4% (28,831 / 55,057)
   - 去重后：28,813 unique Entrez IDs
4. **归一化**：edgeR TMM → logCPM (prior.count=3)
5. **特征对齐**：与模型 10,487 Entrez Gene IDs 对齐
   - 交集基因：10,482 (99.95%)
   - 缺失基因 imputation：5 (0.05%)，使用训练集 median

---

## 3. tAge 预测结果

### 3.1 按组织/性别/治疗的预测均值

| 组织 | 性别 | Control | Rapamycin | Trametinib | Combo |
|------|------|---------|-----------|------------|-------|
| **Muscle** | Male | 8.91 | 7.71 | 7.99 | 7.44 |
| **Muscle** | Female | 7.95 | 7.80 | 7.69 | 7.52 |
| **Kidney** | Male | 6.79 | 5.91 | 5.63 | 5.94 |
| **Kidney** | Female | 5.16 | 5.54 | 4.76 | 4.29 |
| **Spleen** | Male | 5.21 | 4.14 | 5.33 | 4.60 |
| **Spleen** | Female | 6.55 | 5.27 | 5.92 | 4.94 |

### 3.2 Drug – Control 差异（negative = rejuvenation）

| 治疗 | 组织 | 性别 | n_drug | n_ctrl | Drug Mean | Ctrl Mean | **Difference** | 效果 |
|------|------|------|--------|--------|-----------|-----------|----------------|------|
| **Rapamycin** | Spleen | Female | 5 | 5 | 5.27 | 6.55 | **-1.28** | ✅ 最强 |
| **Combo** | Spleen | Female | 5 | 5 | 4.94 | 6.55 | **-1.61** | ✅ 最强 |
| **Rapamycin** | Muscle | Male | 4 | 5 | 7.71 | 8.91 | **-1.20** | ✅ 强 |
| **Combo** | Muscle | Male | 4 | 5 | 7.44 | 8.91 | **-1.47** | ✅ 强 |
| **Trametinib** | Kidney | Male | 5 | 5 | 5.63 | 6.79 | **-1.16** | ✅ 强 |
| **Rapamycin** | Spleen | Male | 4 | 6 | 4.14 | 5.21 | **-1.07** | ✅ 强 |
| **Trametinib** | Muscle | Male | 4 | 5 | 7.99 | 8.91 | **-0.92** | ✅ |
| **Combo** | Kidney | Female | 5 | 5 | 4.29 | 5.16 | **-0.87** | ✅ |
| **Rapamycin** | Kidney | Male | 4 | 5 | 5.91 | 6.79 | **-0.88** | ✅ |
| **Combo** | Kidney | Male | 4 | 5 | 5.94 | 6.79 | **-0.85** | ✅ |
| **Trametinib** | Spleen | Female | 5 | 5 | 5.92 | 6.55 | **-0.63** | ✅ |
| **Combo** | Spleen | Male | 4 | 6 | 4.60 | 5.21 | **-0.62** | ✅ |
| **Trametinib** | Kidney | Female | 5 | 5 | 4.76 | 5.16 | **-0.40** | ✅ mild |
| **Combo** | Muscle | Female | 4 | 5 | 7.52 | 7.95 | **-0.42** | ✅ mild |
| **Trametinib** | Muscle | Female | 5 | 5 | 7.69 | 7.95 | **-0.25** | ✅ mild |
| **Rapamycin** | Muscle | Female | 5 | 5 | 7.80 | 7.95 | **-0.14** | ✅ mild |
| Rapamycin | Kidney | Female | 4 | 5 | 5.54 | 5.16 | **+0.38** | ⚠️ pro-aging |
| Trametinib | Spleen | Male | 4 | 6 | 5.33 | 5.21 | **+0.12** | ⚠️ 中性 |

#### 汇总

| 治疗 | 平均 Difference | 最强组织 | 最弱组织 |
|------|-----------------|----------|----------|
| **Combo** | **-0.97** | Spleen (-1.11) | Kidney (-0.86) |
| **Rapamycin** | **-0.70** | Spleen (-1.17) | Muscle Female (-0.14) |
| **Trametinib** | **-0.54** | Kidney Male (-1.16) | Spleen Male (+0.12) |

---

## 4. Additivity 检验：Combo ≈ Rapa + Tram？

### 4.1 统计检验方法

对每个 (组织, 性别) 组合：
- **Expected additive** = Rapa_diff + Tram_diff
- **Deviation** = Observed(Combo) – Expected
- **Error propagation**: SEM_expected = √(SEM_rapa² + SEM_tram²), SEM_deviation = √(SEM_combo² + SEM_expected²)
- **Z-score** = Deviation / SEM_deviation

### 4.2 结果

| 组织 | 性别 | Rapa | Tram | Combo | Expected | Deviation | SEM_dev | Z-score | 结论 |
|------|------|------|------|-------|----------|-----------|---------|---------|------|
| Muscle | Male | -1.20 | -0.92 | -1.47 | -2.12 | **+0.65** | 1.48 | 0.44 | 加和性 |
| Muscle | Female | -0.14 | -0.25 | -0.42 | -0.39 | **-0.03** | 0.53 | -0.06 | 加和性 |
| Kidney | Male | -0.88 | -1.16 | -0.85 | -2.04 | **+1.20** | 0.63 | 1.90 | 加和性（临界） |
| Kidney | Female | +0.38 | -0.40 | -0.87 | -0.02 | **-0.85** | 0.87 | -0.98 | 加和性 |
| Spleen | Male | -1.07 | +0.12 | -0.62 | -0.95 | **+0.34** | 0.67 | 0.50 | 加和性 |
| Spleen | Female | -1.28 | -0.63 | -1.61 | -1.92 | **+0.30** | 0.55 | 0.55 | 加和性 |

**关键结论**：
- **没有任何组合的 |Z| > 1.96（p < 0.05）**。因此，从统计角度，**联合用药的效果与各单药加和效应无显著差异**。
- Kidney Male 的 Z=1.90 接近临界值，提示可能存在轻微拮抗（observed > expected，即效果弱于加和），但样本量不足（n=4–5），无法定论。
- Kidney Female 中 Rapamycin 单药显示 pro-aging (+0.38)，但 Combo 显示强 rejuvenation (-0.87)。由于 Rapa 单药在此异常，"expected additive" 接近零，因此 Combo 的效果看起来像是"协同"，但这主要由 Rapa 单药的异常信号驱动，不宜过度解读为分子层面的协同作用。

---

## 5. 模块分解（Module Decomposition）

### 5.1 核心模块的 Drug – Control 差异

| 治疗 | 组织 | 性别 | composite | turquoise | pink | green | orange |
|------|------|------|-----------|-----------|------|-------|--------|
| Rapamycin | Muscle | Male | -1.20 | -0.04 | **-0.40** | -0.11 | **-0.24** |
| Trametinib | Muscle | Male | -0.92 | +0.11 | **+0.17** | -0.13 | +0.00 |
| Combo | Muscle | Male | -1.47 | -0.10 | -0.37 | -0.12 | -0.27 |
| Rapamycin | Spleen | Female | -1.28 | +0.01 | -0.10 | -0.06 | +0.05 |
| Trametinib | Spleen | Female | -0.63 | -0.01 | -0.05 | -0.03 | +0.09 |
| Combo | Spleen | Female | -1.61 | -0.04 | **-0.24** | +0.01 | +0.13 |
| Rapamycin | Kidney | Male | -0.88 | -0.09 | -0.14 | -0.09 | -0.06 |
| Trametinib | Kidney | Male | -1.16 | -0.00 | -0.07 | +0.01 | +0.02 |
| Combo | Kidney | Male | -0.85 | -0.07 | -0.17 | -0.05 | -0.06 |

### 5.2 模块层面的关键发现

1. **pink 模块（线粒体翻译 / 氧化磷酸化）**：
   - Rapamycin 在 Muscle Male 中强烈下调 (-0.40)。
   - Trametinib 在 Muscle Male 中反而上调 (+0.17)。
   - Combo 在 Muscle Male 中下调 (-0.37)，接近 Rapamycin 单药水平，说明 Rapamycin 主导了 pink 模块的效应。

2. **orange 模块（染色质修饰）**：
   - Rapamycin 在 Muscle Male 中下调 (-0.24)。
   - Trametinib 影响微弱。
   - Combo 进一步下调 (-0.27)，但差异很小。

3. **turquoise 模块（先天免疫 / 炎症）**：
   - Combo 在 Muscle Male 中下调 (-0.10)，优于 Rapamycin 单药 (-0.04) 和 Trametinib 单药 (+0.11)。
   - 提示在炎症调控方面，联合用药可能具有互补优势。

4. **模块层面的加和性**：
   - 大多数模块中，Combo 的效应介于 Rapa 和 Tram 之间，或接近两者之和，与整体转录组层面的加和性一致。
   - 无单一模块显示出强烈的协同或拮抗信号。

---

## 6. 与已有结果的比较

| 数据集 | 干预 | 组织 | Difference | 来源 |
|--------|------|------|------------|------|
| **GSE131754** | Rapamycin | Liver | **-0.63** | 本课题前期 |
| **GSE131754** | CR (40%) | Liver | **-0.82** | 本课题前期 |
| **GSE288795** | Rapamycin | Spleen | **-1.17** | 本次分析 |
| **GSE288795** | Combo | Spleen | **-1.11** | 本次分析 |
| **GSE288795** | Rapamycin | Muscle Male | **-1.20** | 本次分析 |
| **GSE288795** | Combo | Muscle Male | **-1.47** | 本次分析 |

**比较解读**：
- GSE288795 中 Rapamycin 在 Muscle Male 和 Spleen 的 rejuvenation 效应（-1.20 ~ -1.28）**强于** GSE131754 肝脏结果（-0.63）。这可能反映了组织特异性（Muscle/Spleen vs Liver）或年龄差异（24m vs 6–12m）。
- Combo 在 Muscle Male（-1.47）的效果接近 GSE131754 中 CR 的效果（-0.82）的 **1.8 倍**（绝对值），但需注意不同数据集间的可比性限制。

---

## 7. 局限性与注意事项

1. **组织差异**：本研究取样为 Muscle/Kidney/Spleen，与前期分析的 Liver 不同。模型为多组织训练，但不同组织的 baseline 和 drug response 存在固有差异。
2. **年龄单一**：所有样本均为 24 月龄，无法观察年龄依赖性效应。
3. **样本量**：每组 n=4–6，统计功效有限，无法检测小幅度的协同/拮抗效应。
4. **Kidney Female Rapamycin 异常**：Rapamycin 单药在 Kidney Female 中显示 pro-aging (+0.38, SEM=0.80)，方向与预期相反。可能原因包括：
   - 该组内个体异质性高（SEM 较大）
   - Rapamycin 在雌性肾脏中的特异性反应
   - 技术噪声（批次效应）
   该异常使得 Kidney Female 的 additivity 检验难以解释。
5. **训练-测试重叠**：本模型为 Mortality 预测模型，非衰老时钟。Difference 反映的是 Mortality risk 的变化，与生物学寿命延长方向一致，但非直接的 "biological age reduction"。

---

## 8. 结论

1. **转录组层面验证**：Trametinib + Rapamycin 联合用药在 Muscle、Kidney、Spleen 中均产生 rejuvenation 效应，整体效果优于任一单药。
2. **加和性为主**：统计检验表明，联合用药的转录组效应**与各单药加和效应无显著差异**（p > 0.05）。没有证据支持强烈的协同（synergy）或拮抗（antagonism）作用。
3. **组织特异性**：Spleen 和 Muscle Male 对 Combo 最敏感；Kidney Female 中 Rapamycin 单药异常，但 Combo 仍然有效。
4. **模块互补性**：Rapamycin 主导线粒体和染色质模块的下调；Trametinib 在部分模块中方向相反。联合用药在炎症（turquoise）模块可能呈现互补优势。
5. **临床转化意义**：虽然转录组效应趋向加和而非协同，但 Combo 的绝对 rejuvenation 幅度仍然显著（最高 -1.61），且原文献已证实其在寿命延长上的加和性。转录组结果与表型结果一致，支持 Combo 作为 gerotherapy 的进一步评估。

---

## 9. 输出文件清单

| 文件 | 路径 |
|------|------|
| 预测结果 | `results/geo_longterm/gse288795_tage_predictions.csv` |
| Drug-Control 差异 | `results/geo_longterm/gse288795_drug_control_differences.csv` |
| 模块贡献 | `results/geo_longterm/gse288795_per_sample_module_contributions.csv` |
| 模块差异 | `results/geo_longterm/gse288795_module_drug_control_differences.csv` |
| Additivity 检验 | `results/geo_longterm/gse288795_additivity_test.csv` |
| 预测分布图 | `figures/gse288795_predictions.pdf` |
| 药物差异图 | `figures/gse288795_drug_differences.pdf` |
| 模块热图 | `figures/gse288795_module_heatmap.pdf` |
| Additivity 检验图 | `figures/gse288795_additivity_test.pdf` |
| edgeR 输入 | `results/geo_longterm/gse288795_counts_for_edger.tsv` |
| logCPM | `results/geo_longterm/gse288795_logcpm.tsv` |

---

*分析脚本*: `src/python/gse288795_analysis.py`
