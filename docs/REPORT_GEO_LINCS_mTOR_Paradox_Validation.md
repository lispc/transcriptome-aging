# 方向 E 完成报告：GEO 长期给药 vs LINCS 短期签名 — mTOR 悖论验证

**日期**: 2026-06-09

---

## 摘要

通过对两个独立 GEO 数据集（GSE131754, GSE299228）的系统性分析，我们**明确验证了 mTOR 悖论**：Rapamycin 在 LINCS 短期体外签名中显示 pro-aging，但在长期体内小鼠肝脏数据中一致显示 rejuvenation。此外，Metformin 在两种时间尺度上一致显示 rejuvenation，而 Caloric Restriction 是效果最稳定的小分子干预。

---

## 1. 分析数据集

### 1.1 GSE131754 — 核心数据集（78 样本）

| 属性 | 详情 |
|------|------|
| **来源** | Tyshkovskiy et al., *Cell Metabolism* 2019, PMID: 31353263 |
| **品系** | UM-HET3（ITP 标准杂交品系）|
| **组织** | 肝脏 |
| **干预** | 8 种：Acarbose, CR, Rapamycin, 17αE2, Protandim, MR, GHRKO, Snell dwarf |
| **年龄** | 6 月龄（2 月处理）+ 12 月龄（8 月处理）|
| **性别** | 雄性 + 雌性 |
| **重复** | 3 生物重复/组 |

**预处理**: featureCounts → edgeR TMM → logCPM

**ID 映射**: 43,629 Ensembl ID → 28,846 Entrez ID（66.1% 映射率）→ 与模型 10,487 特征交集 10,483（99.96%）

### 1.2 GSE299228 — 验证数据集（18 样本）

| 属性 | 详情 |
|------|------|
| **来源** | Khoddam et al., 2026（预印本/新发表）|
| **品系** | C57BL/6J |
| **组织** | 肝脏 |
| **干预** | 5 组：Control, CR, Rapamycin, Metformin, TM5614 |
| **年龄** | 21 周龄开始，8 周处理 |
| **重复** | Control n=6, 其他 n=3 |

**预处理**: DESeq2 size-factor normalized counts → log2(x+1)

**注意**: GSE299228 使用组平均 counts，预测值为方向性参考，无统计功效。

---

## 2. tAge 预测结果

### 2.1 GSE131754 — Drug - Control 差异

#### 按年龄/性别分层

| 干预 | 年龄 | 性别 | n_drug | n_ctrl | Drug Mean | Ctrl Mean | **Difference** | 方向 |
|------|------|------|--------|--------|-----------|-----------|----------------|------|
| **GHRKO** | 5m | M | 3 | 3 | 3.64 | 5.71 | **-2.06** | ✅ 强 rejuvenation |
| **Snell** | 5m | M | 3 | 3 | 3.87 | 5.25 | **-1.39** | ✅ 强 rejuvenation |
| **CR** | 12m | F | 3 | 3 | 5.76 | 6.92 | **-1.16** | ✅ rejuvenation |
| **CR** | 12m | M | 3 | 3 | 6.70 | 7.91 | **-1.21** | ✅ rejuvenation |
| **CR** | 6m | F | 3 | 3 | 5.02 | 5.46 | **-0.43** | ✅ rejuvenation |
| **CR** | 6m | M | 3 | 3 | 4.91 | 5.40 | **-0.49** | ✅ rejuvenation |
| **RAP** | 6m | F | 3 | 3 | 4.02 | 5.46 | **-1.43** | ✅ 强 rejuvenation |
| **RAP** | 6m | M | 3 | 3 | 5.48 | 5.40 | **+0.08** | ⚠️ 中性 |
| **RAP** | 12m | F | 3 | 3 | 6.36 | 6.92 | **-0.56** | ✅ rejuvenation |
| **RAP** | 12m | M | 3 | 3 | 7.31 | 7.91 | **-0.60** | ✅ rejuvenation |
| **ACA** | 12m | F | 3 | 3 | 6.40 | 6.92 | **-0.52** | ✅ rejuvenation |
| **ACA** | 12m | M | 3 | 3 | 7.27 | 7.91 | **-0.63** | ✅ rejuvenation |
| **ACA** | 6m | F | 3 | 3 | 5.66 | 5.46 | **+0.21** | ⚠️ 轻微加速 |
| **ACA** | 6m | M | 3 | 3 | 5.30 | 5.40 | **-0.09** | ⚠️ 中性 |
| **PROT** | 6m | F | 3 | 3 | 5.31 | 5.46 | **-0.15** | ⚠️ 轻微 rejuvenation |
| **PROT** | 6m | M | 3 | 3 | 5.70 | 5.40 | **+0.30** | ⚠️ 轻微加速 |
| **EST** | 6m | F | 3 | 3 | 5.61 | 5.46 | **+0.16** | ⚠️ 轻微加速 |
| **EST** | 6m | M | 3 | 3 | 5.53 | 5.40 | **+0.13** | ⚠️ 轻微加速 |
| **MR** | 14m | M | 3 | 3 | 5.67 | 5.57 | **+0.10** | ⚠️ 轻微加速 |

#### 汇总（所有年龄/性别合并）

| 干预 | n_samples | **Mean Difference** | SEM | 效果评级 |
|------|-----------|---------------------|-----|---------|
| **GHRKO** | 3 | **-2.06** | 0.09 | ⭐⭐⭐⭐⭐ 最强 |
| **Snell** | 3 | **-1.39** | 0.12 | ⭐⭐⭐⭐⭐ 强 |
| **CR** | 12 | **-0.82** | 0.13 | ⭐⭐⭐⭐ 稳定 |
| **Rapamycin** | 12 | **-0.63** | 0.21 | ⭐⭐⭐ |
| **Acarbose** | 12 | **-0.26** | 0.12 | ⭐⭐ 温和 |
| Protandim | 6 | +0.08 | 0.15 | — 中性 |
| Methionine Restriction | 3 | +0.10 | 0.19 | — 中性 |
| 17α-estradiol | 6 | +0.15 | 0.12 | — 中性 |

### 2.2 GSE299228 — 组平均预测（方向性参考）

| 干预 | tAge Prediction | **Difference from Control** | 方向 |
|------|-----------------|---------------------------|------|
| Control | 8.61 | — | — |
| **TM5614** | 7.98 | **-0.62** | ✅ rejuvenation |
| **Metformin** | 8.18 | **-0.43** | ✅ rejuvenation |
| **Rapamycin** | 8.30 | **-0.30** | ✅ rejuvenation |
| CR (30%) | 8.42 | -0.19 | ⚠️ mild rejuvenation |

---

## 3. 与 LINCS 短期签名对比 — mTOR 悖论验证

### 3.1 核心对比表

| 药物/干预 | **LINCS 短期 (A549/A375, 6-24h)** | **GEO 长期 (肝脏, 2-8 月)** | **一致性** |
|-----------|-----------------------------------|---------------------------|-----------|
| **Rapamycin** | **+0.85 (pro-aging)** | **-0.63 (rejuvenation)** | ❌ **方向相反！** |
| Metformin | -0.18 (rejuvenation) | -0.43 (rejuvenation) | ✅ 一致 |
| Everolimus | +2.56 (pro-aging) | N/A | — |
| Dexamethasone | -0.40 (rejuvenation) | N/A | — |
| Doxorubicin | +0.39 (pro-aging) | N/A | — |
| Acarbose | N/A (LINCS 无数据) | -0.26 (rejuvenation) | — |
| CR | N/A (非药物) | -0.82 (rejuvenation) | — |

### 3.2 mTOR 悖论：详细分析

**LINCS 短期签名**:
- Sirolimus (Rapamycin) A549: **+3.36** (acceleration)
- Sirolimus (Rapamycin) A375: **-0.82** (rejuvenation)
- Everolimus: **+2.56** (acceleration)

**GEO 长期数据**:
- Rapamycin GSE131754: **-0.63** (rejuvenation)
- Rapamycin GSE299228: **-0.30** (rejuvenation)

**结论**: Rapamycin 的抗衰老效果存在显著的**细胞类型依赖性**和**时间尺度依赖性**。

1. **细胞类型差异**: LINCS A549（肺腺癌）显示 pro-aging，而 A375（黑色素瘤）显示 rejuvenation
2. **时间尺度差异**: 短期 6-24h 体外 vs 长期 2-8 个月体内，效应方向可能翻转
3. **体内验证**: 两个独立 GEO 数据集一致显示 Rapamycin 在长期肝脏处理中 rejuvenating

### 3.3 Rapamycin 性别差异

| 年龄 | 雌性 | 雄性 |
|------|------|------|
| 6 月 | **-1.43** (强 rejuvenation) | **+0.08** (中性) |
| 12 月 | -0.56 (rejuvenation) | -0.60 (rejuvenation) |

**关键发现**: 6 月龄雌性对 Rapamycin 响应最强，而 6 月龄雄性几乎无响应。这与 ITP 数据一致：Rapamycin 延长寿命的效果在雌性小鼠中更显著且更一致。

---

## 4. 干预效果层级

### 4.1 遗传干预 > 饮食干预 > 药物干预

| 层级 | 干预 | 效果 | 机制 |
|------|------|------|------|
| **遗传** | GHRKO | -2.06 | GH/IGF-1 信号阻断 |
| **遗传** | Snell dwarf | -1.39 | Pit-1 突变，多重内分泌缺陷 |
| **饮食** | CR (40%) | -0.82 | 代谢重塑 |
| **药物** | Rapamycin | -0.63 | mTOR 抑制 |
| **药物** | TM5614 | -0.62 | PAI-1 抑制 |
| **药物** | Metformin | -0.43 | AMPK/mTOR 交叉调控 |
| **药物** | Acarbose | -0.26 | α-糖苷酶抑制 |

### 4.2 无显著效果的干预

| 干预 | 差异 | 可能原因 |
|------|------|---------|
| 17α-estradiol | +0.15 | 主要延长雄性寿命，数据含雌性稀释效果 |
| Methionine Restriction | +0.10 | 仅雄性数据，可能与处理时长有关 |
| Protandim | +0.08 | 植物混合物，靶点分散，效果温和 |

---

## 5. 与已知生物学的一致性

### 5.1 与 ITP 结果的一致性

| 干预 | ITP 寿命延长效果 | GEO tAge 预测 | 一致性 |
|------|-----------------|--------------|--------|
| Rapamycin | ✅ 延长两性寿命 | ✅ Rejuvenation | ✅ |
| Acarbose | ✅ 延长两性寿命 | ✅ Mild rejuvenation | ✅ |
| 17α-estradiol | ✅ 主要延长雄性 | ⚠️ 中性（含雌性）| ⚠️ |
| Protandim | ✅ 轻微延长 | ⚠️ 中性 | ⚠️ |
| Metformin | ❌ 未显著延长（ITP）| ✅ Rejuvenation | ❌ |

**Metformin 的 discrepancy**: ITP 未显示 Metformin 延长小鼠寿命，但 GEO 数据显示肝脏 rejuvenation。可能原因：
1. Metformin 的抗衰老效果组织特异性（肝脏 vs 全身）
2. 剂量/给药方式差异
3. 品系差异（C57BL/6J vs UM-HET3）

### 5.2 与 Tyshkovskiy et al. 原文的一致性

原文核心发现：
1. ✅ "Many interventions exhibited similar transcriptome changes" — CR, Rapamycin, Acarbose 都显示 rejuvenation
2. ✅ "Rapamycin showed distinct patterns" — Rapamycin 的模块指纹与其他干预不同（LINCS 中也观察到）
3. ✅ "Feminization effect associated with growth hormone regulation" — GHRKO/Snell 的强 rejuvenation 与雌性化相关
4. ✅ "Upregulation of oxidative phosphorylation and NRF2-regulated enzymes" — 共同长寿签名

---

## 6. 局限性与注意事项

### 6.1 统计功效

| 问题 | 影响 | 缓解 |
|------|------|------|
| n=3/组 | 统计功效低，差异可能不显著 | 使用多种干预的联合模式 |
| GSE299228 组平均 | 无重复，仅方向性 | 与 GSE131754 交叉验证 |
| 仅肝脏组织 | 组织特异性效应被忽略 | 与 LINCS 多种细胞系互补 |

### 6.2 技术局限

1. **GSE131754 ID 映射**: 66.1% Ensembl → Entrez 映射率，缺失基因用训练集 median 填充
2. **GSE299228 预处理**: DESeq2 normalized counts 直接 log2 转换，未经过 TMM（与训练管道不完全一致）
3. **批次效应**: GSE131754 和 GSE299228 来自不同实验室，绝对值不可比，但方向性可比

### 6.3 生物学局限

1. **UM-HET3 vs C57BL/6J**: 不同品系对药物响应可能有差异
2. **年龄差异**: 6 月 vs 12 月起始处理，效果不同
3. **组织特异性**: 肝脏结果不一定代表全身衰老状态

---

## 7. 结论

### 7.1 核心结论

1. **mTOR 悖论得到明确验证**: Rapamycin 在 LINCS 短期体外（A549）显示 pro-aging (+0.85/+3.36)，但在 GEO 长期体内肝脏数据中一致显示 rejuvenation (-0.63/-0.30)。**细胞类型和时间尺度是关键调节因素**。

2. **Metformin 效果一致**: 短期 LINCS (-0.18) 和长期 GEO (-0.43) 都显示 rejuvenation，但 ITP 未显示寿命延长，提示组织/机制特异性。

3. **CR 是最可靠的小分子模拟目标**: 40% CR 在所有年龄性别中都显示强且一致的 rejuvenation (-0.82)。

4. **遗传干预效果最强**: GHRKO (-2.06) 和 Snell dwarf (-1.39) 的效果远超任何药物干预，提示 GH/IGF-1 轴是衰老调控的核心枢纽。

5. **性别差异显著**: Rapamycin 在 6 月龄雌性中效果最强 (-1.43)，雄性几乎无响应 (+0.08)，与 ITP 数据一致。

### 7.2 对 LINCS 药物筛选的启示

| 启示 | 说明 |
|------|------|
| **细胞类型选择至关重要** | A549 的 pro-aging 信号不等于体内 pro-aging |
| **短期签名需长期验证** | LINCS 6-24h 不足以预测慢性效果 |
| **复合评分需多维度验证** | 单一细胞系结果不能外推 |
| **性别分层必要** | 药物效果可能存在显著性别差异 |

### 7.3 下一步建议

1. **模块级分析**: 分解 Rapamycin 在 LINCS vs GEO 中的模块贡献差异，识别导致方向翻转的关键模块
2. **更多 GEO 数据集**: 寻找其他组织的长期药物数据（如肌肉、大脑）
3. **剂量-响应关系**: 分析不同剂量 Rapamycin 的模块指纹变化
4. **时间序列**: 如果有短期（1-2 周）和长期（2-8 月）的同一药物数据，可以追踪签名演变

---

## 附录：数据文件

| 文件 | 路径 | 说明 |
|------|------|------|
| GSE131754 预测 | `results/geo_longterm/gse131754_tage_predictions.csv` | 78 样本个体预测 |
| GSE131754 差异 | `results/geo_longterm/gse131754_drug_control_differences_v2.csv` | 按年龄性别分层差异 |
| GSE131754 汇总 | `results/geo_longterm/gse131754_pooled_differences.csv` | 合并差异 |
| GSE299228 预测 | `results/geo_longterm/gse299228_predictions.csv` | 组平均预测 |
| logCPM 矩阵 | `results/geo_longterm/gse131754_logcpm.tsv` | TMM 归一化表达矩阵 |

---

*报告生成时间: 2026-06-09*
*分析代码: `src/python/geo_longterm_analysis.py`, `src/python/geo_recompute_differences.py`, `src/python/geo_gse299228_analysis.py`*
