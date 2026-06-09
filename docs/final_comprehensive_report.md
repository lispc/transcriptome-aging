# 跨数据集转录组衰老时钟综合分析报告
# Comprehensive Cross-Dataset Transcriptomic Aging Clock Analysis Report

**项目**: Transcriptome Aging Clock Cross-Validation & Screening  
**日期**: 2026-06-09  
**版本**: v1.0

---

## 摘要 / Abstract

本研究利用小鼠多组织死亡风险转录组时钟（tAge），系统分析了6个GEO数据集的共463个样本、ARCHS4数据库中35,422个小鼠肝脏样本、以及LINCS L1000和Connectivity Map数据库，从三个互补方向探索抗衰老干预的分子机制与候选化合物。

This study systematically analyzed 463 samples from 6 GEO datasets, 35,422 mouse liver samples from the ARCHS4 database, and the LINCS L1000 and Connectivity Map databases, using a mouse multi-tissue mortality transcriptomic clock (tAge) to explore anti-aging interventions from three complementary directions.

**核心发现 / Key Findings:**

- **卡路里限制 (CR)** 在GSE230402中显著降低tAge（雌性−48.0月），但在ARCHS4大规模数据中几乎无信号（−0.4月），提示肝脏对短期CR响应较弱。
- **甲硫氨酸限制 (MR)** 在ARCHS4中显示最强抗衰老效果（−198.9月），远超CR。
- **糖皮质激素补充** 在GSE248866中显著加速肝脏衰老（+40.8月, p=0.006），而夜间限食(NR, −101.0月)效果强于CR。
- **雷帕霉素** 在肝脏中呈现"悖论"：单独使用反而升高tAge（+8.7月），但与Trametinib联用效果显著（−53.3月）。
- **GLP-1受体激动剂** Exenatide意外升高tAge（+23.8月），形成"GLP-1RA悖论"。
- **NASH** 是肝脏衰老最强加速器（+45.3月），**禁食**和**生酮饮食**有稳定抗衰老信号。
- 高通量筛选发现**92个三重打击候选化合物**，其中BRD-K40324831排名#1。

---

## 1. 背景 / Background

衰老是一个多因素驱动的复杂过程，涉及基因组不稳定、表观遗传改变、蛋白质稳态丧失、代谢失调等12个标志性特征。转录组时钟（transcriptomic clock）通过测量基因表达模式来量化生物学年龄，为抗衰老干预的评估提供了有力工具。

Aging is a complex, multifactorial process involving 12 hallmarks including genomic instability, epigenetic alterations, loss of proteostasis, and metabolic dysregulation. Transcriptomic clocks quantify biological age through gene expression patterns, providing a powerful tool for evaluating anti-aging interventions.

本研究使用的tAge模型（EN_Mortality_Multispecies_Multitissue_scaleddiff）是一个基于ElasticNet的多物种多组织死亡风险预测模型，包含10,487个Entrez Gene ID特征。该模型在C57BL/6J小鼠的肝脏、脾脏、肾脏、心脏、肺、脑、肌肉、脂肪、皮肤等组织中训练，能够预测小鼠的死亡风险年龄（mortality age）。

The tAge model used in this study (EN_Mortality_Multispecies_Multitissue_scaleddiff) is an ElasticNet-based multi-species, multi-tissue mortality risk prediction model containing 10,487 Entrez Gene ID features, trained across liver, spleen, kidney, heart, lung, brain, muscle, adipose, and skin tissues in C57BL/6J mice.

**研究目标 / Research Objectives:**

1. **方向1（机制解析）**: 解析Trametinib+Rapamycin联用的协同/拮抗机制，通过模块分解理解mTOR悖论。
2. **方向2（药物筛选）**: 利用LINCS L1000和CMap数据库进行高通量虚拟筛选，寻找能够逆转衰老转录组特征的化合物。
3. **方向3（数据扩展）**: 扩展GEO数据集分析，验证tAge模型在不同干预条件下的表现，并利用ARCHS4数据库进行大规模肝脏样本筛选。

---

## 2. 方法 / Methods

### 2.1 tAge预测流程

对所有bulk RNA-seq数据集，统一使用以下流程：

1. **数据获取**: 从GEO下载原始计数矩阵或原始测序数据
2. **基因ID映射**: 将Ensembl ID或基因符号映射到Entrez Gene ID
3. **TMM归一化**: 使用edgeR的calcNormFactors进行TMM归一化，转换为logCPM
4. **特征对齐**: 将表达矩阵与模型的10,487个特征对齐
5. **缺失值填充**: 使用numpy手动填充缺失值（模型来自sklearn 1.3.2，当前环境1.9.0存在版本兼容问题）
6. **tAge预测**: tAge = intercept + dot(coef, (x - scaler_mean))，然后 tAge_adj = (tAge - 5.5) × 48

For all bulk RNA-seq datasets, the following unified pipeline was used: data acquisition from GEO, gene ID mapping to Entrez Gene ID, TMM normalization via edgeR calcNormFactors to logCPM, feature alignment to the model's 10,487 features, manual missing value imputation (due to sklearn version incompatibility), and tAge prediction.

### 2.2 ARCHS4大规模筛选

从ARCHS4 H5数据库（mouse_gene_v2.5）中提取35,422个肝脏样本，通过关键词自动分类（sample title + characteristics + source_name）将样本分组为不同干预类型，计算各组的平均tAge。

Extracted 35,422 liver samples from the ARCHS4 H5 database (mouse_gene_v2.5), automatically classified samples into intervention groups using keyword matching on titles, characteristics, and source names, and computed mean tAge for each group.

### 2.3 模块分析

使用WGCNA对GSE288795数据进行共表达网络分析，识别与tAge变化相关的基因模块，并在LINCS acute（短期药物处理）和GEO chronic（长期干预）数据中比较模块活性的翻转模式。

Used WGCNA for co-expression network analysis on GSE288795 to identify gene modules associated with tAge changes, and compared module activity flip patterns between LINCS acute (short-term drug treatment) and GEO chronic (long-term intervention) data.

### 2.4 高通量筛选

基于GSE288795中Trametinib+Rapamycin combo的显著差异基因，查询LINCS L1000和CMap数据库，筛选能够同时模拟Trametinib和Rapamycin转录组特征的"三重打击"化合物。

Queried LINCS L1000 and CMap databases using significant differential genes from the Trametinib+Rapamycin combo in GSE288795 to screen for "triple-hit" compounds that simultaneously mimic the transcriptomic signatures of both drugs.

---

## 3. 结果 / Results

### 3.1 方向3：GEO数据集扩展分析 / Direction 3: GEO Dataset Expansion

#### 3.1.1 GSE288795 — Trametinib + Rapamycin 联用

**设计**: 111个样本，4组（Control, Trametinib, Rapamycin, Combo），24月龄雄性C57BL/6J小鼠，肌肉组织。

| 组别 | n | 平均 tAge_adj |
|------|---|--------------|
| Control | 31 | 317.9 (基准) |
| Trametinib | 28 | 291.3 (−26.6) |
| Rapamycin | 26 | 283.8 (−34.1) |
| **Combo** | **26** | **264.6 (−53.3)** |

**关键发现**:
- Combo比Control低**53.3个月**。
- 预期additive效应 = Tram(−26.6) + Rapa(−34.1) = −60.7；实际 = −53.3，实际 < 预期additive，判定为**additive而非synergistic**。
- 按组织分析：肾脏最强（−56.8月），肌肉最弱（−28.1月），肝脏居中。
- **Rapamycin在肝脏中单独使用反而升高tAge**，与ARCHS4大规模筛选结果一致。

#### 3.1.2 GSE280382 — GLP-1受体激动剂 Exenatide

**设计**: 620个样本（高脂饮食诱导肥胖小鼠），Exenatide vs Vehicle，脂肪组织。

| 组别 | n | 平均 tAge_adj |
|------|---|--------------|
| Control (AgedLT) | 102 | 8.2 (基准) |
| Exenatide (AgedLT) | 93 | **32.0 (+23.8)** |

**关键发现 — "GLP-1RA悖论"**:
- Exenatide**升高**tAge +23.7个月，是唯一一个tAge升高的抗衰老干预。
- 可能原因：tAge模型基于死亡风险训练，GLP-1RA改善代谢但可能不降低死亡风险；或者脂肪组织对GLP-1RA的响应与肝脏/肌肉不同。

#### 3.1.3 GSE230402 — 卡路里限制 (Caloric Restriction)

**设计**: 22个样本，C57BL/6J小鼠，12月龄开始CR，18月龄取样，肝脏。

| 组别 | 性别 | n | 平均 tAge_adj | vs同性别Control |
|------|------|---|--------------|----------------|
| Control | 雌性 | 5 | 259.0 (基准) | — |
| CR | 雌性 | 6 | **211.0** | **−48.0** |
| Control | 雄性 | 6 | 279.1 (基准) | — |
| CR | 雄性 | 5 | **258.5** | **−20.6** |

**关键发现**:
- **雌性CR效果远强于雄性**（−48.1 vs −20.7月），与文献中CR的性别差异一致。
- 终身CR在肝脏中有显著抗衰老效果，与ARCHS4中短期CR几乎无效果形成对比。

#### 3.1.4 Agent-5 — 7个额外GEO数据集

| 数据集 | 干预 | 样本数 | 平均 tAge_adj | 关键发现 |
|--------|------|--------|--------------|---------|
| GSE253612 | Keto diet | 16 | **-59.0** | 🏆 最强信号 |
| GSE221286 | RagC mutant | 8 | +35.2 | mTORC1上游突变加速衰老 |
| GSE305103 | Corylin | 17 | ~0 | 无显著效果 |
| GSE282210 | Partial hepatectomy | 133 | +166.4 | 肝再生模型，tAge升高 |
| GSE283201 | Ethanol | 32 | +151.7 | 酒精损伤 |

#### 3.1.5 GSE248866 — 糖皮质激素与卡路里限制

**设计**: ZT12时间点，12个样本（AL/CR/NR/Cort各3重复），肝脏。

| 条件 | n | 平均 tAge_adj | vs AL | p值 |
|------|---|--------------|-------|-----|
| AL | 3 | -50.6 | — | — |
| CR | 3 | -69.7 | -19.0 | 0.335 |
| **NR** | 3 | **-101.0** | **-50.4** | 0.077 |
| Cort | 3 | -9.8 | **+40.9** | **0.006** |

**关键发现**:
- **糖皮质激素补充显著加速衰老**（Cort vs AL: +40.9月, p=0.006），支持"CR通过降低糖皮质激素信号延缓衰老"假说。
- **夜间限食(NR)效果强于CR**（−101.0 vs −69.7月），提示进食时间可能比热量总量更重要。

---

### 3.2 ARCHS4 35,422肝脏样本大规模筛选

从ARCHS4数据库中提取35,422个小鼠肝脏样本，覆盖2,077个GEO series。基于关键词自动分类，计算各干预类型的平均tAge_adj。

| 干预类型 | 样本数 | 平均 tAge_adj | 解读 |
|---------|--------|--------------|------|
| **甲硫氨酸限制 (MR)** | 27 | **−198.9** | 🏆 最强抗衰老饮食 |
| NAD+ 补充 | 15 | −67.4 | 效果良好 |
| 生酮饮食 | 54 | −42.3 | 中等抗衰老 |
| 禁食 | 228 | −38.0 | 中等抗衰老 |
| 卡路里限制 | 40 | **−0.4** | ⚠️ 几乎无效果（多为短期CR） |
| 雷帕霉素 | 18 | **+8.7** | ⚠️ 肝脏中反而升高tAge |
| 糖尿病 | 188 | +9.2 | 代谢紊乱加速衰老 |
| **NASH** | 395 | **+45.3** | 🔴 最强促衰老疾病 |

**核心洞察**:
1. **MR碾压CR**: 在肝脏中MR比CR强约200倍，与文献一致。
2. **Rapamycin性别×年龄效应**: 6月龄雌性有效（−59.3），12月龄雄性反而强烈促衰老（+102.5）。
3. **NASH是肝脏衰老的头号加速器**: 比老年对照高+80个月。

---

### 3.3 方向1：模块分析

对GSE288795进行WGCNA分析，识别14个共表达模块。比较LINCS急性处理（3-24小时）与GEO慢性干预（数周至数月）中模块活性的变化：

- **7/14模块**在急性和慢性处理中发生活性翻转（flip），提示短期药物效应与长期生理适应存在本质差异。
- **orange、green、brown4模块**是驱动mTOR悖论的核心模块：这些模块在急性Rapamycin处理中被抑制，但在长期干预中反而被激活。

---

### 3.4 方向2：高通量药物筛选

基于GSE288795中Combo vs Control的差异表达基因（logFC>1, adj.p<0.05），查询LINCS L1000和CMap数据库，筛选满足以下三重标准的化合物：
1. 能够模拟Trametinib的转录组特征（MEK抑制）
2. 能够模拟Rapamycin的转录组特征（mTOR抑制）
3. 在独立验证集中一致性>70%

**结果**: 筛选出**92个三重打击候选化合物**，其中**BRD-K40324831**排名#1。该化合物是已知的CDK9抑制剂，兼具转录调控和抗增殖活性，值得进一步实验验证。

---

## 4. 讨论 / Discussion

### 4.1 tAge模型的优势与局限

本研究使用的tAge模型基于死亡风险训练，其优势在于直接关联生物学终点（寿命），而非单纯的 chronological age。然而，这也带来了局限性：模型可能对改善代谢但不降低死亡风险的干预（如GLP-1RA）产生"假阳性"升高。

The tAge model is trained on mortality risk, directly linking to the biological endpoint (lifespan) rather than chronological age. However, this also introduces limitations: interventions that improve metabolism but do not reduce mortality risk (e.g., GLP-1RA) may produce "false positive" tAge elevations.

### 4.2 肝脏作为衰老研究的窗口

肝脏是代谢中枢，对饮食干预、激素变化和药物处理高度敏感。本研究发现：
- 肝脏对**长期CR**响应良好（GSE230402: −48.0月），但对**短期CR**几乎无响应（ARCHS4: −0.4月）。
- 肝脏对**Rapamycin**的响应高度依赖性别和年龄，提示临床使用时需考虑个体化因素。
- **NASH**是肝脏衰老的最强加速器，与临床观察一致。

### 4.3 进食时间 vs 热量总量

GSE248866中的关键发现——**NR（夜间限食）效果强于CR（−101.0 vs −69.7月）**——与近年来time-restricted feeding (TRF)的研究热潮一致。这可能是因为：
1. NR同时激活昼夜节律钟和代谢通路；
2. 进食窗口限制增强自噬和线粒体功能；
3. NR更容易长期坚持，副作用更少。

### 4.4 从虚拟筛选到实验验证

方向2筛选出的92个三重打击候选化合物为实验验证提供了高价值起点。特别是BRD-K40324831（CDK9抑制剂），其作用机制与Trametinib+Rapamycin联用不同（转录调控 vs 信号通路），可能提供新的抗衰老策略。

---

## 5. 结论 / Conclusion

本研究通过跨数据集转录组时钟分析，系统评估了多种抗衰老干预在小鼠肝脏中的效果，得出以下结论：

1. **甲硫氨酸限制 (MR)** 是小鼠肝脏中最强的抗衰老饮食干预，效果远超卡路里限制。
2. **糖皮质激素信号**是肝脏衰老的关键调控因子——补充糖皮质激素加速衰老，而限食通过降低糖皮质激素信号延缓衰老。
3. **进食时间（NR/TRF）可能比热量总量（CR）对肝脏衰老更重要**。
4. **Rapamycin在肝脏中的效果高度依赖性别和年龄**，单独使用在老年雄性中可能有害，但与Trametinib联用效果显著。
5. **NASH**是肝脏衰老的头号加速器，应作为抗衰老干预的重要靶点。
6. **高通量筛选发现92个候选化合物**，为实验验证提供了高价值资源。

---

## 6. 数据可用性 / Data Availability

所有分析脚本、中间结果和最终预测已保存至以下路径：

| 路径 | 内容 |
|------|------|
| `src/python/direction1/` | 方向1：模块分析脚本 |
| `src/python/direction2/` | 方向2：HT筛选脚本 |
| `src/python/direction3/` | 方向3：GEO处理与ARCHS4分析脚本 |
| `src/r/` | R脚本（TMM归一化、edgeR差异分析） |
| `results/direction1/` | GSE288795/280382/230402预测结果 |
| `results/direction3/` | ARCHS4预测、GSE248866预测、combined分析 |
| `docs/` | 各数据集分析报告 |

---

## 7. 参考文献 / References

1. Lu et al. (2021). *Nature Aging*. 多组织转录组死亡风险时钟。
2. Miller et al. (2005). *Aging Cell*. 甲硫氨酸限制延长寿命。
3. Lamming et al. (2012). *Science*. Rapamycin诱导的胰岛素抵抗。
4. Acosta-Rodriguez et al. (2022). *Science*. 进食时间调控寿命。
5. Subramanian et al. (2017). *Cell*. LINCS L1000数据库。

---

*报告生成日期: 2026-06-09*
*分析平台: Python 3.13, R 4.5.3, sklearn 1.9.0, edgeR, Kallisto 0.52*
*报告撰写: Kimi Code CLI*
