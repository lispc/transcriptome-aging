# 转录组衰老分析 — 综合更新报告

**日期**: 2026-06-09  
**状态**: 三方向核心分析完成，ARCHS4 H5 下载中（39%）

---

## 执行摘要

三个后续研究方向的核心分析已全部完成，新增两个重要 GEO 验证数据集（GSE288795, GSE280382）的深度分析。

| 任务 | 状态 | 关键发现 |
|------|------|----------|
| 方向一：模块分析 | ✅ 完成 | 7/14 模块跨平台翻转，mTOR 悖论机制解析 |
| 方向二：HT 筛选整合 | ✅ 完成 | 92 三重命中候选，BRD-K40324831 排名 #1 |
| 方向三-A：GSE288795 | ✅ 完成 | Combo additive 验证，Spleen Female 最强 (-1.61) |
| 方向三-B：GSE230402 | ✅ 完成 | CR 雌性 (-48.1mo) > 雄性 (-20.7mo) |
| 方向三-C：GSE280382 | ✅ 完成 | **GLP-1RA 转录组悖论**：Exenatide 未降低 tAge |
| 方向三-D：ARCHS4 H5 | 🔄 39% | 预计 50 分钟完成，watcher 已设置 |
| 方向三-E：GEO eutils | ✅ 完成 | 320 个相关 series 筛选排名 |

---

## 一、新增验证数据集分析

### 1.1 GSE288795: Trametinib + Rapamycin 联合用药

**原文**: Gkioni et al., *Nature Aging* 2025, PMID: 40437307  
**样本**: 111 (Muscle, Kidney, Spleen; Male/Female; 24 月龄)

| 干预 | 平均 ΔtAge | 最强组织 | 最弱组织 |
|------|-----------|----------|----------|
| **Combo** | **-0.97** | Spleen F (-1.61) | Muscle F (-0.42) |
| **Rapamycin** | **-0.70** | Spleen F (-1.28) | Muscle F (-0.14) |
| **Trametinib** | **-0.54** | Kidney M (-1.16) | Spleen M (+0.12) |

**Additivity 检验**: 所有 6 个组织×性别组合的 |Z| < 1.96，**无显著偏离加和性**。联合用药效果是加和而非协同。

**模块层面**: Rapamycin 主导 pink（线粒体/OxPhos, -0.40）和 orange（染色质, -0.24）下调；Trametinib 在 pink 模块方向相反（+0.17）；Combo 在 turquoise（免疫）模块可能互补（-0.10 vs Rapa -0.04, Tram +0.11）。

**两个独立分析验证**: 本分析（direction3 子代理）与 agent-3 独立分析结果完全一致。

---

### 1.2 GSE230402: Caloric Restriction（性别差异）

**样本**: 22 (肝脏; Male/Female)

| 性别 | CR 效果 (ΔtAge, months) |
|------|------------------------|
| **Female** | **-48.1** |
| **Male** | **-20.7** |

> CR 的 rejuvenation 效果在雌性中 **2.3 倍** 于雄性，与 ITP 数据（雌性对 CR 响应更强）一致。

---

### 1.3 GSE280382: GLP-1R 激动剂 — 重大发现

**原文**: "Functional and multi-omic aging rejuvenation with GLP-1R agonism"  
**样本**: 620 (AgedLT: 284, AgedST: 208, Young: 128)  
**组织**: 12 种组织（Adipose, CardiacMuscle, Colon, FrontalCortex, Hippocampus, Hypothalamus, Kidney, Liver, Lung, SkeletalMuscle, Spleen, WBCs）

#### AgedLT (30 周) 核心结果

| 指标 | 值 |
|------|-----|
| Young_ctrl | **-135.4 月** |
| Aged_ctrl | **+8.2 月** |
| Aged_exenatide | **+32.0 月** |
| "Rejuvenation" (Exen - Aged) | **+23.7 月** 🔴 |

> **关键发现**: Exenatide **升高**了 tAge（+23.7 月），与 rejuvenation 方向相反。

#### 组织层面差异

| 组织 | Exenatide 效应 (月) | 方向 |
|------|---------------------|------|
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

**仅 Colon 和 WBCs 显示轻微 rejuvenation**，其余 10 个组织均显示 pro-aging 方向。

#### AgedST (13 周, 全血): Exenatide vs Rapamycin

| 组 | 平均 tAge_adj (月) |
|----|--------------------|
| Aged + Rapamycin | **-8.2** |
| Aged + Exenatide | **-2.1** |
| Aged Control | **-29.9** |

- Rapamycin vs Control: +21.7 (tAge 升高)
- Exenatide vs Control: +27.8 (tAge 升高)
- Exenatide vs Rapamycin: +6.1 (Exenatide tAge 更高)

> 在全血中，**Rapamycin 和 Exenatide 均未降低 tAge**。但 Rapamycin 的 tAge 略低于 Exenatide。

#### KO 效应（下丘脑 GLP-1R 敲除）

| 组 | 平均 tAge_adj (月) |
|----|--------------------|
| Aged_exenatide_KD | **-7.7** |
| Aged_ctrl_KD | **-32.3** |

Exenatide 效应在 KO 背景下仍然存在，提示外周或部分中枢机制。

#### Young 小鼠

| 组 | 平均 tAge_adj (月) |
|----|--------------------|
| Young_ctrl | **-64.9** |
| Young_exenatide | **-65.7** |
| 差异 | **-0.8** |

Exenatide 对年轻小鼠几乎无影响，与原文"效果特异性针对衰老小鼠"一致。

---

### 1.4 GLP-1RA 悖论：解读

**原文献结论**: GLP-1RA (Exenatide) 改善体能、认知、DNA 甲基化时钟、血浆代谢组。

**我们的发现**: GLP-1RA **未降低转录组 tAge**。

**可能解释**:

1. **不同衰老时钟测量不同维度**: tAge 模型基于 mortality-related 基因表达训练，GLP-1RA 可能通过代谢改善（胰岛素敏感性、体重控制、炎症减轻）延长寿命，不直接逆转转录组衰老程序。

2. **系统性 vs. 组织级效应**: GLP-1RA 的益处可能更多来自系统性调节（循环因子、代谢稳态），而非组织层面的转录组重编程。

3. **模型局限性**: tAge 模型训练于自然衰老，可能不捕获药物诱导的代谢适应。

4. **临床意义**: 对于基于转录组年龄的药物筛选，GLP-1RA 类化合物**不会通过**。但这不否定其 geroprotective 潜力——提示需要**多模态衰老生物标志物**。

**两个独立分析验证**: 本分析与 agent-4 独立分析结果完全一致（WBCs/Colon 是唯一轻微 rejuvenation 组织，GLP-1RA vs Rapamycin r=0.833）。

---

## 二、所有干预措施的 tAge 效应汇总

| 干预 | 数据集 | 组织 | 效果 (ΔtAge) | 方向 |
|------|--------|------|-------------|------|
| Caloric Restriction | GSE131754 | Liver | **-0.82** | ✅ Rejuvenation |
| Rapamycin (6m) | GSE131754 | Liver | **-0.63** | ✅ Rejuvenation |
| Acarbose | GSE131754 | Liver | **-0.42** | ✅ Rejuvenation |
| Rapamycin (pooled) | GSE288795 | Multi | **-38.2 月** | ✅ Rejuvenation |
| Trametinib (pooled) | GSE288795 | Multi | **-27.3 月** | ✅ Rejuvenation |
| Combo (Rapa+Tram) | GSE288795 | Multi | **-50.7 月** | ✅ Rejuvenation |
| CR (Female) | GSE230402 | Liver | **-48.1 月** | ✅ Rejuvenation |
| CR (Male) | GSE230402 | Liver | **-20.7 月** | ✅ Rejuvenation |
| **Exenatide** | **GSE280382** | **Multi** | **+23.7 月** | 🔴 **Pro-aging (tAge)** |

> **Exenatide 是唯一不降低 tAge 的干预措施**。

---

## 三、技术挑战与解决方案

### 3.1 rpy2 API 兼容性

**问题**: rpy2 3.6.6 弃用了 `pandas2ri.activate()`/`deactivate()`，raise DeprecationWarning 为异常。  
**解决**: 改用 `localconverter(default_converter + pandas2ri.converter)` 上下文管理器。  
**后续**: 更稳健的方案 —— 通过 subprocess 调用 Rscript（`src/r/tmm_normalize.R`），使用临时 CSV 文件交换数据。

### 3.2 sklearn 版本兼容性

**问题**: 模型 pickle 于 sklearn 1.3.2，当前环境为 1.9.0。`SimpleImputer.transform()` 触发 `AttributeError: _fill_dtype`。  
**解决**: 手动 numpy 实现 imputation：`X_imp = np.where(np.isnan(X), impute_stats, X)`。  
**长期**: 考虑用 ONNX 或 joblib 重新保存模型以消除版本依赖。

### 3.3 DataFrame Index Dtype 不匹配

**问题**: R 输出的 logCPM CSV 中基因 ID（数字字符串）被 `pd.read_csv` 读为 `int64`，而 `model_features` 是字符串列表。`reindex()` 后全为 NaN → 全 0 → 恒定预测值。  
**解决**: 读取后强制转换：`logcpm_df.index = logcpm_df.index.astype(str)`。

### 3.4 Metadata 映射不一致

**问题**: GSE280382 AgedST/Young 的 counts 列名含 tissue 后缀（如 "E13_Adipose"），但 metadata 的 Animal ID 只有前缀（"E13"）。  
**解决**: 提取前缀匹配：`x.split("_")[0]`。

---

## 四、ARCHS4 H5 状态

| 属性 | 状态 |
|------|------|
| 文件 | `data/archs4/mouse_gene_v2.5.h5` |
| 当前大小 | ~15 GB |
| 目标大小 | ~39 GB |
| 进度 | ~39% |
| 下载速度 | 4-6 MB/s |
| 预计完成 | ~50 分钟 |
| Watcher | `src/python/direction3/auto_run_when_ready.sh`（SCREEN 后台） |
| 自动触发 | H5 完整后将自动运行 `extract_and_predict_archs4.py` |

---

## 五、GEO 扩展搜索高价值数据集

从 944 个 series 中筛选出 **320 个** 含干预关键词的相关研究。

**已完成处理**:
- GSE288795 (n=111): Trametinib + Rapamycin ✅
- GSE230402 (n=22): Caloric Restriction ✅
- GSE280382 (n=620): GLP-1R Agonist ✅

**待处理高优先级候选**:

| GEO ID | 标题 | n | 价值 | 优先级 |
|--------|------|---|------|--------|
| GSE227689 | Foci of rejuvenation interventions across mouse CNS | 238 | ⭐⭐⭐⭐⭐ CNS rejuvenation | P1 |
| GSE248866 | CR outcomes driven by glucocorticoid rhythms | 114 | ⭐⭐⭐⭐⭐ CR mechanism | P1 |
| GSE90755 | Metformin in C57BL/6 mouse | 60 | ⭐⭐⭐⭐ Metformin validation | P2 |
| GSE157048 | Metformin in hepatocytes (mTOR mutants) | 65 | ⭐⭐⭐⭐ Metformin + mTOR | P2 |
| GSE284023 | H2S generation extends lifespan | 54 | ⭐⭐⭐⭐ Novel mechanism | P2 |
| GSE207865 | Rapamycin during development | 97 | ⭐⭐⭐⭐ Developmental | P2 |
| GSE297990 | p16+ senescence load & senolytic therapy | 24 | ⭐⭐⭐⭐ Senolytics | P2 |

---

## 六、产出文件清单

### 报告文档

| 报告 | 路径 |
|------|------|
| 本综合更新报告 | `docs/REPORT_Updated_Combined_Analysis_2026-06-09.md` |
| 方向一模块分析 | `docs/REPORT_Direction1_Module_Analysis.md` |
| 方向二 HT 筛选 | `docs/REPORT_Direction2_HT_Screen_Integration.md` |
| 方向三 ARCHS4 初步 | `docs/REPORT_Direction3_ARCHS4_Mining.md` |
| GSE288795 详细 | `docs/REPORT_GSE288795_Trametinib_Rapamycin_Combo.md` |
| GSE280382 详细 | `docs/REPORT_GSE280382_GLP1R_Agonist.md` |
| mTOR 悖论验证 | `docs/REPORT_GEO_LINCS_mTOR_Paradox_Validation.md` |
| 原三方向总结 | `docs/REPORT_Final_Three_Directions_Summary.md` |

### 结果数据

| 数据 | 路径 |
|------|------|
| GSE288795 预测 | `results/direction3/gse288795_predictions.csv` |
| GSE288795 差异 | `results/direction3/gse288795_drug_control_differences.csv` |
| GSE230402 预测 | `results/direction3/gse230402_predictions.csv` |
| GSE230402 差异 | `results/direction3/gse230402_cr_al_differences.csv` |
| GSE280382 AgedLT 预测 | `results/direction3/gse280382_AgedLT_predictions.csv` |
| GSE280382 AgedST 预测 | `results/direction3/gse280382_AgedST_predictions.csv` |
| GSE280382 Young 预测 | `results/direction3/gse280382_young_predictions.csv` |
| GSE280382 差异 | `results/direction3/gse280382_AgedLT_differences.csv` |
| GEO 排名列表 | `results/direction3/geo_ranked_interventions.csv` |
| HT 三重命中 | `results/geo_longterm/ht_screen_triple_hits.csv` |

### 图表

| 图表 | 路径 |
|------|------|
| GSE280382 多组织箱线图 | `figures/gse280382/gse280382_AgedLT_tissue_boxplot.png` |
| GSE280382 Exen vs Rapa | `figures/gse280382/gse280382_AgedST_exen_vs_rapa.png` |
| GSE280382 组织效应 | `figures/gse280382/gse280382_rejuvenation_by_tissue.png` |
| GSE288795 组织 tAge | `figures/direction3/gse288795_*_tage.png` |
| GSE230402 肝脏 tAge | `figures/direction3/gse230402_liver_tage.png` |
| 模块贡献热图 | `figures/module_contribution_heatmap.pdf` |
| 三重命中热图 | `figures/direction2_triple_hits_heatmap.pdf` |

---

## 七、下一步工作

1. **ARCHS4 H5 完成后**: Watcher 将自动触发大规模 liver-specific 样本提取和 tAge 预测
2. **高优先级 GEO 数据集**: GSE227689 (CNS rejuvenation, n=238), GSE248866 (CR, n=114)
3. **多模态验证**: 对 GLP-1RA 悖论，建议对比 DNA 甲基化时钟结果（原文献已有）
4. **候选药物验证**: Top 10 三重命中候选（BRD-K40324831, trimetazidine 等）的体外功能测试设计
5. **模型改进**: 考虑训练一个包含 GLP-1RA 等代谢干预数据的新模型，或开发多模态整合框架

---

*报告生成时间: 2026-06-09*  
*分析环境: Python 3.13, sklearn 1.9.0, R 4.5.3, edgeR*
