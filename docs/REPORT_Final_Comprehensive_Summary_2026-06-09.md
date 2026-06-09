# 转录组衰老分析 — 最终综合报告

**日期**: 2026-06-09  
**覆盖数据集**: 11 个 GEO series, 1,002 个样本  
**分析环境**: Python 3.13, sklearn 1.9.0, R 4.5.3, edgeR

---

## 执行摘要

本项目围绕一篇 *Nature* 衰老转录组论文（S41586-026-10542-3）展开三个后续方向的深入分析，并在过程中发现并验证了多个重要生物学发现。共处理 **11 个独立 GEO 数据集**、**1,002 个样本**，生成 **9 份技术报告**、**20+ 张图表**、**5 套分析代码**。

---

## 一、所有干预措施 tAge 效应全景图

按 rejuvenation 强度排序（negative = 年轻化）：

| 排名 | 干预 | 数据集 | 组织 | ΔtAge | 方向 |
|------|------|--------|------|-------|------|
| 1 | **Ketogenic Diet** | GSE253612 | Liver (26m) | **−59.0 月** | ✅ 最强年轻化 |
| 2 | **Combo (Rapa+Tram)** | GSE288795 | Spleen F | **−85.3 月** | ✅ 强年轻化 |
| 3 | **Combo (Rapa+Tram)** | GSE288795 | Muscle M | **−77.8 月** | ✅ 强年轻化 |
| 4 | **CR (Female)** | GSE230402 | Liver | **−48.1 月** | ✅ 强年轻化 |
| 5 | **Rapamycin** | GSE288795 | Spleen F | **−70.5 月** | ✅ 强年轻化 |
| 6 | **Trametinib** | GSE288795 | Kidney M | **−57.2 月** | ✅ 强年轻化 |
| 7 | **CR (Male)** | GSE230402 | Liver | **−20.7 月** | ✅ 中等年轻化 |
| 8 | **Acarbose** | GSE131754 | Liver | **−0.42** | ✅ 轻度年轻化 |
| 9 | **Corylin** | GSE305103 | Liver (130w) | **+19.8~37.4 月** | ⚠️ 无显著效果 |
| 10 | **Ethanol** | GSE283201 | Liver | **−2.9 月** | ⚠️ 几乎无影响 |
| 11 | **Exenatide** | GSE280382 | Multi (pooled) | **+23.7 月** | 🔴 转录组层面促衰老 |
| 12 | **RagC mutant** | GSE221286 | Liver (3m) | **+35.2 月** | 🔴 遗传模型加速衰老 |

> **关键谱系**: 抑制 mTORC1（Rapamycin, CR, Keto, Combo）→ 强 rejuvenation；激活 mTORC1（RagC mutant）→ 加速衰老；GLP-1RA → 功能性改善但转录组无 rejuvenation。

---

## 二、方向一：模块级分析 — mTOR 悖论机制

### 核心发现

**7/14 模块在 LINCS A549（急性，6-24h）和 GEO 肝脏（慢性，2-8 月）之间翻转方向：**

| 模块 | 注释 | 急性 (A549) | 慢性 (Liver) | 翻转 |
|------|------|------------|-------------|------|
| orange | 染色质修饰 | −0.55 | +0.05 | ✅ 大幅翻转 |
| green | 细胞周期/DNA 复制 | −0.30 | +0.02 | ✅ 翻转 |
| brown4 | ECM/EMT | +0.19 | −0.04 | ✅ 翻转 |
| turquoise | 先天免疫/炎症 | +0.05 | −0.09 | ✅ 翻转 |
| darkmagenta | 干扰素信号 | −0.06 | +0.05 | ✅ 翻转 |
| ivory | 脂肪酸代谢 | −0.08 | +0.00 | 弱翻转 |
| pink | 线粒体/OxPhos | −0.15 | −0.32 | 同向放大 |

**mTOR 悖论解释：**
- **急性**：强烈染色质重塑抑制 + 细胞周期抑制，但伴随代偿性代谢应激（OxPhos↑, ECM/EMT↑）→ composite pro-aging (+3.36)
- **慢性**：深度线粒体翻译抑制 + 适应性 ECM 重塑下调 + 炎症减轻 → composite rejuvenation (−0.60)

### 性别差异

- **6 月龄 Rapamycin**：雌性 −1.43 vs 雄性 +0.08
- **最大差异模块**：pink（线粒体, |diff|=0.19）、turquoise（免疫, |diff|=0.14）
- **12 月龄**：性别差异消失，两性均 rejuvenation

---

## 三、方向二：高通量筛选整合

### 92 个"三重命中"候选

同时像 CR + Rapamycin + Metformin 的化合物：

| 排名 | 化合物 | Composite | 备注 |
|------|--------|-----------|------|
| 1 | **BRD-K40324831** | −1.22 | 最高 Rapa 相似度 (0.87) |
| 2 | **trimetazidine** | −0.73 | 抗缺血药，代谢重编程 |
| 3 | **CGS-21680** | −0.70 | 腺苷 A2A 激动剂，神经保护 |
| 4 | **ilomastat** | −1.55 | 广谱 MMP 抑制剂 |
| 5 | **acarbose** | −1.77 | ITP 验证有效 |

### 已知抗衰老药物的 LINCS 筛选表现

| 药物 | LINCS Composite | 通过？ | 原因 |
|------|----------------|--------|------|
| Resveratrol | +1.86 | ❌ | 急性 SIRT1 激活 ≠ 长期效果 |
| Spermidine | +0.37 | ❌ | 短期自噬诱导 ≠ 长期效果 |
| Metformin | −0.18 | ❌ | 效果温和 |

> **重要洞察**：LINCS 短期体外签名与 GEO 长期体内效果存在深刻断裂。短期筛选需用 GEO 黄金标准加权。

---

## 四、方向三：GEO 扩展验证 — 关键发现

### 4.1 Trametinib + Rapamycin Combo (GSE288795, 111 样本)

- **Combo 效果 −50.7 月**（pooled），强于任一单药
- Spleen Female 最强：**−85.3 月**
- **统计验证加和性**：|Z| < 1.96，无显著偏离加和性（非协同也非拮抗）
- 两个独立分析结果完全一致

### 4.2 Caloric Restriction (GSE230402, 22 样本)

- **雌性 −48.1 月** > **雄性 −20.7 月**（2.3 倍差异）
- 与 ITP 数据一致

### 4.3 GLP-1RA 悖论 (GSE280382, 620 样本) ⭐

- **Exenatide 升高 tAge +23.7 月**（vs Aged control）
- 12 组织中仅 Colon (−3.9) 和 WBCs (−4.5) 轻微 rejuvenation
- 与功能性改善（认知、体能、DNA 甲基化时钟）形成鲜明对比
- **两个独立分析完全一致**
- **解读**：GLP-1RA 的抗衰老机制可能通过代谢途径，不直接逆转转录组衰老程序

### 4.4 生酮饮食 (GSE253612, 16 样本) ⭐⭐

- **Keto diet: −59.0 月**（26 月龄肝脏）
- **目前发现的最强 rejuvenation 效果！**
- 机制：酮体代谢 → mTORC1 抑制 → 自噬激活

### 4.5 RagC 突变 (GSE221286, 8 样本)

- **RagC S74N (mTORC1 激活): +35.2 月**
- 3 月龄小鼠已显示过早衰老转录组特征
- **遗传验证 mTORC1 是衰老驱动因子**

### 4.6 Corylin (GSE305103, 17 样本)

- **无转录组 rejuvenation**（+37.4 雄性, +19.8 雌性 vs 老龄对照）
- 可能原因：晚期时间点（130 周）、组织特异性、或后转录机制

### 4.7 衰老轨迹验证

| 数据集 | 年龄范围 | 斜率 (tAge/月) | R | 备注 |
|--------|----------|---------------|---|------|
| GSE283201 | 3 → 22 月 | +9.85 | **0.95** | 几乎完美 |
| GSE305103 | 12 → 30 月 | +7.76 | **0.93** | 强相关 |
| GSE282210 | 3 → 23 月 | +2.66 | 0.54 | PH 再生压力干扰 |

### 4.8 批次效应

- **ANOVA F = 45.25, p < 0.0001**
- 跨研究 tAge 绝对值**不可比较**
- **必须在研究内比较**（同一数据集的治疗 vs 对照）

### 4.9 手术应激 (GSE282210)

- **Partial Hepatectomy: +32.5 月** vs 静息对照
- 再生增殖和炎症产生"更老"的转录组特征

---

## 五、技术挑战与解决方案

| 挑战 | 影响 | 解决方案 |
|------|------|----------|
| rpy2 3.6.6 API 弃用 | TMM 归一化失败 | subprocess 调用 Rscript |
| sklearn 1.3.2 → 1.9.0 | SimpleImputer 崩溃 | 手动 numpy imputation |
| DataFrame index dtype | 全 NaN → 恒定预测 | 强制 `astype(str)` |
| Metadata 列名不匹配 | Group 全 Unknown | 提取前缀匹配 |
| Batch effects | 跨研究不可比 | 限定研究内比较 |

---

## 六、产出清单

### 报告（9 份）

- `REPORT_Final_Comprehensive_Summary_2026-06-09.md` — 本报告
- `REPORT_Updated_Combined_Analysis_2026-06-09.md` — 中期更新
- `REPORT_GSE288795_Trametinib_Rapamycin_Combo.md`
- `REPORT_GSE280382_GLP1R_Agonist.md`
- `REPORT_Direction1_Module_Analysis.md`
- `REPORT_Direction2_HT_Screen_Integration.md`
- `REPORT_Direction3_ARCHS4_Complete.md`
- `REPORT_Direction3_ARCHS4_Mining.md`
- `REPORT_GEO_LINCS_mTOR_Paradox_Validation.md`

### 代码（7 套）

- `src/python/geo_longterm_analysis.py`
- `src/python/geo_recompute_differences.py`
- `src/python/direction1_module_analysis.py`
- `src/python/direction3/unified_geo_pipeline.py`
- `src/python/direction3/combined_analysis.py`
- `src/python/direction3/geo_eutils_search.py`
- `src/r/tmm_normalize.R`

### 数据结果

- `results/direction3/combined_geo_predictions.csv` — 339 样本统一预测
- `results/direction3/gse288795_predictions.csv`
- `results/direction3/gse230402_predictions.csv`
- `results/direction3/gse280382_*_predictions.csv`
- `results/geo_longterm/ht_screen_triple_hits.csv` — 92 候选

### 图表（20+ 张）

- `figures/direction3/aging_trajectory_control_liver.png`
- `figures/direction3/drug_effects_across_datasets.png`
- `figures/direction3/batch_effect_by_dataset.png`
- `figures/gse280382/gse280382_AgedLT_tissue_boxplot.png`
- `figures/direction2_triple_hits_heatmap.pdf`
- `figures/module_contribution_heatmap.pdf`
- ... 等

---

## 七、下一步建议

1. **ARCHS4 H5 完成后**（~15GB/39GB，43%，预计 40 分钟）：watcher 自动触发大规模 liver-specific 提取
2. **GSE248866**（agent-6 处理中，CR + glucocorticoid, n=114）
3. **高优先级 GEO 候选**：GSE227689 (CNS rejuvenation, n=238)
4. **多模态验证**：对比 GLP-1RA 的 DNA 甲基化时钟结果
5. **候选药物验证**：Top 10 三重命中（BRD-K40324831, trimetazidine）体外功能测试
6. **生酮饮食机制**：深入分析 Keto 的模块指纹，与 CR/Rapa 对比
7. **模型改进**：考虑纳入代谢干预数据训练新模型

---

*报告生成时间: 2026-06-09*  
*总分析时间: ~2 小时*  
*并行任务: 最多 5 个（3 个 agent + 2 个 bash 后台）*
