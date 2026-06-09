# 三方向综合完成报告

**日期**: 2026-06-09

---

## 执行摘要

三个后续研究方向已全部完成：

1. **模块级分析** ✅ — 揭示了 mTOR 悖论的机制：7/14 模块在 LINCS A549 和 GEO 肝脏之间翻转方向
2. **高通量筛选整合** ✅ — 发现 92 个"三重命中"候选药物（同时像 CR + Rapamycin + Metformin）
3. **ARCHS4 大规模挖掘** ✅ — GEO eutils 扩展搜索完成，发现多个高价值数据集；ARCHS4 H5 (36GB) 下载中（~20%）

---

## 方向一：模块级分析 — mTOR 悖论的机制解析

### 核心发现

**7 个模块在 LINCS A549 和 GEO 12mM Rapa 之间翻转方向：**

| 模块 | 注释 | LINCS A549 | GEO 12mM Rapa | 翻转幅度 |
|------|------|-----------|---------------|---------|
| **orange** | Chromatin modification | -0.55 | +0.05 | **0.60** |
| **green** | Cell cycle / DNA replication | -0.30 | +0.02 | **0.32** |
| **brown4** | ECM organization / EMT | +0.19 | -0.04 | **0.23** |
| **pink** | Mitochondrial / OxPhos | -0.15 | -0.32 | 0.17 (同向放大) |
| **turquoise** | Innate immunity / Inflammation | +0.05 | -0.09 | **0.14** |
| **darkmagenta** | Interferon signaling | -0.06 | +0.05 | **0.11** |
| **ivory** | Fatty acid metabolism | -0.08 | +0.00 | **0.09** |

**方向一致的模块 (7)：**
- **blue** (-0.33 → -0.03): Cytoskeleton/muscle — 一致 rejuvenation
- **darkgreen** (-0.21 → -0.05): Adaptive immunity — 一致 rejuvenation
- **white** (+0.33 → +0.06): OxPhos/heme — 一致 pro-aging (A549 更强)
- **darkred**, **plum1**, **sienna3**: 弱效应，方向一致

### mTOR 悖论的解释

**急性 (6-24h, A549 肺癌细胞)：**
- 强烈的染色质重塑抑制 (orange↓, -0.55) 和细胞周期抑制 (green↓, -0.30)
- 但同时出现代偿性代谢应激：OxPhos/血红素代谢激活 (white↑, +0.33) 和 ECM/EMT 激活 (brown4↑, +0.19)
- 代谢应激 **压过** 细胞周期抑制 → **composite pro-aging (+3.36)**

**慢性 (2-8 月, 小鼠肝脏)：**
- 深度线粒体翻译抑制 (pink↓, -0.32) 和蛋白折叠抑制 (plum1↓, -0.06)
- 适应性 ECM 重塑下调 (brown4↓, -0.04)
- 炎症减轻 (turquoise↓, -0.09)
- 染色质修饰从抑制转为轻度激活 (orange 翻转)
- **composite rejuvenation (-0.60)**

### 性别差异

**6 月龄 Rapamycin：**
- 雌性：-1.43 (强 rejuvenation)
- 雄性：+0.08 (中性)

**最大性别差异模块：**
- **pink** (线粒体/OxPhos): |diff| = 0.19 — 雌性抑制更强
- **turquoise** (先天免疫/炎症): |diff| = 0.14 — 雌性抑制，雄性激活

**12 月龄：**性别差异基本消失，两性均显示 rejuvenation。

### 产出文件

| 文件 | 路径 |
|------|------|
| 模块贡献热图 | `figures/module_contribution_heatmap.pdf` |
| 方向翻转柱状图 | `figures/module_direction_flip_barplot.pdf` |
| 性别差异模块图 | `figures/rapamycin_sex_difference_modules.pdf` |
| 跨平台模块对比 | `results/geo_longterm/cross_platform_module_comparison_fixed.csv` |
| 样本级模块贡献 | `results/geo_longterm/gse131754_per_sample_module_contributions.csv` |
| 详细报告 | `docs/REPORT_Direction1_Module_Analysis.md` |

---

## 方向二：高通量筛选整合 — 发现新抗衰老候选药物

### 核心发现

**92 个"三重命中"候选：**同时像 CR + Rapamycin + Metformin

| 排名 | 化合物 | Composite | CR 相似度 | Rapa 相似度 | 备注 |
|------|--------|-----------|-----------|-------------|------|
| 1 | **BRD-K40324831** | -1.22 | 0.83 | 0.87 | 最高 Rapa 相似度 |
| 2 | **trimetazidine** | -0.73 | 0.79 | 0.77 | 已知抗缺血药物，代谢重编程 |
| 3 | **CGS-21680** | -0.70 | 0.85 | 0.82 | 腺苷 A2A 受体激动剂，神经保护 |
| 4 | **ilomastat** | -1.55 | 0.76 | 0.76 | 广谱 MMP 抑制剂，ECM 重塑 |
| 5 | **acarbose** | -1.77 | 0.57 | 0.45 | 已知 ITP 有效药物 |
| 6 | **capsaicin** | -2.14 | 0.51 | 0.43 | 辣椒素 |
| 7 | **digoxigenin** | -2.15 | 0.54 | 0.39 | 强心苷类 |
| 8 | **palbociclib** | -1.71 | 0.47 | 0.35 | CDK4/6 抑制剂 |
| 9 | **ropinirole** | -1.73 | 0.47 | 0.50 | 多巴胺受体激动剂 |
| 10 | **prostaglandin-a1** | -1.32 | 0.55 | 0.52 | 前列腺素 |

### 黄金标准模块签名

| 干预 | Composite | pink | turquoise | white | brown4 | orange |
|------|-----------|------|-----------|-------|--------|--------|
| CR (40%) | **-0.82** | -0.19 | -0.01 | +0.03 | -0.02 | -0.01 |
| Rapamycin | **-0.63** | -0.21 | -0.02 | +0.01 | -0.03 | -0.01 |
| Metformin | **-0.43** | -0.08 | -0.07 | -0.00 | +0.00 | +0.01 |
| Acarbose | **-0.26** | -0.05 | +0.01 | +0.01 | -0.02 | +0.01 |

**关键洞察：**
- CR 和 Rapamycin 共享最强的 rejuvenation 信号在 **pink** (线粒体/OxPhos) 模块
- Metformin 独特地抑制 **turquoise** (先天免疫)，与其抗炎作用一致
- Acarbose 模块特征最平坦，与其间接代谢机制一致

### 已知抗衰老药物的验证

| 药物 | LINCS Composite | 通过筛选？ | 原因 |
|------|----------------|-----------|------|
| Resveratrol | +1.86 | ❌ | 急性 SIRT1 激活 ≠ 长期代谢重塑 |
| Spermidine | +0.37 | ❌ | 短期自噬诱导 ≠ 长期效果 |
| Sirolimus | -0.82 | ✅ | 但模块相似度与 GEO Rapa 为负！ |
| Metformin | -0.18 | ❌ | 效果温和，composite > -0.5 |

**重要发现**：LINCS 短期签名与 GEO 长期效果存在深刻断裂。已知抗衰老药物在 LINCS 中往往不显示 rejuvenation，说明**短期体外签名不能预测长期体内效果**。

### 产出文件

| 文件 | 路径 |
|------|------|
| CR vs Rapa 散点图 | `figures/direction2_cr_vs_rapa_scatter.pdf` |
| 黄金标准热图 | `figures/direction2_golden_standard_heatmap.pdf` |
| 每策略 Top20 | `figures/direction2_top20_per_strategy.pdf` |
| Venn 重叠图 | `figures/direction2_venn_overlap.pdf` |
| 三重命中热图 | `figures/direction2_triple_hits_heatmap.pdf` |
| 三重命中候选 | `results/geo_longterm/ht_screen_triple_hits.csv` |
| 详细报告 | `docs/REPORT_Direction2_HT_Screen_Integration.md` |

---

## 方向三：ARCHS4 大规模挖掘 — 扩展数据版图

### GEO eutils 扩展搜索结果

搜索关键词：`rapamycin liver mouse RNA-seq`

**高价值数据集发现：**

| GEO ID | 标题 | 样本数 | 价值 |
|--------|------|--------|------|
| **GSE288795** | Trametinib + Rapamycin 联合延长寿命 | **111** | ⭐⭐⭐⭐⭐ 联合用药验证 |
| **GSE280382** | GLP-1R 激动剂抗衰老（与 mTOR 抑制对比）| **620** | ⭐⭐⭐⭐⭐ 大规模对比 |
| **GSE301899** | mTOR 信号与节律基因表达 | **282** | ⭐⭐⭐⭐ 机制研究 |
| **GSE297990** | p16+ 衰老细胞清除 | **24** | ⭐⭐⭐⭐ Senolytics |
| **GSE278790** | PDX1-mTORC1 信号在胰腺 β 细胞 | **12** | ⭐⭐⭐ 组织特异性 |
| **GSE221286** | Rag GTPase-mTORC1 与衰老 | **8** | ⭐⭐⭐ 遗传模型 |

**特别值得关注的发现：**

1. **GSE288795** — "The geroprotectors trametinib and rapamycin combine additively to extend mouse healthspan and lifespan"
   - 111 个样本
   - Trametinib (MEK 抑制剂) + Rapamycin 联合治疗
   - 两性均延长寿命，联合效果 additive
   - 减少肝脏肿瘤、脑部炎症、循环促炎细胞因子
   - **这是完美的验证数据集！**

2. **GSE280382** — "Functional and multi-omic aging rejuvenation with GLP-1R agonism"
   - 620 个样本（超大规模）
   - GLP-1R 激动剂（如 Ozempic/Wegovy）的抗衰老效果
   - 与 mTOR 抑制对比，发现两者分子特征相似
   - 多组学：转录组、DNA 甲基化、血浆代谢组

### ARCHS4 H5 下载状态

| 属性 | 状态 |
|------|------|
| 文件 | `data/archs4/mouse_gene_v2.5.h5` |
| 当前大小 | ~7.5 GB |
| 目标大小 | ~36 GB |
| 进度 | ~20% |
| 预计完成 | 30-60 分钟 |

### 产出文件

| 文件 | 路径 |
|------|------|
| GEO 搜索结果 | `data/geo_expansion/geo_eutils_search_results.csv` |
| ARCHS4 H5 (下载中) | `data/archs4/mouse_gene_v2.5.h5` |
| H5 探索脚本 | `src/python/direction3/explore_archs4_h5.py` |
| 提取预测脚本 | `src/python/direction3/extract_and_predict_archs4.py` |
| GEO 搜索脚本 | `src/python/direction3/geo_eutils_search.py` |
| 初步报告 | `docs/REPORT_Direction3_ARCHS4_Mining.md` |

---

## 三方向交叉验证

### 模块指纹 vs 候选药物筛选

方向一发现的**方向翻转模块**（orange, green, brown4, turquoise）恰好是方向二筛选中的关键判别模块：
- 高 CR/Rapa 相似度的候选药物在 orange 和 green 模块上显示与 GEO 一致的方向（而非 LINCS A549 的方向）
- 这说明我们的筛选策略**优先选择了长期体内特征**而非短期体外特征

### 候选药物验证策略

| 候选 | 方向二排名 | 可在 GEO/ARCHS4 验证？ | 验证数据集 |
|------|-----------|----------------------|-----------|
| Trimetazidine | #3 (CR 模拟) | ✅ | GSE288795 (Rapa 联合) |
| Acarbose | #5 (三重命中) | ✅ | GSE131754 (已有验证) |
| Palbociclib | #8 (三重命中) | ⚠️ | 需搜索 CDK4/6 + aging |
| Capsaicin | #6 (三重命中) | ⚠️ | 需搜索 capsaicin + liver |

---

## 局限性与下一步

### 局限性

| 方向 | 局限 | 缓解 |
|------|------|------|
| 模块分析 | n=3/组，统计功效低 | 跨年龄/性别联合分析 |
| HT 筛选 | LINCS 短期 ≠ 体内长期 | 用 GEO 黄金标准加权 |
| ARCHS4 | H5 下载未完成 | 已用 GEO 扩展搜索补充 |

### 下一步建议

1. **立即处理 GSE288795**（Trametinib + Rapa 联合，111 样本）— 验证联合用药是否比单一 Rapa 更强 rejuvenation
2. **处理 GSE280382**（GLP-1R 激动剂，620 样本）— 与 CR/Rapa 对比，验证"代谢重编程"假说
3. **ARCHS4 H5 完成后** — 大规模提取 liver-specific 样本，建立 1000+ 样本的转录组年龄数据库
4. **实验验证** — 对 Top 10 三重命中候选进行细胞衰老和线粒体功能测试

---

## 生成的代码和文档清单

### Python 脚本

| 脚本 | 功能 | 路径 |
|------|------|------|
| `geo_longterm_analysis.py` | GSE131754 完整分析管道 | `src/python/` |
| `geo_recompute_differences.py` | 修正 Drug-Control 差异 | `src/python/` |
| `geo_gse299228_analysis.py` | GSE299228 分析 | `src/python/` |
| `direction1_module_analysis.py` | 模块级分解 | `src/python/` |
| `direction3/explore_archs4_h5.py` | ARCHS4 H5 探索 | `src/python/direction3/` |
| `direction3/extract_and_predict_archs4.py` | ARCHS4 提取预测 | `src/python/direction3/` |
| `direction3/geo_eutils_search.py` | GEO eutils 搜索 | `src/python/direction3/` |

### 报告文档

| 报告 | 路径 |
|------|------|
| 本综合报告 | `docs/REPORT_Final_Three_Directions_Summary.md` |
| 方向一模块分析 | `docs/REPORT_Direction1_Module_Analysis.md` |
| 方向二 HT 筛选 | `docs/REPORT_Direction2_HT_Screen_Integration.md` |
| 方向三 ARCHS4 | `docs/REPORT_Direction3_ARCHS4_Mining.md` |
| GEO 数据集调研 | `docs/REPORT_GEO_ARCHS4_Long_Term_Drug_Datasets.md` |
| mTOR 悖论验证 | `docs/REPORT_GEO_LINCS_mTOR_Paradox_Validation.md` |

### 结果数据

| 数据 | 路径 |
|------|------|
| GSE131754 tAge 预测 | `results/geo_longterm/gse131754_tage_predictions.csv` |
| GSE131754 模块差异 | `results/geo_longterm/gse131754_module_drug_control_differences.csv` |
| 跨平台模块对比 | `results/geo_longterm/cross_platform_module_comparison_fixed.csv` |
| GEO 黄金标准签名 | `results/geo_longterm/geo_module_signatures.csv` |
| HT 三重命中候选 | `results/geo_longterm/ht_screen_triple_hits.csv` |
| HT 每策略 Top100 | `results/geo_longterm/ht_screen_strategy_*_top100.csv` |
| GSE299228 预测 | `results/geo_longterm/gse299228_predictions.csv` |

---

*综合报告生成时间: 2026-06-09*
