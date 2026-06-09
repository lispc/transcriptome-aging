# 方向 E：GEO/ARCHS4 长期给药数据挖掘 — 数据集调研报告

**日期**: 2026-06-08

---

## 摘要

经过对 GEO (Gene Expression Omnibus) 和 ARCHS4 数据库的系统搜索，我们发现了一个**近乎理想的数据集**——GSE131754（Tyshkovskiy et al., *Cell Metabolism* 2019），它包含了 ITP（Intervention Testing Program）中几乎所有关键抗衰老干预的小鼠肝脏转录组数据。此外，我们还鉴定了多个补充数据集，覆盖 Rapamycin、Metformin、Caloric Restriction 等干预。本报告详细评估了每个数据集的数据类型、实验设计、样本量、与我们的 tAge 预测管道的兼容性，以及推荐使用策略。

---

## 1. 核心发现：GSE131754 — "近乎完美"的 ITP 干预数据集

### 1.1 基本信息

| 属性 | 详情 |
|------|------|
| **GEO ID** | GSE131754 |
| **标题** | RNA sequencing of mouse hepatic response to lifespan-extending interventions |
| **发表期刊** | *Cell Metabolism* (2019) |
| **PubMed** | [PMID: 31353263](https://pubmed.ncbi.nlm.nih.gov/31353263/) |
| **样本数** | **78** (3 生物重复/组) |
| **数据类型** | **RNA-seq**, 配对端 100bp |
| **组织** | 肝脏 (liver) |
| **品系** | **UM-HET3** (ITP 标准杂交品系!) |

### 1.2 实验设计 — 8 种干预 + 对照

| 干预 | 处理组 | 对照组 | 年龄 | 处理时长 | 剂量 |
|------|--------|--------|------|----------|------|
| **Acarbose** | 3 (M+F) | 3 (M+F) | 6 月 | 2 月 | 1000 ppm |
| **Acarbose** | 3 (M+F) | 3 (M+F) | 12 月 | 8 月 | 1000 ppm |
| **Caloric Restriction (CR)** | 3 (M+F) | 3 (M+F) | 6 月 | 2 月 | 40% 限制 |
| **Caloric Restriction (CR)** | 3 (M+F) | 3 (M+F) | 12 月 | 8 月 | 40% 限制 |
| **Rapamycin** | 3 (M+F) | 3 (M+F) | 6 月 | 2 月 | 42 ppm |
| **Rapamycin** | 3 (M+F) | 3 (M+F) | 12 月 | 8 月 | 14 ppm |
| **17α-estradiol (17aE2)** | 3 (M+F) | 3 (M+F) | 6 月 | 2 月 | 14.4 ppm |
| **Protandim** | 3 (M+F) | 3 (M+F) | 6 月 | 2 月 | 1200 ppm |
| **Methionine Restriction (MR)** | 3 (M) | 3 (M) | 2 月 | 12 月 | 0.12% w/w |
| **GHRKO** | 3 (M) | 3 (M) | 5 月 | — | 基因敲除 |
| **Snell dwarf** | 3 (M) | 3 (M) | 5 月 | — | Pit1 -/- |

**总计**: 78 个样本，每组 3 个生物重复。

### 1.3 为什么这是"近乎完美"的数据集？

1. **品系匹配 ITP**: UM-HET3 是 ITP 使用的标准四向杂交品系 (C57BL/6J × BALB/cByJ × C3H/HeJ × DBA/2J)，我们的 Phase 0 rodent ITP 验证也使用此品系数据。
2. **干预全面**: 覆盖 ITP 中所有经证实延寿的小分子药物（Acarbose、Rapamycin、17aE2）和 CR。
3. **年龄分层**: 6 月龄（短期 2 月）和 12 月龄（长期 8 月），可以比较短期 vs 长期效应。
4. **性别分层**: 雄性和雌性都有，可以检测性别差异（该论文的核心发现之一）。
5. ** littermate controls**: 每组都有年龄、性别、品系匹配的同窝对照。
6. **RNA-seq 数据类型**: 与我们的 tAge 预测管道完全兼容（TMM normalization + logCPM）。
7. **发表的高质量**: *Cell Metabolism* 期刊，已被广泛引用，数据可靠性高。

### 1.4 已知发现（来自原文）

- 许多干预措施表现出相似的转录组变化
- **Rapamycin 显示出独特的模式**（与其他干预不同）
- 检测到与生长激素调控相关的**女性化效应**
- 氧化磷酸化和 NRF2 调控酶的上调是寿命延长的共同特征
- 葡萄糖代谢和免疫功能相关基因的变化与长寿有定量和定性关联

### 1.5 数据可用性

- GEO 页面: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131754
- SRA 原始数据: 可通过 GSM 样本 ID 下载 FASTQ
- 可能的 processed counts: 需要检查 GEO supplementary files

---

## 2. 补充数据集

### 2.1 GSE299228 — 多干预比较（2026 年新发表）

| 属性 | 详情 |
|------|------|
| **GEO ID** | GSE299228 |
| **标题** | RNA Sequencing of Mouse Liver Post Caloric Resrtiction and Anti-Anging Interventions |
| **样本数** | **18** |
| **数据类型** | **RNA-seq** |
| **组织** | 肝脏 |
| **品系** | C57BL/6J |
| **年龄** | 21 周龄开始 |
| **处理时长** | 8 周 |
| **PubMed** | 未收录（2026 年 3 月公开） |

**实验设计**（5 组 × 3-6 重复）：

| 组 | n | 干预 |
|----|---|------|
| Control (ad libitum) | 6 | 对照 |
| Caloric Restriction | 3 | 30% CR |
| TM5614 | 3 | PAI-1 抑制剂 20 mg/kg/day |
| **Rapamycin** | 3 | 42 mg/kg/day |
| **Metformin** | 3 | 37.5 mg/kg/day |

**亮点**:
- **已提供 DESeq2 归一化 counts** (`GSE299228_CR_normalized_counts_DESeq2.csv.gz`)
- 可以直接下载使用，无需重新比对
- 比较了 CR 的药物模拟物
- TM5614 被鉴定为最接近 CR 的模拟物（共享 ~50% DEGs）
- **Rapamycin 和 Metformin 与 CR 的重叠最小**

**与 GSE131754 的互补性**:
- GSE299228 有 Metformin（GSE131754 没有）
- GSE131754 有更多干预和年龄组
- 两者都有 Rapamycin 和 CR，可以交叉验证

### 2.2 GSE48333 — Rapamycin 长期处理（微阵列）

| 属性 | 详情 |
|------|------|
| **GEO ID** | GSE48333 |
| **标题** | The effect of chronic Rapamycin on the transcriptome of old mouse liver |
| **样本数** | **69** (19 Control M, 13 Rapa M, 16 Control F, 21 Rapa F) |
| **数据类型** | **微阵列** (Affymetrix Mouse Gene 1.0 ST Array, GPL6885) |
| **组织** | 肝脏 |
| **品系** | C57BL6/J |
| **年龄** | 25 月龄（从 4 月龄开始处理） |
| **处理时长** | **21 个月** |
| **PubMed** | [PMID: 24409289](https://pubmed.ncbi.nlm.nih.gov/24409289/) |

**亮点**:
- **真正的长期处理**: 21 个月，从 4 月龄到 25 月龄
- **大样本量**: 69 个样本
- **性别分层**: 雄性和雌性
- Wilson et al. (2014)，较早期但重要

**⚠️ 数据类型兼容性问题**:

这是**微阵列**数据，不是 RNA-seq。我们的 tAge 预测管道基于 RNA-seq 数据训练（TMM normalization + edgeR + logCPM）。微阵列数据需要不同的预处理流程：

| 步骤 | RNA-seq (我们的管道) | 微阵列 |
|------|---------------------|--------|
| 输入 | Raw counts | Fluorescence intensities |
| 背景校正 | 不需要 | RMA / MAS5 |
| 归一化 | TMM (edgeR) | Quantile normalization |
| 转换 | logCPM | log2 |
| 基因 ID | Entrez Gene | Probe ID → Gene Symbol 映射 |
| 动态范围 | 宽 (0 ~ 10^6 counts) | 窄 (log2 ~ 5-15) |

**结论**: GSE48333 的数据**不能直接用我们的 tAge 管道**，需要：
1. 使用 `oligo` 或 `affy` 包进行 RMA 背景校正和归一化
2. 探针到基因 Symbol 的映射（可能需要 `mogene10sttranscriptcluster.db`）
3. 验证微阵列 log2 表达值与 RNA-seq logCPM 的分布匹配
4. 或者，使用交叉平台标准化方法（如 ComBat、quantile normalization）

**替代方案**: 如果微阵列数据太麻烦，可以跳过 GSE48333，因为 GSE131754 已经有 Rapamycin 的 RNA-seq 数据。

### 2.3 GSE319230 — Metformin 在 SCA8 模型

| 属性 | 详情 |
|------|------|
| **GEO ID** | GSE319230 |
| **标题** | Metformin reduces RAN proteins and rescues molecular and behavioral phenotypes in SCA8 mice |
| **样本数** | 12 |
| **数据类型** | **RNA-seq** |
| **组织** | 小脑 |
| **PubMed** | [PMID: 41771688](https://pubmed.ncbi.nlm.nih.gov/41771688/) |

**⚠️ 问题**: SCA8（脊髓小脑共济失调 8 型）是疾病模型，存在显著的病理背景噪音，可能掩盖 Metformin 的抗衰老效应。

### 2.4 GSE293164 — Metformin 肠道机制

| 属性 | 详情 |
|------|------|
| **GEO ID** | GSE293164 |
| **标题** | Metformin inhibits mitochondrial complex I in intestinal epithelium to promote glycemic control |
| **样本数** | 10 |
| **数据类型** | **RNA-seq** |
| **组织** | 肠道上皮 |
| **PubMed** | [PMID: 41256385](https://pubmed.ncbi.nlm.nih.gov/41256385/) |

**问题**: 组织是肠道而非肝脏，且聚焦于血糖控制机制，与系统性衰老关联较弱。

### 2.5 GSE284988 — CR 对造血干细胞的影响

| 属性 | 详情 |
|------|------|
| **GEO ID** | GSE284988 |
| **标题** | Epigenetic profiling of hematopoietic stem cells identifies KDR and PU.1 as regulators of aging transcriptome and caloric restriction response |
| **样本数** | 58（混合 ChIP-seq + RNA-seq）|
| **数据类型** | ChIP-seq + RNA-seq |
| **组织** | 造血干细胞 |
| **PubMed** | [PMID: 41720793](https://pubmed.ncbi.nlm.nih.gov/41720793/) |

**问题**: 组织特殊（HSCs），RNA-seq 样本数可能较少（需要从 58 中分离）。

---

## 3. 缺失的药物数据

| 药物 | GEO 数据集数 | 备注 |
|------|-------------|------|
| **Acarbose** | ✅ GSE131754 | 有，但仅肝脏 |
| **Rapamycin** | ✅ GSE131754, GSE299228, GSE48333 | 充足 |
| **Metformin** | ⚠️ GSE299228, GSE319230, GSE293164 | 有，但组织/模型不理想 |
| **17α-estradiol** | ✅ GSE131754 | 有，但仅肝脏 |
| **Protandim** | ✅ GSE131754 | 有，但仅肝脏 |
| **Canagliflozin** | ❌ 0 | **无 GEO 数据集** |
| **NR (Nicotinamide Riboside)** | ❌ 0 | 无小鼠肝脏 RNA-seq |
| **NDGA** | ❌ 0 | 无 |
| **Astaxanthin** | ❌ 0 | 无 |
| **Metolazone** | ❌ 0 | 无 |

**Canagliflozin** 在 GEO 上完全没有小鼠 RNA-seq 数据。这是 ITP 中最新发现有效的药物之一（2023 年报告延长雄性寿命），但缺乏公开的转录组数据。

---

## 4. ARCHS4 资源评估

### 4.1 什么是 ARCHS4？

ARCHS4 (Massive Mining of Publicly Available RNA-seq Data) 是 Ma'ayan Lab 维护的资源，将 GEO/SRA 中的所有可处理 RNA-seq 数据统一用 Kallisto 比对到参考基因组，提供标准化的基因级别 pseudocounts。

### 4.2 数据规格

| 属性 | 详情 |
|------|------|
| **小鼠数据文件** | `mouse_gene_v2.5.h5` |
| **文件大小** | **36 GB** |
| **更新日期** | 2024-08-24 |
| **比对工具** | Kallisto (pseudoalignment) |
| **参考基因组** | GRCm38 (Ensembl 107) |
| **基因注释** | Entrez Gene Symbol |
| **计数类型** | Pseudocounts（取整为整数）|

### 4.3 与我们的管道兼容性

| 方面 | 评估 |
|------|------|
| **数据类型** | ✅ Kallisto pseudocounts 可以作为 raw counts 输入 edgeR/TMM |
| **基因 ID** | ✅ Entrez Gene Symbol，与我们的 .pkl 模型匹配 |
| **批次效应** | ⚠️ **严重问题**: 跨实验、跨平台、跨实验室的批次效应 |
| **组织混杂** | ⚠️ ARCHS4 包含所有组织的混合数据，需要仔细筛选 |

### 4.4 使用策略

ARCHS4 的价值在于**补充 GEO 的缺口**，特别是当我们需要：
1. 更多重复样本以增加统计功效
2. 不同组织的数据
3. 特定药物在其他模型中的数据

**但 ARCHS4 的主要问题**：

1. **文件太大**: 36GB H5 文件下载和存储成本高
2. **查询复杂**: 需要通过 GEO 元数据筛选样本，然后从 H5 提取
3. **批次效应**: 不同实验的数据混杂，需要仔细的批次校正

**可行的使用方式**:

```python
# 伪代码：从 ARCHS4 H5 提取特定样本
import h5py

with h5py.File('mouse_gene_v2.5.h5', 'r') as f:
    # 获取所有样本的 GSM ID
    samples = f['meta']['Sample_geo_accession'][:]
    # 找到目标 GSM
    target_idx = [i for i, s in enumerate(samples) if s in target_gsms]
    # 提取表达矩阵
    expr = f['data']['expression'][:, target_idx]
    genes = f['meta']['genes'][:]
```

### 4.5 ARCHS4 查询脚本

ARCHS4 网站提供自动生成的 R 脚本，可以根据元数据关键词搜索并提取样本。例如，搜索 "rapamycin" 会返回所有相关 GSM，然后生成 R 脚本提取这些样本的表达矩阵。

**但是**，这仍然需要下载 36GB 的 H5 文件才能运行脚本。

---

## 5. 数据类型兼容性详细评估

### 5.1 可直接使用的 RNA-seq 数据集

以下数据集**可以直接用我们的 tAge 预测管道**处理：

| 数据集 | 组织 | 需要的预处理 |
|--------|------|-------------|
| GSE131754 | 肝脏 | SRA → FASTQ → STAR/Kallisto → counts → TMM/logCPM |
| GSE299228 | 肝脏 | **已提供 DESeq2 counts**，可直接 TMM/logCPM |
| GSE319230 | 小脑 | SRA → FASTQ → counts → TMM/logCPM |
| GSE293164 | 肠道 | SRA → FASTQ → counts → TMM/logCPM |

### 5.2 需要适配的数据集

| 数据集 | 问题 | 解决方案 |
|--------|------|----------|
| GSE48333 | 微阵列 (Affymetrix) | RMA → probe→gene 映射 → quantile normalization → 验证分布 |

### 5.3 微阵列适配的技术细节

如果决定使用 GSE48333，需要以下步骤：

```r
# R 代码示例
library(oligo)
library(mogene10sttranscriptcluster.db)

# 读取 CEL 文件
celFiles <- list.celfiles("GSE48333_RAW/", full.names=TRUE)
data <- read.celfiles(celFiles)

# RMA 背景校正和归一化
eset <- rma(data)

# 探针到基因 Symbol 映射
probe_ids <- featureNames(eset)
gene_symbols <- mapIds(mogene10sttranscriptcluster.db, 
                       keys=probe_ids, 
                       column="SYMBOL", 
                       keytype="PROBEID")

# 聚合探针到基因（取中位数）
expr <- exprs(eset)
# ... 聚合和转换
```

**关键问题**: 微阵列的 log2 表达值与 RNA-seq 的 logCPM 分布不同。需要验证：
1. 微阵列数据的基因表达排名是否与 RNA-seq 一致
2. 微阵列的缺失探针（未覆盖的基因）在我们的 10,487 个模型基因中的比例
3. 是否可以用 cross-platform normalization（如 ComBat）校正

---

## 6. 推荐的使用策略

### 6.1 优先级排序

| 优先级 | 数据集 | 理由 |
|--------|--------|------|
| **P0 (立即)** | GSE131754 | 最全面、品系匹配、RNA-seq、多干预、多年龄 |
| **P1 (高)** | GSE299228 | 有 Metformin，已提供 counts，可与 GSE131754 交叉验证 |
| **P2 (中)** | GSE48333 | 仅用于 Rapamycin 长期效应验证，但需微阵列适配 |
| **P3 (低)** | GSE319230, GSE293164 | Metformin 补充，但组织/模型不理想 |

### 6.2 建议的下游分析流程

```
Step 1: 下载 GSE131754 的 SRA 原始数据
        ↓
Step 2: 使用 STAR 或 Kallisto 比对到 GRCm38
        ↓
Step 3: edgeR TMM normalization → logCPM
        ↓
Step 4: 提取 10,487 个模型基因（Entrez ID）
        ↓
Step 5: tAge 预测（Python 管道）
        ↓
Step 6: 药物 vs 对照的差值（Drug - Control）
        ↓
Step 7: 与 LINCS 短期签名比较
        ↓
Step 8: 模块指纹分解，验证 mTOR 悖论
```

### 6.3 具体研究问题

使用这些数据集，我们可以回答：

1. **LINCS 短期 vs 长期一致性**: 
   - Rapamycin 在 A549 中显示 pro-aging (+3.36)，但在小鼠肝脏长期处理中是 rejuvenating？
   - GSE131754 的 Rapamycin 签名是否与我们的 LINCS 模块指纹一致？

2. **药物特异性 vs 共同机制**:
   - 不同干预（CR, Rapamycin, Acarbose）是否有共同的抗衰老模块？
   - 原文说 Rapamycin 有独特模式，模块分解能否验证？

3. **性别差异**:
   - GSE131754 有 M/F 分层，可以验证抗衰老干预是否存在性别差异
   - 这与 LINCS 数据（无性别信息）形成互补

4. **年龄依赖性**:
   - 6 月 vs 12 月的处理效应是否不同？
   - 长期处理（GSE48333 的 21 个月）vs 短期（2-8 个月）

5. **新药物发现**:
   - 使用 GSE131754 的数据训练/验证，预测 GEO 中其他药物的效果

---

## 7. 数据获取计划

### 7.1 GSE131754 数据下载

```bash
# 使用 SRA Toolkit 下载
# 首先获取 SRA run list
esearch -db gds -query "GSE131754[Accession]" | efetch -format docsum | grep -o 'SRR[0-9]*' > srr_list.txt

# 并行下载
prefetch --option-file srr_list.txt
fasterq-dump --split-files SRR*
```

**预计数据量**: 78 样本 × ~2GB (PE100) = ~156GB 原始 FASTQ

**替代方案**: 直接从 GEO 下载 processed counts（如果有提供）

### 7.2 GSE299228 数据下载

```bash
# 已提供 DESeq2 归一化 counts，直接下载
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299228/suppl/GSE299228_CR_normalized_counts_DESeq2.csv.gz
```

### 7.3 ARCHS4（可选）

```bash
# 下载 mouse H5 文件（36GB）
wget https://s3.amazonaws.com/mssm-seq-matrix-updates/macOS_gene_v2.5.h5
# 或
wget https://maayanlab.cloud/archs4/download.html
```

**注意**: 36GB 下载可能需要较长时间，建议仅在必要时使用。

---

## 8. 风险评估

| 风险 | 可能性 | 影响 | 缓解措施 |
|------|--------|------|----------|
| GSE131754 数据质量不佳 | 低 | 高 | 原始论文已验证，可检查 QC 指标 |
| 样本量太小（n=3/组） | 中 | 中 | 使用所有组联合分析，或与 GSE299228 合并 |
| 批次效应（跨 GEO 数据集） | 高 | 中 | 使用 ComBat-seq 或 surrogate variable analysis |
| 微阵列适配失败 | 中 | 低 | 优先使用 RNA-seq 数据，微阵列作为补充 |
| ARCHS4 下载/存储问题 | 中 | 低 | 可跳过 ARCHS4，GEO 数据已足够 |
| 药物剂量/处理时间不匹配 | 中 | 中 | 关注定性模式（方向/模块）而非绝对值 |

---

## 9. 结论

### 9.1 总体评估

| 方面 | 评分 | 说明 |
|------|------|------|
| **数据可用性** | ⭐⭐⭐⭐⭐ | GSE131754 近乎完美，多干预覆盖 |
| **数据类型匹配** | ⭐⭐⭐⭐☆ | 主要是 RNA-seq，仅一个微阵列 |
| **样本量** | ⭐⭐⭐☆☆ | n=3/组偏小，但多组联合可弥补 |
| **实验设计** | ⭐⭐⭐⭐⭐ | 年龄、性别、品系匹配 |
| **与 LINCS 互补性** | ⭐⭐⭐⭐⭐ | 长期体内 vs 短期体外，完美互补 |

### 9.2 关键建议

1. **立即启动 GSE131754 的分析**。这是目前发现的最好数据集，几乎可以回答我们所有的研究问题。

2. **使用 GSE299228 作为验证和补充**，特别是其 Metformin 数据和已提供的 counts 可以加速分析。

3. **GSE48333 作为可选补充**。如果时间允许，可以尝试微阵列适配；否则跳过，GSE131754 的 Rapamycin 数据已足够。

4. **ARCHS4 暂缓**。36GB 的下载成本较高，且批次效应问题需要额外处理。仅在需要更多重复样本时考虑。

5. **Canagliflozin 等缺失药物**需要通过其他方式获取数据（如联系 ITP 合作者或文献中查找未上传的数据）。

---

## 附录 A: 数据集汇总表

| GEO ID | 药物/干预 | 组织 | 类型 | 样本数 | 品系 | 年龄 | 时长 | 匹配对照 | PMID |
|--------|----------|------|------|--------|------|------|------|----------|------|
| **GSE131754** | CR/MR/Rapa/Aca/17aE2/Pro/GHRKO/Snell | 肝脏 | RNA-seq | 78 | UM-HET3 | 6/12 月 | 2-12 月 | ✅ Littermate | 31353263 |
| **GSE299228** | CR/Rapa/Met/TM5614 | 肝脏 | RNA-seq | 18 | C57BL/6J | 21 周 | 8 周 | ✅ | — |
| **GSE48333** | Rapamycin | 肝脏 | 微阵列 | 69 | C57BL6/J | 25 月 | 21 月 | ✅ | 24409289 |
| GSE319230 | Metformin | 小脑 | RNA-seq | 12 | — | — | — | ❌ (SCA8 模型) | 41771688 |
| GSE293164 | Metformin | 肠道 | RNA-seq | 10 | — | — | — | ✅ | 41256385 |
| GSE284988 | CR | HSCs | RNA-seq+ChIP | 58 | — | — | — | ✅ | 41720793 |
| GSE315851 | CR | 胃黏膜 | RNA-seq | 14 | — | — | — | ✅ | — |
| GSE289993 | CR | 肝脏 | RNA-seq | 10 | — | — | — | ✅ | — |

---

## 附录 B: 与 tAge 管道的兼容性检查清单

- [x] GSE131754: RNA-seq → TMM/logCPM → 10,487 genes → tAge ✅
- [x] GSE299228: RNA-seq (counts provided) → TMM/logCPM → tAge ✅
- [ ] GSE48333: 微阵列 → RMA → probe mapping → 分布验证 → tAge ⚠️
- [x] GSE319230: RNA-seq → TMM/logCPM → tAge ✅ (但 SCA8 噪音)
- [x] GSE293164: RNA-seq → TMM/logCPM → tAge ✅ (但肠道组织)

---

*报告生成时间: 2026-06-08*
*分析工具: NCBI E-utilities (eutils), GEO query, ARCHS4 documentation*
