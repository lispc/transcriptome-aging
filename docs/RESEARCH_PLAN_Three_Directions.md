# 三个后续研究方向：详细方案

**日期**: 2026-06-09

---

## 方向一：ARCHS4 扩展 — 从 78 样本到数万样本

### 1.1 核心思路

ARCHS4 包含 ~36GB 的小鼠 RNA-seq 数据（`mouse_gene_v2.5.h5`），覆盖数万个样本。我们的目标是：

1. **不下载完整 36GB H5 文件**，而是通过 ARCHS4 API / 网页查询提取特定 GSM 样本
2. 或者，**选择性下载**（如仅下载 liver-specific 样本子集）
3. 对提取的样本运行 tAge 预测，建立大规模的"转录组年龄数据库"
4. 识别哪些 GEO 实验中的药物处理显示 rejuvenation/acceleration 模式

### 1.2 ARCHS4 技术规格

| 属性 | 详情 |
|------|------|
| 文件 | `mouse_gene_v2.5.h5` (36GB) |
| 数据 | Kallisto pseudocounts, 基因级别 (Entrez Symbol) |
| 参考 | GRCm38, Ensembl 107 |
| 更新 | 2024-08-24 |
| 样本数 | ~100,000+ 小鼠样本 |

### 1.3 可行策略（无需下载 36GB）

#### 策略 A：ARCHS4 网页查询 + 自动生成 R 脚本

ARCHS4 网站支持：
1. 按关键词搜索样本（如 "rapamycin", "caloric restriction", "liver"）
2. 自动生成 R 脚本，从 H5 文件中提取匹配的 GSM 样本
3. 运行脚本生成 TSV 表达矩阵

**问题**：仍需下载 36GB H5 文件才能运行脚本。

#### 策略 B：直接使用 GEO 原始数据（我们已经在做）

我们已经从 GEO 直接下载了 GSE131754 和 GSE299228 的 counts。对于更多数据，可以继续：
1. 用 eutils 搜索特定关键词的 GEO 数据集
2. 下载 processed counts 或 SRA FASTQ
3. 走相同的 TMM → logCPM → tAge 管道

**优势**：数据质量可控，实验设计清晰，无需处理 ARCHS4 的批次效应

#### 策略 C：选择性下载 ARCHS4 子集

ARCHS4 也提供按**组织/细胞类型**预提取的数据。例如：
- liver samples only (~2-3GB?)
- 或者通过 ARCHS4 API 批量查询特定 GSM 列表

**实施步骤**：

```python
# 伪代码：查询 ARCHS4 API
import requests

# 搜索包含 "rapamycin" 的样本
search_url = "https://maayanlab.cloud/archs4/search/search"
payload = {
    "search": "rapamycin",
    "species": "mouse",
    "type": "sample"
}
resp = requests.post(search_url, json=payload)
gsms = resp.json()["results"]

# 生成提取脚本
# ... 下载 H5 并提取
```

### 1.4 最有价值的数据挖掘方向

| 搜索词 | 预期样本数 | 价值 |
|--------|-----------|------|
| "rapamycin" + "liver" | ~50-100 | 验证剂量/时间响应 |
| "metformin" + "liver" | ~30-50 | 补充 GSE299228 |
| "caloric restriction" + "liver" | ~100+ | 剂量-响应曲线 |
| "aging" + "liver" | ~500+ | 建立年龄基线 |
| "dexamethasone" + "liver" | ~20-30 | LINCS rejuvenation 验证 |
| "doxorubicin" + "heart" | ~30-50 | 心脏毒性 vs 衰老加速 |

### 1.5 风险评估

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| 36GB 下载耗时 | 高 | 中 | 只下载一次，夜间进行 |
| 批次效应严重 | 高 | 高 | 按实验分组分析，不跨实验比较绝对值 |
| 元数据不完整 | 中 | 中 | 手动验证关键样本的实验设计 |
| 存储空间不足 | 中 | 高 | 当前 155GB 可用，36GB 可承受 |

### 1.6 预期成果

- 建立包含 1000+ 样本的"小鼠转录组年龄数据库"
- 识别除 ITP 药物外的其他潜在抗衰老干预
- 发现新的衰老加速因子（如化疗药物、环境毒素）

---

## 方向二：模块级分析 — Rapamycin 的方向翻转机制

### 2.1 核心问题

为什么 Rapamycin 在 LINCS A549 中 **pro-aging (+3.36)**，但在 GEO 肝脏中 **rejuvenation (-0.63)**？

模块分解可以回答：
1. 哪些模块贡献了方向翻转？
2. 是同一模块从 pro-aging 变为 rejuvenation，还是不同模块主导？
3. 哪些模块是"体外特异性"、哪些是"体内特异性"？

### 2.2 方法

我们已有 14 个 WGCNA 模块的基因列表。对于每个数据集，计算每个模块的局部贡献：

```python
module_contrib = sum(coef[module_genes] * x_centered[module_genes])
```

#### 步骤 1：LINCS Rapamycin 模块分解

从 `module_fingerprint_v2_per_drug.csv` 中提取 Sirolimus 的各模块得分：

| 模块 | Sirolimus A549 | Sirolimus A375 |
|------|---------------|----------------|
| composite | +0.85 | -0.82 |
| blue | -0.18 | +0.14 |
| brown4 | +0.10 | -0.10 |
| darkgreen | -0.16 | +0.13 |
| darkred | +0.01 | -0.04 |
| green | **-0.22** | **+0.15** |
| ivory | -0.12 | +0.12 |
| orange | **-0.39** | **+0.24** |
| pink | +0.01 | -0.04 |
| plum1 | -0.04 | +0.04 |
| sienna3 | -0.01 | +0.05 |
| turquoise | **-0.11** | **+0.11** |
| white | +0.23 | +0.03 |

#### 步骤 2：GEO Rapamycin 模块分解

对于 GSE131754 的 Rapamycin 样本：
1. 计算每个样本的各模块局部贡献
2. 与匹配对照比较，得到模块级别的 Drug-Control 差异

#### 步骤 3：跨平台对比

制作"模块贡献热图"：
- 行：14 个模块
- 列：LINCS A549 / LINCS A375 / GEO 肝脏 6mF / GEO 肝脏 6mM / GEO 肝脏 12mF / GEO 肝脏 12mM
- 值：模块级别的 Drug-Control 差异

### 2.3 预期发现

| 假设 | 验证方法 |
|------|---------|
| A549 的 pro-aging 由 "green" 和 "orange" 模块驱动 | 检查这两个模块在 A549 中是否显著正向贡献 |
| 肝脏 rejuvenation 由不同模块（如 "turquoise"）驱动 | 检查 turquoise 在肝脏中是否负向贡献 |
| 雌性响应强与 "免疫/炎症" 模块相关 | 检查免疫相关模块在 6mF vs 6mM 中的差异 |

### 2.4 代码实现

需要修改 `geo_longterm_analysis.py`，保存 per-gene contributions，然后按模块聚合：

```python
def compute_module_contributions(contributions, feature_names, module_genes):
    """Aggregate per-gene contributions into module scores."""
    module_scores = {}
    for module_name, genes in module_genes.items():
        mask = np.isin(feature_names, genes)
        module_scores[module_name] = contributions[:, mask].sum(axis=1)
    return pd.DataFrame(module_scores)
```

### 2.5 风险评估

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| 模块定义与数据不完全匹配 | 中 | 中 | 使用论文提供的模块基因列表 |
| 批次效应掩盖模块信号 | 中 | 高 | 只做 Drug-Control 差值，消除基线 |
| 样本量小导致模块估计不准 | 中 | 中 | 合并所有 Rapamycin 样本增加 n |

---

## 方向三：结合高通量筛选 — 寻找 CR/Rapamycin 的"化学模拟物"

### 3.1 核心思路

我们已有 13,072 个 LINCS 化合物的模块指纹。现在可以：

1. **定义"黄金标准" rejuvenation 签名**：
   - CR 的 GEO 模块指纹（最稳定、最强的小分子干预）
   - Rapamycin 的 GEO 模块指纹（体内验证有效）
   - 两者的共识模块模式

2. **计算每个 LINCS 化合物与黄金标准的相似度**：
   - 模块级别的 cosine similarity
   - 加权相似度（按模块重要性加权）

3. **筛选新候选药物**：
   - 排除已知抗衰老药物（避免重复发现）
   - 排除 cytotoxic 化合物
   - 优先选择模块模式接近 CR/Rapamycin 但 composite 得分也显示 rejuvenation 的化合物

### 3.2 黄金标准签名的构建

#### 选项 A：仅用 GEO 数据（推荐）

从 GSE131754 计算每个干预的模块级 Drug-Control 差异：

```python
golden_standard = {
    "CR": module_diff_cr,      # 最可靠
    "Rapamycin": module_diff_rapa,  # 第二可靠
    "Acarbose": module_diff_aca,    # 第三
}
```

#### 选项 B：GEO + LINCS 共识

对同时在 GEO 和 LINCS 中显示 rejuvenation 的干预（如 Metformin），构建共识签名。

### 3.3 筛选算法

```python
def screen_for_mimics(lincs_df, golden_standard, top_n=50):
    """
    Find compounds whose module fingerprint mimics golden standard.
    
    Parameters:
    -----------
    lincs_df : DataFrame
        Module fingerprint per compound (from HT screening)
    golden_standard : dict
        Module name -> desired change (Drug - Control)
    top_n : int
        Number of top candidates to return
    
    Returns:
    --------
    candidates : DataFrame
        Ranked list of compounds with similarity scores
    """
    
    # Compute cosine similarity for each compound
    similarities = []
    for idx, row in lincs_df.iterrows():
        compound_modules = row[module_names].values
        similarity = cosine_similarity(
            compound_modules.reshape(1, -1),
            golden_vector.reshape(1, -1)
        )[0, 0]
        
        # Also require composite < -0.5 (rejuvenation threshold)
        if row['composite'] < -0.5:
            similarities.append({
                'compound': row['drug'],
                'similarity': similarity,
                'composite': row['composite'],
                'cytotoxic': row['orange'] < -0.5 and row['green'] < -0.5
            })
    
    candidates = pd.DataFrame(similarities)
    candidates = candidates[~candidates['cytotoxic']]  # exclude cytotoxic
    return candidates.nlargest(top_n, 'similarity')
```

### 3.4 最有前景的筛选策略

| 策略 | 黄金标准 | 预期发现 |
|------|---------|---------|
| **CR 模拟物** | CR 的模块指纹 | 找到"口服 CR"候选 |
| **Rapamycin 模拟物** | Rapa 的模块指纹 | 找到非 mTOR 机制的 mTOR 效应模拟物 |
| **CR+Rapa 共识** | 两者共同模块变化 | 找到多重机制抗衰老候选 |
| **性别特异性** | 6mF Rapa 的指纹 | 找到雌性特异性抗衰老药物 |

### 3.5 与现有筛选结果的整合

我们已有 `results/ht_screen_metformin_like_candidates.csv`。可以：

1. 计算这些 Metformin-like 候选与 CR/Rapamycin 指纹的相似度
2. 识别"双重命中"：既像 Metformin 又像 CR 的化合物
3. 这些可能是比单一机制药物更优的抗衰老候选

### 3.6 预期成果

- 50-100 个新的 CR/Rapamycin 模拟物候选
- 按机制分类（mTOR-like, AMPK-like, 混合机制）
- 可验证的预测：在 GEO 中搜索这些候选药物的已有数据

### 3.7 风险评估

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| LINCS 体外签名不能预测体内效果 | 高 | 高 | 用 GEO 验证已有数据的候选 |
| 模块相似度不代表生物学相似 | 中 | 中 | 结合 pathway enrichment 分析 |
| 候选药物多为已研究化合物 | 中 | 低 | 关注新适应症而非新分子 |

---

## 三个方向的优先级建议

| 优先级 | 方向 | 工作量 | 风险 | 预期价值 |
|--------|------|--------|------|---------|
| **P0** | 模块级分析 | 2-3 小时 | 低 | **直接回答 mTOR 悖论机制** |
| **P1** | 高通量筛选整合 | 3-4 小时 | 中 | **发现新候选药物** |
| **P2** | ARCHS4 扩展 | 1-2 天 | 高 | **大规模验证/发现** |

### 推荐顺序

1. **先做模块分析**（最快，风险最低，直接回答核心问题）
2. **再做筛选整合**（利用已有数据，快速产出候选列表）
3. **最后做 ARCHS4**（需要更多时间和存储，但潜在回报最大）

---

*方案生成时间: 2026-06-09*
