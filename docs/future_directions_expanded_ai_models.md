# 转录组 × 衰老：AI 基础模型与虚拟细胞方向深入展开

> **撰写日期**: 2026-06-08
> **主题**: 方向四（基础模型）与方向八（虚拟细胞）的深度技术展开
> **定位**: 纯计算实现路径，包含具体代码框架与模型选择建议

---

## 目录

1. [方向四：基础模型驱动的衰老研究](#方向四基础模型驱动的衰老研究)
   - 1.1 三大模型架构对比
   - 1.2 GenePT：零成本的文本语义分析
   - 1.3 Geneformer：从预训练到 In Silico Perturbation
   - 1.4 scGPT：生成式扰动预测与衰老微调
   - 1.5 整合策略：三模型协同分析框架
2. [方向八：虚拟细胞与 In Silico 扰动筛选](#方向八虚拟细胞与-in-silico-扰动筛选)
   - 2.1 生成式模型技术路线对比
   - 2.2 推荐路线：scGen → CPA → Lingshu-Cell 渐进式升级
   - 2.3 与 tAge 模型的闭环验证管线
   - 2.4 从虚拟细胞到虚拟组织：空间扩展
3. [实施路线图与资源需求](#实施路线图与资源需求)
4. [风险与应对](#风险与应对)

---

## 方向四：基础模型驱动的衰老研究

### 1.1 三大模型架构对比

| 特性 | **GenePT** | **Geneformer** | **scGPT** |
|------|-----------|----------------|-----------|
| **核心思想** | LLM 文本嵌入基因语义 | Transformer 学习基因调控网络 | GPT 风格生成式细胞表示 |
| **预训练数据** | NCBI 基因描述文本 (~20K 基因) | ~104M 人类单细胞转录组 | ~33M 细胞，441 个研究 |
| **输入格式** | 基因文本描述 → 嵌入向量 | 基因排名编码 (rank-value encoding) | 离散化表达值 + 基因标识 |
| **架构** | text-embedding-3-large / GPT-3.5 | BERT-style Transformer (12-20 层) | GPT-style Transformer (12 层) |
| **最新版本** | GenePT-Large | Geneformer-V2-316M (2024.12) | checkpoint-0623 |
| **GPU 需求** | ❌ 无 (API 调用或本地轻量模型) | ✅ 需要 (316M 参数) | ✅ 需要 (512-dim embedding) |
| **关键优势** | 零训练成本；捕获文献级基因关系 | 网络层次结构编码在注意力权重中 | 生成能力 + 多组态整合 |
| **关键局限** | 无表达量信息；仅语义 | 非生成式；静态表示 | 需要大量微调数据 |
| **GitHub** | [GenePT](https://github.com/yiqunchen/GenePT) | [Geneformer](https://huggingface.co/ctheodoris/Geneformer) | [scGPT](https://github.com/bowang-lab/scGPT) |

**核心洞察**：这三个模型不是竞争关系，而是**互补**的——GenePT 提供"文献知识"，Geneformer 提供"网络结构"，scGPT 提供"生成能力"。最佳策略是**三管齐下**。

---

### 1.2 GenePT：零成本的文本语义分析

#### 技术原理

GenePT (Chen & Zou, 2024) 的核心 idea 极其简单但有效：
1. 对每个基因，从 NCBI 获取其功能描述文本
2. 用 OpenAI 的 `text-embedding-3-large` 将文本转为 3072-dim 向量
3. 对单个细胞，其嵌入 = 所有表达基因的加权平均 (权重 = 表达量)

**为什么这能工作？** LLM 在预训练时已经"读过"海量生物学文献，基因描述的语义相似度反映了功能相似度。

#### 衰老方向的零样本应用

**应用 1：tAge 特征的语义聚类**
```python
# 伪代码框架
import numpy as np
import openai

# 1. 加载 GenePT 基因嵌入 (预计算)
gene_embeddings = load_genept_embeddings()  # shape: (n_genes, 3072)

# 2. 取 tAge 模型 top 500 特征基因
top_genes = tage_model.get_top_features(k=500)

# 3. 对这些基因做文本语义聚类
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

embeddings = np.stack([gene_embeddings[g] for g in top_genes])
clusters = KMeans(n_clusters=10).fit_predict(embeddings)

# 4. 对每个 cluster，查看基因名称 → 人工/LLM 标注通路
for c in range(10):
    cluster_genes = [top_genes[i] for i, cl in enumerate(clusters) if cl == c]
    print(f"Cluster {c}: {cluster_genes[:10]}")
    # 预期发现: 某些 cluster 自动对应 "mitochondrial", "immune", "sirtuin" 等
```

**应用 2：衰老 vs 年轻细胞的语义距离**
```python
# 将 bulk 样本转换为 GenePT 细胞嵌入
def sample_to_genept_embedding(expression_vector, gene_names):
    """expression_vector: 样本的基因表达值 (TPM/logCPM)"""
    # 加权平均
    weights = expression_vector / expression_vector.sum()
    cell_embedding = np.average(
        [gene_embeddings[g] for g in gene_names],
        axis=0,
        weights=weights
    )
    return cell_embedding

# 计算年轻 vs 年老样本的语义距离
young_emb = sample_to_genept_embedding(young_expr, gene_names)
old_emb = sample_to_genept_embedding(old_expr, gene_names)
cosine_distance = 1 - cosine_similarity(young_emb, old_emb)
# 高 cosine_distance = 年轻和年老在语义空间中相距很远
```

**应用 3：药物扰动的语义预测 (零样本)**
```python
# 给定药物靶点基因，预测其对衰老的影响
drug_target_embedding = gene_embeddings["MTOR"]  # Rapamycin 靶点
aging_signature_embedding = average(tAge_top_genes_embeddings)

# 计算药物靶点与衰老特征的语义相关性
similarity = cosine_similarity(drug_target_embedding, aging_signature_embedding)
# 如果 MTOR 的语义与衰老特征高度相关 → 抑制 MTOR 可能影响衰老
```

#### 数据可用性

- **GenePT 预计算嵌入**: 可从 GitHub 仓库下载 (`genept_embeddings.npy`)
- **成本**: 如果下载预计算版本，**零 API 成本**；如果自己重新生成，OpenAI API 费用约 $5-10
- **时间**: 从下载到首次分析，< 1 小时

---

### 1.3 Geneformer：从预训练到 In Silico Perturbation

#### 技术架构详解

Geneformer 使用 **Rank Value Encoding**：
1. 对每个细胞，按表达量对基因排序 (不保存具体表达值，只保存排名)
2. 将排名序列输入 BERT-style Transformer
3. 预训练目标：mask 15% 的基因位置，预测被 mask 的基因是什么

**关键设计**：通过排名编码，模型学会了**基因之间的相对关系**而非绝对表达量——这使其对批次效应更鲁棒。

Geneformer-V2-316M (2024年12月发布) 的重大升级：
- 参数量从 10M → 316M
- 训练数据从 30M → 104M 细胞
- 输入大小从 2048 → 4096 基因
- 词汇表从 ~25K → ~20K 蛋白编码基因 (更干净)

#### 衰老方向的实施路径

**Phase 1：Zero-shot 细胞嵌入分析**
```python
# 安装
# git lfs install
# git clone https://huggingface.co/ctheodoris/Geneformer
# cd Geneformer && pip install .

from geneformer import EmbExtractor
import anndata as ad

# 将我们的 bulk 数据转换为伪单细胞格式 (Geneformer 要求单细胞输入)
# 策略: 把 bulk 样本当作一个"大细胞"
adata = ad.AnnData(X=expression_matrix)  # shape: (n_samples, n_genes)
adata.obs["sample_id"] = sample_ids
adata.obs["age_group"] = ["young" if a < 30 else "old" for a in ages]

# 使用预训练模型提取嵌入
emb_extractor = EmbExtractor(
    model_type="Geneformer",
    num_layers=12,
    emb_mode="cell",
    forward_batch_size=8,
    nproc=4
)

embeddings = emb_extractor.extract_embs(
    model_directory="geneformer-316M",  # 预训练权重路径
    input_data_file="our_data.tokenized.dataset",  # 需要先做 tokenization
    output_directory="./geneformer_outputs/",
    output_prefix="aging_analysis"
)

# 分析: 年轻 vs 年老在嵌入空间中的分离度
from sklearn.manifold import UMAP
import umap
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1)
embedding_2d = reducer.fit_transform(embeddings)

# 可视化: 年轻和年老是否形成明显 cluster？
# 如果是 → 模型 zero-shot 就捕捉到了衰老信号
```

**Phase 2：衰老分类微调**
```python
from geneformer import Classifier

# 微调任务: 预测细胞年龄组 (young / middle / old)
classifier = Classifier(
    classifier="cell",
    cell_state_dict={"age_group": {"young": 0, "middle": 1, "old": 2}},
    filter_data=None,
    num_crossval_splits=5,
    training_args={
        "learning_rate": 5e-5,
        "num_train_epochs": 10,
        "per_device_train_batch_size": 8,
        "seed": 42,
    }
)

classifier.prepare_data(
    input_data_file="our_aging_data.tokenized.dataset",
    output_directory="./geneformer_finetuned/",
    output_prefix="age_classifier"
)

# 关键决策: 冻结多少层？
# - 冻结前 6 层 (保留通用网络知识) + 微调后 6 层 (学习衰老特异性)
# - 或全部微调 (如果数据量 > 10K 细胞)
```

**Phase 3：In Silico Perturbation (ISP) — 核心应用**

这是 Geneformer 最独特的功能：模拟基因敲除/过表达对细胞状态的影响。

```python
from geneformer import InSilicoPerturber, InSilicoPerturberStats

# 1. 使用微调后的衰老分类器
isp = InSilicoPerturber(
    model_type="Geneformer",
    num_layers=12,
    emb_mode="cell",
    model_directory="./geneformer_finetuned/age_classifier/",  # 微调后的模型
)

# 2. 对一个"年老"细胞，模拟抑制 MTOR 的效果
perturbed_emb = isp.perturb(
    input_data_file="old_cells.tokenized.dataset",
    perturb_type="knockdown",  # 或 "overexpress"
    perturb_rank=500,  # 将 MTOR 的排名移到 500 位 (模拟低表达)
    genes_to_perturb=["MTOR"],
    output_directory="./isp_outputs/",
    output_prefix="mtor_kd_in_old_cells"
)

# 3. 量化扰动效果: 细胞嵌入向"年轻"方向移动了多少？
original_emb = isp.extract_embs("old_cells.tokenized.dataset")
young_centroid = np.mean(isp.extract_embs("young_cells.tokenized.dataset"), axis=0)
old_centroid = np.mean(original_emb, axis=0)
perturbed_centroid = np.mean(perturbed_emb, axis=0)

# 计算向年轻方向移动的比例
original_distance = np.linalg.norm(old_centroid - young_centroid)
perturbed_distance = np.linalg.norm(perturbed_centroid - young_centroid)
rejuvenation_ratio = (original_distance - perturbed_distance) / original_distance
# rejuvenation_ratio > 0 → MTOR 抑制使细胞"更年轻"
```

**关键洞察**: ISP 不是预测具体表达值，而是预测**细胞嵌入的移动方向**。通过与微调后的衰老分类器结合，我们可以将嵌入移动"翻译"为衰老状态变化。

**Phase 4：大规模筛选"虚拟抗衰老靶点"**
```python
# 对所有 ~20K 基因做 ISP，筛选使年老细胞向年轻方向移动的基因
isp_stats = InSilicoPerturberStats(
    model_directory="./geneformer_finetuned/age_classifier/"
)

# 批量扰动分析
results = isp_stats.perturb_all_genes(
    input_data_file="old_cells.tokenized.dataset",
    perturb_type="knockdown",
    output_directory="./genome_wide_isp/"
)

# 筛选标准:
# 1. perturbation 使细胞嵌入向 young centroid 移动
# 2. 统计显著性 (p < 0.05, 基于多次采样)
# 3. 与我们的 tAge 模型系数方向一致 (交叉验证)
anti_aging_targets = results.filter(
    direction="toward_young",
    p_value < 0.05,
    consistency_with_tage=True
)
```

#### 资源需求

| 组件 | 需求 | 说明 |
|------|------|------|
| GPU | A100 40GB 或 V100 32GB | 316M 模型需要大显存 |
| 存储 | ~20GB | 预训练权重 + tokenized 数据 |
| 时间 | Zero-shot: 数小时; 微调: 1-2 天 | 取决于数据量 |
| 数据 | 单细胞数据最佳; bulk 需伪单细胞化 | 衰老标注的 scRNA-seq |

---

### 1.4 scGPT：生成式扰动预测与衰老微调

#### 技术架构详解

scGPT 采用 **GPT 风格生成式架构**：
1. 将基因表达值离散化为 51 个 bins (类似 NLP 的 token)
2. 使用 causal mask 的自注意力 (只 attend 到前面的基因)
3. 预训练目标：给定部分基因，预测下一个基因的表达 bin

**与 Geneformer 的关键区别**：
- Geneformer 是 bidirectional (BERT)，scGPT 是 autoregressive (GPT)
- scGPT 明确建模了**零表达概率** (zero-inflation)，这对 scRNA-seq 至关重要
- scGPT 支持**多组态** (RNA + ATAC + protein)

#### 衰老方向的实施路径

**Phase 1：Zero-shot 细胞嵌入与参考映射**

scGPT 提供了一个强大的 zero-shot 功能：将查询细胞映射到参考图谱 (CellXGene 图谱)。

```python
# 安装
# pip install scgpt

import scgpt
from scgpt.tasks import embed_data
import scanpy as sc

# 加载预训练模型 (whole-human model)
model_dir = "./scgpt_checkpoint/"

# 准备数据: 将我们的 bulk 数据转换为 AnnData
# scGPT 要求 raw counts (不要 normalized!)
adata = sc.AnnData(X=raw_count_matrix)
adata.var["gene_name"] = gene_symbols

# 提取嵌入
embedded_adata = embed_data(
    adata,
    model_dir,
    gene_col="gene_name",
    batch_size=64,
    device="cuda"
)

# embedded_adata.obsm["X_scGPT"] 包含细胞嵌入 (512-dim)
```

**Phase 2：衰老连续值回归微调**

与 Geneformer 的分类不同，scGPT 更适合做**连续年龄回归** (chronological age / biological age)。

```python
from scgpt import TransformerModel
from scgpt.trainer import SeqTrainer
import torch

# 加载预训练模型
model = TransformerModel(
    ntokens=len(vocab),
    embsize=512,
    nhead=8,
    d_hid=512,
    nlayers=12,
    vocab=vocab,
    # ... 其他参数
)

# 加载预训练权重
checkpoint = torch.load(f"{model_dir}/best_model.pt")
model.load_state_dict(checkpoint["model_state_dict"])

# 修改输出头: 从分类改为回归
# 在 scGPT 的 [CLS] token 上接一个线性层，输出预测年龄
model.age_regressor = torch.nn.Linear(512, 1)  # 输出连续年龄

# 训练配置
training_args = {
    "learning_rate": 1e-4,
    "batch_size": 32,
    "epochs": 20,
    "device": "cuda",
    "freeze_encoder": True,  # 冻结 Transformer 编码器，只训练回归头
}

# 训练数据: scRNA-seq 细胞 + 真实年龄标签
trainer = SeqTrainer(
    model=model,
    train_dataset=aging_train_dataset,
    val_dataset=aging_val_dataset,
    args=training_args
)
trainer.train()
```

**Phase 3：扰动预测 (Perturbation Prediction)**

scGPT 官方提供了**扰动预测微调**的 tutorial。

```python
# 参考: scGPT/tutorials/Tutorial_Perturbation.ipynb

from scgpt.tasks import perturbation_prediction

# 数据格式要求:
# adata.obs["condition"] = "ctrl" or "perturbation_name"
# adata.obs["cell_type"] = cell_type_label

# 微调模型预测扰动响应
perturbation_prediction(
    model=model,
    adata=perturbation_train_data,  # 包含 ctrl + perturbation 的 scRNA-seq
    config={
        "max_epochs": 15,
        "batch_size": 64,
        "lr": 1e-4,
        "save_dir": "./scgpt_perturbation_finetuned/"
    }
)

# 预测新扰动 (zero-shot)
# 输入: 未扰动的细胞 + 扰动条件
# 输出: 预测的扰动后表达谱
predicted_perturbed = model.predict_perturbation(
    control_cells=ctrl_adata,
    perturbation="Rapamycin",  # 或基因名如 "MTOR"
    cell_type="hepatocyte"
)
```

**关键应用：预测药物组合效应**

我们发现 Trametinib + Rapamycin 组合有协同效应。scGPT 可以**预测未测试的组合**：

```python
# 训练数据: 单药扰动 (Trametinib, Rapamycin)
# 测试: 预测组合扰动 (Trametinib + Rapamycin)

predicted_combo = model.predict_combinatorial_perturbation(
    control_cells=ctrl_adata,
    perturbations=["Trametinib", "Rapamycin"],
    cell_type="hepatocyte"
)

# 与真实 combo 数据 (GSE288795) 比较验证
# 计算预测 vs 真实的 Pearson 相关性
```

**Phase 4：基因网络推断 (GRN Inference)**

scGPT 的注意力权重可以提取基因调控网络：

```python
from scgpt.tasks import grn_inference

# 从微调后的模型提取注意力权重
attention_matrix = grn_inference(
    model=model,
    adata=aging_scRNA_data,
    gene_names=gene_list,
    top_k=20  # 每个基因取 top-20 注意力权重邻居
)

# 构建衰老特异性 GRN
# 对比: 年轻细胞 vs 年老细胞的 attention pattern 差异
# 识别"衰老 rewiring"事件——即调控关系随年龄变化
```

---

### 1.5 整合策略：三模型协同分析框架

单一模型有局限，但三者结合可形成强大的分析 pipeline：

```
┌─────────────────────────────────────────────────────────────────────┐
│                    三模型协同衰老分析 Pipeline                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  Step 1: GenePT 快速筛选 (零成本)                                    │
│  ───────────────────────────────                                    │
│  • 对 tAge top 500 基因做语义聚类                                     │
│  • 识别"文献已知"的衰老通路 vs "未知"cluster                          │
│  • 对候选药物靶点做语义相关性快速评估                                   │
│                                                                     │
│                    ↓ 输出: 优先级排序的假设列表                         │
│                                                                     │
│  Step 2: Geneformer 机制验证 (中等成本)                                │
│  ────────────────────────────────────                                │
│  • 在衰老 scRNA-seq 数据上微调分类器                                   │
│  • In Silico Perturbation 验证 Step 1 的假设                          │
│  • 量化扰动使细胞"年轻化"的程度                                       │
│                                                                     │
│                    ↓ 输出: 机制验证 + 效应量化                          │
│                                                                     │
│  Step 3: scGPT 生成式预测 (高成本)                                     │
│  ─────────────────────────────────                                   │
│  • 微调生成模型预测扰动后的完整转录组                                   │
│  • 用 tAge 模型评分预测转录组 → 计算 ΔtAge                            │
│  • 预测未测试的药物组合                                               │
│                                                                     │
│                    ↓ 输出: 可验证的定量预测                             │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

**协同示例**：筛选抗衰老靶点

1. **GenePT 发现**："SIRT1" 和 "FOXO3" 的文本嵌入与 tAge 特征高度相关 (cosine sim > 0.7)
2. **Geneformer 验证**：ISP 显示 SIRT1 过表达使年老肝细胞嵌入向年轻方向移动 15%
3. **scGPT 量化**：生成 SIRT1-OE 后的预测转录组 → tAge 评分 → 预测 tAge 降低 23 个月

---

## 方向八：虚拟细胞与 In Silico 扰动筛选

### 2.1 生成式模型技术路线对比

| 模型 | 架构 | 训练数据 | 药物预测 | 基因预测 | 组合预测 | 代码成熟度 |
|------|------|---------|---------|---------|---------|-----------|
| **scGen** | VAE + latent arithmetic | 单细胞 | ✅ | ✅ | ❌ (纯相加) | ⭐⭐⭐⭐⭐ |
| **CPA** | 条件 VAE | 单细胞/ bulk | ✅ | ✅ | ✅ (相加) | ⭐⭐⭐⭐⭐ |
| **CellOT** | Neural Optimal Transport | 单细胞 | ✅ | ✅ | ❌ | ⭐⭐⭐⭐ |
| **scDiffusion** | Diffusion Model | 单细胞 | ❌ | ❌ | ❌ | ⭐⭐⭐ |
| **Lingshu-Cell** | Masked Discrete Diffusion | 单细胞 (18K 基因) | ❌ | ✅ | ✅ | ⭐⭐ (2026.03 新发) |
| **scDFM** | Flow Matching | 单细胞 | ✅ | ✅ | ✅ | ⭐⭐⭐ (2026.02) |

**推荐渐进式路线**：
- **入门**: scGen (最简单，有完整 tutorial)
- **进阶**: CPA (支持组合预测，有 dose 信息)
- **前沿**: Lingshu-Cell (最新 SOTA，但代码可能不稳定)

---

### 2.2 推荐路线：scGen → CPA → Lingshu-Cell 渐进式升级

#### Level 1: scGen — 快速验证概念

**核心思想**：VAE 将细胞映射到低维 latent space，扰动效应 = latent space 中的向量加法。

```python
# 安装: pip install scgen (或通过 pertpy)
import pertpy as pt
import scanpy as sc

# 准备数据: 需要单细胞水平的 ctrl + perturbation
# 数据来源: 公开 scRNA-seq 扰动数据集 (如 Kang 2018, Norman 2019)
# 或我们的 bulk 数据经过去卷积后的伪单细胞

adata = sc.read_h5ad("aging_perturbation_scRNA.h5ad")
# adata.obs 需要包含:
#   - "condition": "ctrl" or "Rapamycin" or "CR_simulated"
#   - "cell_type": "hepatocyte", "fibroblast", etc.
#   - "age_group": "young" or "old"

# 训练 scGen
pt.tl.SCGEN.setup_anndata(adata, batch_key="condition", labels_key="cell_type")
model = pt.tl.SCGEN(adata)
model.train(max_epochs=100, batch_size=32, early_stopping=True)

# 预测: 将年轻细胞的扰动效应转移到年老细胞
# 先学习"年轻 hepatocyte + Rapamycin"的 latent shift
pred_young, delta_young = model.predict(
    ctrl_key="ctrl",
    stim_key="Rapamycin",
    celltype_to_predict="hepatocyte"
)

# 将同样的 delta 应用到年老细胞
old_hepatocytes = adata[
    (adata.obs["cell_type"] == "hepatocyte") & 
    (adata.obs["age_group"] == "old") &
    (adata.obs["condition"] == "ctrl")
]
predicted_old_treated = model.predict_similar(
    control_adata=old_hepatocytes,
    delta=delta_young  # 使用年轻细胞的扰动向量
)

# 用 tAge 评分
old_tAge = tage_predict(old_hepatocytes.X)
predicted_tAge = tage_predict(predicted_old_treated.X)
print(f"预测 tAge 变化: {predicted_tAge - old_tAge} months")
```

**scGen 的局限**：
- 假设扰动效应是**可加的**且**与细胞状态无关**
- 不支持组合预测 (组合 = 简单向量相加，不考虑非线性交互)
- 实际上我们发现 Trametinib + Rapamycin 有非线性协同效应

#### Level 2: CPA — 组合预测与剂量建模

**核心升级**：CPA 显式建模了扰动条件 (药物/基因)、细胞类型、剂量作为**可组合的条件嵌入**。

```python
# CPA 支持: perturbation embedding + cell_type embedding + dose embedding
# 组合预测: 将两个条件嵌入相加

# 关键应用: 预测药物组合对衰老的影响
from cpa import CPA

# 训练数据需要包含剂量信息
model = CPA(
    adata=adata,
    perturbation_key="condition",  # "ctrl", "Rapamycin", "Trametinib", "Combo"
    dose_key="dose",               # 0.0, 1.0, 5.0, 10.0 (μM)
    cell_type_key="cell_type",
    # ...
)

model.train(max_epochs=200)

# 预测组合效应 (zero-shot)
combo_prediction = model.predict(
    cell_type="hepatocyte",
    perturbations=["Rapamycin", "Trametinib"],  # 组合!
    doses=[1.0, 0.5],
    age_group="old"
)

# 与 scGen 的简单相加比较:
# CPA 学习到的组合嵌入可能 ≠ Rapamycin_emb + Trametinib_emb
# 这捕捉了非线性协同/拮抗效应
```

#### Level 3: Lingshu-Cell — 最新 SOTA

**Lingshu-Cell** (2026.03, 阿里巴巴达摩院) 是目前最先进的虚拟细胞模型：

**技术亮点**：
- **Masked Discrete Diffusion**: 在离散 token 空间中操作，完美匹配 scRNA-seq 的稀疏性
- **18,080 基因全基因组**: 不需要预选 HVG 或排名
- **条件生成**: 支持细胞类型 + 扰动的联合条件
- **组合泛化**: 在 VCC H1 benchmark 上取得 SOTA

```python
# Lingshu-Cell (概念性代码，实际 API 可能不同)
from lingshu_cell import LingshuCell

model = LingshuCell.from_pretrained("lingshu-cell-base")

# 条件生成: 给定细胞和扰动，预测扰动后状态
predicted_cells = model.generate(
    cell_type="hepatocyte",
    donor_age="25_months",  # 小鼠
    perturbation="Mtor_knockdown",  # 基因扰动
    n_cells=1000
)

# 也可以做药物扰动
predicted_cells_drug = model.generate(
    cell_type="hepatocyte",
    perturbation="Rapamycin_1uM",  # 药物扰动
    n_cells=1000
)

# 评分
for pred in [predicted_cells, predicted_cells_drug]:
    tAge_score = tage_predict(pred.to_bulk())
    print(f"预测 tAge: {tAge_score}")
```

**当前问题**：Lingshu-Cell 2026 年 3 月刚发预印本，代码和权重可能尚未公开。建议**密切关注其 GitHub 仓库**。

---

### 2.3 与 tAge 模型的闭环验证管线

虚拟细胞的核心价值不是生成"逼真的"转录组，而是**生成后能被 tAge 评分，且评分变化与真实数据一致**。

```
┌─────────────────────────────────────────────────────────────────────┐
│                  虚拟细胞 - tAge 闭环验证管线                          │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  输入: 未扰动细胞 (年轻/年老, 组织, 细胞类型)                          │
│                    ↓                                                │
│  ┌─────────────────────────────────────┐                            │
│  │    生成式模型 (scGen/CPA/Lingshu)    │                            │
│  │                                     │                            │
│  │  条件: "模拟药物 X 处理 24h"          │                            │
│  │  输出: 预测的扰动后单细胞转录组        │                            │
│  └─────────────────────────────────────┘                            │
│                    ↓                                                │
│  步骤 A: Pseudobulk 聚合 → 模拟 bulk 表达谱                            │
│  (或保持单细胞，构建单细胞 tAge 评分)                                   │
│                    ↓                                                │
│  步骤 B: tAge 模型评分                                                │
│  tAge_pred = tage_predict(predicted_pseudobulk)                     │
│                    ↓                                                │
│  步骤 C: 与真实数据比较 (如果有的话)                                    │
│  ┌─────────────────────────────────────────────────────────────┐    │
│  │  验证指标:                                                    │    │
│  │  1. 方向一致性: sign(ΔtAge_pred) == sign(ΔtAge_real)?        │    │
│  │  2. 幅度相关性: corr(ΔtAge_pred, ΔtAge_real) across drugs    │    │
│  │  3. 排名一致性: Spearman rank corr between predicted and real │    │
│  └─────────────────────────────────────────────────────────────┘    │
│                    ↓                                                │
│  步骤 D: 模型校准 (如果系统性偏差)                                     │
│  如果预测总是高估/低估 → Isotonic Regression 校准                     │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

**关键验证数据集**（我们已有）：

| 数据集 | 真实 ΔtAge | 可用于验证 |
|--------|-----------|-----------|
| GSE288795 (Trametinib+Rapa) | Combo -53.3mo | ✅ |
| GSE230402 (CR) | Female -48.0mo, Male -20.6mo | ✅ |
| GSE280382 (GLP-1RA) | +23.8mo (悖论) | ✅ |
| ARCHS4 MR studies | -198.9mo | ⚠️ 作为极端值测试 |

**预期性能基准**：
- 方向一致性 > 80% (预测抗衰老的药物真实也抗衰老)
- 排名相关性 > 0.6 (预测的抗衰老强度排名与真实排名一致)
- 如果达不到 → 说明生成模型没有准确捕捉衰老相关基因的变化模式

---

### 2.4 从虚拟细胞到虚拟组织：空间扩展

终极愿景：不仅预测单个细胞的响应，还预测**组织中衰老微环境的空间演化**。

```python
# 概念性框架 (目前技术前沿)

# 步骤 1: 用 Stereo-seq/Visium 数据构建"空间衰老图谱"
spatial_adata = sc.read_h5ad("stereo_seq_liver_aging.h5ad")
# spatial_adata.obsm["spatial"] = (x, y) 坐标

# 步骤 2: 对每个 spatial spot 生成虚拟细胞群体
# 每个 spot 包含多个细胞，模型生成混合群体
spot_predictions = []
for spot_id in spatial_adata.obs_names:
    spot_cells = model.generate(
        cell_type_mixture=infer_from_spot(spot_id),  # 去卷积推断的细胞组成
        perturbation="Rapamycin",
        spatial_context=spot_id  # 考虑邻域信号
    )
    spot_predictions.append(spot_cells)

# 步骤 3: 计算空间 tAge 变化热图
spatial_tAge_change = compute_spatial_tAge(spot_predictions) - compute_spatial_tAge(spatial_adata)
# 识别: 药物是否优先逆转某些空间区域的衰老？
```

**当前限制**：空间虚拟细胞仍是研究前沿，但可以先用已发表的 Stereo-seq 衰老数据 (Ma et al., 2024, Cell) 做概念验证。

---

## 实施路线图与资源需求

### 推荐的两阶段实施计划

#### 第一阶段 (4–6 周)：快速验证概念

| 周 | 任务 | 模型 | 产出 |
|----|------|------|------|
| 1 | GenePT 语义分析 tAge 特征 | GenePT | 衰老基因语义聚类报告 |
| 1–2 | 下载 Geneformer/scGPT 权重，跑通 zero-shot 嵌入 | Geneformer + scGPT | 年轻 vs 年老嵌入分离度分析 |
| 2–3 | 在公开衰老 scRNA-seq 数据上微调 Geneformer 分类器 | Geneformer | 衰老分类模型 |
| 3–4 | In Silico Perturbation 筛选 top 候选基因 | Geneformer | 50 个虚拟抗衰老靶点 |
| 4–5 | scGen 入门：用 Kang 2018 数据跑通 tutorial | scGen | 扰动预测基础能力 |
| 5–6 | 将 scGen 应用于衰老数据，与 tAge 交叉验证 | scGen + tAge | 虚拟细胞-tAge 闭环验证报告 |

**资源需求 (第一阶段)**：
- GPU: 1× A100 40GB (或 2× RTX 4090 24GB)
- 存储: ~100GB
- 数据: 全部公开下载
- 人力: 1 人全职

#### 第二阶段 (6–8 周)：深度建模

| 周 | 任务 | 模型 | 产出 |
|----|------|------|------|
| 7–9 | CPA 训练：支持组合预测和剂量建模 | CPA | 药物组合效应预测模型 |
| 9–11 | scGPT 微调：衰老回归 + 扰动预测 | scGPT | 生成式衰老预测模型 |
| 11–13 | 大规模虚拟筛选：1000+ 化合物 | scGPT/CPA + tAge | 排名前 20 的虚拟抗衰老候选 |
| 13–14 | 与 HT 筛选结果整合，重排序 | 全部 | 整合优先级候选列表 |

**资源需求 (第二阶段)**：
- GPU: 2–4× A100 40GB (多卡并行)
- 存储: ~500GB
- 人力: 1–2 人

---

## 风险与应对

| 风险 | 概率 | 影响 | 应对策略 |
|------|------|------|---------|
| GPU 显存不足 (316M 模型) | 中 | 高 | 使用梯度累积 + 混合精度 (fp16); 或先用 104M 版本 |
| 预训练权重下载失败 | 低 | 高 | 提前用 `git-lfs` 克隆; 准备 HuggingFace mirror |
| scRNA-seq 衰老数据质量差 | 中 | 中 | 用 Tabula Muris Senis 作为主要数据源; 做严格 QC |
| 生成模型预测不准确 | 高 | 中 | 与简单基线 (线性模型) 比较; 不过度解读单个预测 |
| 代码/依赖版本冲突 | 高 | 低 | 用 Docker/conda 隔离环境; 锁定版本 |
| Lingshu-Cell 代码未公开 | 中 | 中 | 以 CPA 作为备选方案; 保持关注 |

---

## 附录：关键代码片段

### A. Geneformer Tokenization (Bulk → Pseudo-single-cell)

```python
from geneformer import TranscriptomeTokenizer
import anndata as ad

# Bulk 数据 → 伪单细胞 (每个样本 = 一个"细胞")
adata = ad.AnnData(X=expression_matrix.T)  # genes x samples → samples x genes
adata.obs["sample_id"] = sample_ids
adata.obs["age"] = ages

# Tokenization (将表达转为排名编码)
tk = TranscriptomeTokenizer(
    {"age": "age"},  # 元数据列
    nproc=4
)
tk.tokenize_data(
    adata,
    output_directory="./tokenized_data/",
    output_prefix="aging_bulk",
    file_format="loom"  # 或 "h5ad"
)
```

### B. scGPT 数据预处理

```python
import scanpy as sc
from scgpt.preprocess import Preprocessor

# scGPT 要求 raw counts
adata = sc.read_h5ad("raw_counts.h5ad")

preprocessor = Preprocessor(
    use_key="X",  # raw counts
    filter_gene_by_counts=3,
    filter_cell_by_counts=False,
    normalize_total=1e4,
    result_normed_key="X_normed",
    log1p=True,
    result_log1p_key="X_log1p",
    subset_hvg=3000,
    hvg_flavor="seurat_v3",
    binning=51,  # scGPT 关键：离散化为 51 bins
    result_binned_key="X_binned",
)
preprocessor(adata, batch_key=None)
```

### C. 虚拟细胞生成后的 tAge 评分

```python
import numpy as np

# 假设 predicted_cells 是 scGen/CPA 生成的 AnnData (cells x genes)
# 需要聚合为 pseudobulk 才能输入 tAge 模型

def sc_to_pseudobulk(adata, groupby="condition"):
    """将单细胞数据聚合为 pseudobulk"""
    pseudobulk = []
    for group in adata.obs[groupby].unique():
        cells = adata[adata.obs[groupby] == group]
        bulk_expr = np.array(cells.X.mean(axis=0)).flatten()
        pseudobulk.append(bulk_expr)
    return np.vstack(pseudobulk)

# 预测扰动前后的 pseudobulk
ctrl_bulk = sc_to_pseudobulk(ctrl_cells)
pred_bulk = sc_to_pseudobulk(predicted_cells)

# tAge 评分
ctrl_tAge = tage_predict(ctrl_bulk)
pred_tAge = tage_predict(pred_bulk)

delta_tAge = pred_tAge - ctrl_tAge
print(f"预测抗衰老效应: {delta_tAge} months")
```

---

## 参考文献

1. Chen, Y. & Zou, J. GenePT: a simple but effective foundation model for genes and cells built from ChatGPT. *bioRxiv* (2024).
2. Theodoris, C.V. et al. Transfer learning enables predictions in network biology. *Nature* (2023).
3. Cui, H. et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods* (2024).
4. Yuan, G.H. et al. Lingshu-Cell: A generative cellular world model for transcriptome modeling toward virtual cells. *arXiv* (2026).
5. Lotfollahi, M. et al. scGen predicts single-cell perturbation responses. *Nature Methods* (2019).
6. Roohani, Y. et al. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nature Biotechnology* (2024).
7. Bunne, C. et al. Learning single-cell perturbation responses using neural optimal transport. *Nature Methods* (2023).
8. Istrate, A.M. et al. scGenePT: Is language all you need for modeling single-cell perturbations? *bioRxiv* (2024).
9. Bicklab. Closing the loop: Teaching single-cell foundation models to learn from perturbations. *GitHub* (2024).
10. Yu, T. et al. scDFM: Distributional Flow Matching for Robust Single-Cell Perturbation Prediction. *arXiv* (2026).
