# TranscriptFormer × Tabula Muris Senis: 纯小鼠衰老分析方案

> 方案版本: v1.0 | 状态: 待用户最终确认

---

## 1. 为什么考虑 TranscriptFormer

### 1.1 核心优势

| 属性 | TranscriptFormer | UCE（对比）|
|------|-----------------|-----------|
| 训练数据 | **112M cells, 12 species** | 36M cells, 8 species |
| 小鼠训练 | ✅ TF-Exemplar 包含 110M（人+鼠+3 模式生物） | ✅ 36M 中的一部分 |
| CGE | ✅ **Contextual Gene Embeddings**（每个基因在每个细胞中的表示） | ❌ 只有 cell-level |
| 代码状态 | 2025年4月开源，活跃维护 | 2023年，稳定 |
| CLI 体验 | PyPI 包，一行命令 | 需要 clone repo |
| 多 GPU | ✅ DDP 原生支持 | 需要手动处理 |

### 1.2 模型选择建议

对于小鼠 aging 任务，推荐 **TF-Exemplar**：
- 包含小鼠（in-distribution，不需要 OOD protein embeddings）
- 824M 总参数（542M trainable + 282M non-trainable）
- Vocabulary: 110,290
- 比 TF-Metazoa 更快、更聚焦

---

## 2. 关键风险（必须知情）

### 2.1 🔴 高风险

#### R1: PyTorch 版本冲突（确定发生）

- **TranscriptFormer 要求**: PyTorch <= 2.5.1（明确说 2.6.0+ 可能有 pickle errors）
- **当前环境**: PyTorch 2.12.0 + CUDA 12.6
- **影响**: 必须新建 conda 环境，单独安装 PyTorch 2.5.1 + CUDA 12.1/12.4
- **额外风险**: 新建环境可能与现有 `aging-fm` 环境的包产生冲突；需要管理两个环境

#### R2: Fine-tuning 不支持（确定）

- **官方状态**: Helical 文档明确标注 "Example Fine-Tuning: Coming soon!"
- **论文做法**: 只做了 linear probing（固定 embedding + LR）
- **影响**: 如果 zero-shot embedding 没有 capture aging signal，**没有退路**
- **与 UCE 对比**: UCE 也是不支持 fine-tuning，但 UCE 是设计选择（论文论证充分），TranscriptFormer 是还没实现

#### R3: Perturbation 非原生支持（确定）

- **官方声明**: "TranscriptFormer is not specialized for zero-shot perturbation prediction"
- **影响**: ISP 需要手动修改输入 counts → 重新 inference，没有专门 API
- **与 UCE 对比**: UCE 有 issue #57（其他细胞 embedding 会改变），TranscriptFormer 没有文档说明这个问题，但因为是 generative autoregressive model，同样存在 batch processing 副作用的可能

### 2.2 🟡 中风险

#### R4: 基因 ID 映射（高概率发生）

- **TranscriptFormer 要求**: `adata.var` 必须有 `ensembl_id` 列（ENSMUSG ID）
- **Tabula Muris Senis 实际**: MSigDB 确认使用 **MOUSE_GENE_SYMBOL**（如 `Actb`, `Cd3d`）
- **影响**: 需要 mygene.info / BioMart 映射，额外工作量；可能有 5-15% 基因无法映射
- **缓解**: 可以先下载数据确认实际格式

#### R5: 显存需求高（中概率）

- **TF-Exemplar**: 824M 参数，FP16 约 1.6GB 权重 + activations + KV cache
- **官方建议**: A100 40GB；16GB GPU 可用 batch_size 1-4
- **我们的硬件**: 4× RTX 3090 (24GB each)
- **估算**: batch_size 4-8 应该可行，但需要验证
- **影响**: 50万细胞推理时间可能较长（估算 6-12 小时，vs UCE 4-layer 可能 2-4 小时）

#### R6: 推理速度慢（高概率）

- **复杂度**: O(L(n²d + nd²)) per cell，每细胞 10¹¹-10¹² FLOPs
- **估算**: 50万细胞 × 4 GPU ≈ 6-12 小时（vs UCE 4-layer 可能 2-4 小时）
- **影响**: 快速验证 aging signal 的时间成本更高

### 2.3 🟢 低风险

#### R7: 安装复杂度
- PyPI 包可用，`pip install transcriptformer` 即可
- 但需要先解决 PyTorch 版本问题

#### R8: 数据下载
- Tabula Muris Senis droplet 3.5GB，Figshare 可 wget

---

## 3. 风险评估矩阵

| 风险 | 严重度 | 可能性 | 综合等级 | 缓解措施 |
|------|--------|--------|---------|---------|
| R1 PyTorch 冲突 | 高 | 确定 | **🔴 极高** | 新建 conda 环境 `transcriptformer` |
| R2 Fine-tuning 不支持 | 高 | 确定 | **🔴 极高** | 先小规模验证 aging signal；信号弱则汇报 |
| R3 Perturbation 不支持 | 中 | 确定 | **🟡 高** | 手动实现；结果保守解读 |
| R4 基因 ID 映射 | 中 | 高 | **🟡 高** | mygene.info 映射；先做 mapping rate 评估 |
| R5 显存 OOM | 中 | 中 | **🟡 中** | batch_size 4；OOM-safe dataloader |
| R6 推理速度慢 | 中 | 高 | **🟡 中** | 4 GPU DDP；夜间跑 |

---

## 4. Pipeline 设计

### Phase 0: 环境准备（0.5 天）

```bash
# 新建 conda 环境（避免和 aging-fm 冲突）
conda create -n transcriptformer python=3.11
conda activate transcriptformer

# 安装兼容的 PyTorch
pip install torch==2.5.1 torchvision==0.20.1 torchaudio==2.5.1 \
    --index-url https://download.pytorch.org/whl/cu124

# 安装 TranscriptFormer
pip install transcriptformer

# 验证
transcriptformer --version
```

### Phase 1: 数据准备 + 基因映射（0.5-1 天）

```bash
# 下载 Tabula Muris Senis
cd data/
wget https://figshare.com/ndownloader/files/44704373 -O tms_droplet.h5ad
```

```python
# 基因 symbol → ENSEMBL ID 映射
import mygene
import scanpy as sc

adata = sc.read_h5ad('data/tms_droplet.h5ad')
mg = mygene.MyGeneInfo()
result = mg.querymany(list(adata.var_names), scopes='symbol', species='mouse', 
                       fields='ensembl.gene', as_dataframe=True)
adata.var['ensembl_id'] = result['ensembl.gene'].reindex(adata.var_names).values

# 检查 mapping rate
mapping_rate = adata.var['ensembl_id'].notna().mean()
print(f"Mapping rate: {mapping_rate:.1%}")
# 如果 <85%，需要考虑替代方案
```

### Phase 2: TranscriptFormer Zero-Shot Embedding（1-2 天）

```bash
# 下载 TF-Exemplar 权重
transcriptformer download tf-exemplar --checkpoint-dir ./checkpoints

# 分组织推理（以脾为例）
transcriptformer inference \
  --checkpoint-path ./checkpoints/tf_exemplar \
  --data-file data/tms_spleen.h5ad \
  --output-path results/tf_embeddings/ \
  --output-filename spleen_cell_emb.h5ad \
  --batch-size 4 \
  --num-gpus 4 \
  --precision 16-mixed \
  --oom-dataloader \
  --emb-type cell

# 同时提取 CGE（用于后续 interpretation）
transcriptformer inference \
  --checkpoint-path ./checkpoints/tf_exemplar \
  --data-file data/tms_spleen.h5ad \
  --output-path results/tf_embeddings/ \
  --output-filename spleen_cge.h5ad \
  --batch-size 4 \
  --num-gpus 4 \
  --precision 16-mixed \
  --oom-dataloader \
  --emb-type cge
```

**输出**: 
- Cell embeddings: `obsm['embeddings']`
- CGE: `uns['cge_embeddings']`, `uns['cge_cell_indices']`, `uns['cge_gene_names']`

### Phase 3: Aging Classifier（0.5-1 天）

与 UCE 方案相同，在固定 embedding 上训练：
- 二分类：young (1m+3m) vs old (18m+)
- 连续回归：预测月龄
- 模型：ElasticNet / MLP / XGBoost
- CV：leave-one-tissue-out, leave-one-age-out

### Phase 4: Aging Axis + CGE Interpretation（0.5 天）

**独特优势（UCE 没有）**:
```python
# CGE 分析：哪些基因在 young vs old 中的 contextual embedding 差异最大？
import numpy as np

adata = sc.read_h5ad('results/tf_embeddings/spleen_cge.h5ad')
cge = adata.uns['cge_embeddings']      # (n_gene_instances, emb_dim)
cells = adata.uns['cge_cell_indices']   # 每个 gene instance 属于哪个细胞
genes = adata.uns['cge_gene_names']     # 每个 gene instance 的基因名

# 计算每个基因在 young vs old 中的平均 CGE
young_mask = adata.obs['age'] <= 3
old_mask = adata.obs['age'] >= 18

for gene in top_aging_genes:
    gene_mask = genes == gene
    young_emb = cge[gene_mask & young_mask[cells]].mean(axis=0)
    old_emb = cge[gene_mask & old_mask[cells]].mean(axis=0)
    delta = np.linalg.norm(old_emb - young_emb)
    print(f"{gene}: CGE delta = {delta:.3f}")
```

### Phase 5: Perturbation（探索性，0.5 天）

```python
# 手动 ISP：修改输入 counts → 重新 inference
adata_pert = adata.copy()
adata_pert[:, "GeneX"].X = 0

# 保存并重新 inference
adata_pert.write('data/tms_spleen_perturbed.h5ad')
# ... 运行 transcriptformer inference ...

# 比较 perturbed vs baseline 的 aging score
```

**注意**: 因为没有专门的 perturbation API，结果解读需保守。

### Phase 6: 与已知 Signatures 比较（0.5 天）

同 UCE 方案。

---

## 5. 时间线

| 阶段 | 时间 | 关键里程碑 |
|------|------|-----------|
| **Day 1 上午** | 0.5 天 | 新建 conda 环境，安装 TranscriptFormer，验证能跑 |
| **Day 1 下午** | 0.5 天 | 下载 TMS 数据，完成基因 symbol → ENSMUSG 映射 |
| **Day 2** | 1 天 | 1-2 个组织的 TF-Exemplar embedding，快速验证 aging signal |
| **Day 3** | 1 天 | 全部组织 embedding（4 GPU 并行，夜间跑） |
| **Day 4** | 1 天 | Aging classifier + CGE interpretation |
| **Day 5** | 0.5-1 天 | Perturbation + 与已知 signatures 比较 |

**总计: 4-5 天**（比 UCE 多 0.5-1 天，因为环境准备和基因映射）

---

## 6. 与 UCE 方案的直接对比

| 维度 | TranscriptFormer (TF-Exemplar) | UCE (4-layer → 33-layer) |
|------|-------------------------------|--------------------------|
| **训练数据量** | ✅ 110M (人+鼠) | 36M (8物种) |
| **环境兼容性** | ❌ **需要新建 conda 环境** PyTorch 降级 | ✅ 当前环境直接用 |
| **基因 ID 映射** | ❌ **需要 ENSMUSG** | ✅ 不需要（蛋白质序列） |
| **推理速度** | ⚠️ 慢（824M 参数） | ✅ 快（4-layer 很轻） |
| **显存需求** | ⚠️ 高（24GB batch_size 4-8） | ✅ 低（4-layer batch_size 50+） |
| **Fine-tuning** | ❌ Coming soon | ❌ 设计如此 |
| **Perturbation** | ❌ 不是专门设计 | ⚠️ issue #57 |
| **CGE** | ✅ **有** | ❌ 没有 |
| **快速验证** | ⚠️ 需要半天环境准备 | ✅ 立即可启动 |
| **aging signal 弱时的退路** | ❌ 无（不能 fine-tune） | ❌ 无（不能 fine-tune） |

**核心结论**: 两者在 aging 任务上的范式几乎相同（zero-shot embedding + 固定 classifier）。TranscriptFormer 的优势（更大训练数据、CGE）被其劣势（环境冲突、基因映射、更慢推理）抵消。对于"快速验证 aging signal 是否存在"这个首要目标，UCE 更稳妥。

---

## 7. 我的建议

### 7.1 如果坚持 TranscriptFormer

可以执行，但请知情以下额外成本：
1. **0.5 天环境准备**（PyTorch 降级 + 新建 conda 环境）
2. **0.5 天基因映射**（symbol → ENSMUSG，可能 5-15% loss）
3. **更高显存压力**（batch_size 4-8，推理更慢）
4. **如果 zero-shot aging signal 弱 → 没有 fine-tuning 退路**（和 UCE 一样）

### 7.2 如果接受 UCE

立即可启动，无需环境准备和基因映射，4-layer 半天就能验证 aging signal。

### 7.3 折中方案（推荐）

**UCE 为主，TranscriptFormer 为辅**：
1. 先用 UCE 4-layer 快速验证 aging signal（半天）
2. 如果信号强，继续用 UCE 33-layer 做完整分析
3. **同时**，在后台准备 TranscriptFormer 环境（新建 conda、下载权重、基因映射）
4. 如果 UCE 结果好，可以用 TranscriptFormer 的 CGE 做补充 interpretation
5. 如果 UCE 结果不理想，切换到 TranscriptFormer（此时环境已准备好）

这样既不耽误时间，又保留了 TranscriptFormer 的选项。

---

## 8. 决策点

请确认：

- [ ] **A. 坚持 TranscriptFormer** → 我立即开始新建 conda 环境、下载数据、基因映射
- [ ] **B. 改用 UCE** → 立即启动 UCE 4-layer 快速验证
- [ ] **C. 折中方案**（UCE 为主 + 后台准备 TranscriptFormer）→ 推荐

**你的选择？**
