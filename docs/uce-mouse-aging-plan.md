# UCE × Tabula Muris Senis: 纯小鼠衰老分析方案

> 方案版本: v1.0 | 目标: 用户确认后启动

---

## 1. 为什么选 UCE

### 1.1 UCE 是什么

**UCE (Universal Cell Embedding)** 是 Stanford 2023 年提出的单细胞 foundation model。

| 属性 | 详情 |
|------|------|
| 架构 | 33-layer Transformer (~650M 参数) |
| 训练数据 | 3600 万细胞，8 个物种 |
| **小鼠** | ✅ 训练数据明确包含小鼠 |
| 基因表示 | ESM2 蛋白质语言模型编码氨基酸序列 |
| 物种支持 | **任何有蛋白质序列的物种**，无需同源映射 |
| 代码 | https://github.com/snap-stanford/UCE |
| 推理脚本 | `eval_single_anndata.py` |

### 1.2 与其他模型的对比

| 模型 | 小鼠支持 | 需要基因映射 | 支持 fine-tune | 显存需求 |
|------|---------|-------------|---------------|---------|
| Geneformer | ❌ 人类 ENSG 专属 | 必须映射 | ✅ | 24GB 可跑 |
| scGPT | ⚠️ vocab 小写匹配 | 不需要 | ✅ | 24GB 可跑 |
| **UCE** | **✅ 原生支持** | **不需要** | **❌ 仅 zero-shot** | **24GB 可跑 (batch_size 调小)** |

### 1.3 关键限制（必须知情）

1. **不支持 fine-tuning**：UCE 设计为 as-is 使用。论文明确说 "we have designed UCE to be used as-is, without any model fine-tuning"。
   - 这意味着 aging classifier 是在 UCE 的**固定 embedding** 上训练的，不是 end-to-end。
   - 如果 UCE 的 zero-shot embedding 没有 capture aging signal，我们无法通过微调 UCE 来改善。

2. **GitHub issue #57**：修改单个细胞的基因表达后，**其他细胞的 embedding 也可能改变**。这是因为 batch processing 和数据归一化中的信息共享。
   - 影响：perturbation 分析（ISP）的结果解读需要谨慎。
   - 缓解：`--batch_size 1 --no-filter`。

3. **输出是 cell-level**：UCE 输出 1280-dim cell embedding，不是 gene-level 表示。perturbation 需要修改输入重新 embedding。

---

## 2. 数据：Tabula Muris Senis

| 参数 | 详情 |
|------|------|
| 来源 | Figshare (h5ad, 3.52 GB) |
| 细胞数 | 50 万+ |
| 组织 | 18+ (脾、骨髓、肝、肺、肾、心、肌肉、脑、脂肪等) |
| 年龄 | **1, 3, 18, 21, 24, 30 月** — 6 个连续时间点 |
| 技术 | 10x Genomics droplet |
| 基因名 | 小鼠 gene symbol (如 `Actb`, `Cd3d`) |

### 2.1 数据获取

```bash
# Figshare 直接下载
wget https://figshare.com/ndownloader/files/44704373 -O tms_droplet.h5ad

# 或使用 Bioconductor R 包
# TabulaMurisSenisDroplet(tissues="All")
```

### 2.2 预处理策略

1. **质控**：每个组织单独过滤（避免组织间差异淹没 aging signal）
2. **目标组织**：建议先聚焦 **免疫相关组织**（脾、骨髓、胸腺、肺、肝），因为免疫细胞 aging 信号最强
3. **最小细胞数**：每个 age × tissue × cell_type 组合至少 50 个细胞

---

## 3. 完整 Pipeline

### Phase 1: UCE Zero-Shot Embedding (1-2 天)

```bash
# 1. 克隆 UCE
git clone https://github.com/snap-stanford/UCE.git
cd UCE
pip install -r requirements.txt

# 2. 运行 embedding（以单个组织为例）
python eval_single_anndata.py \
  --adata_path data/tms_spleen.h5ad \
  --dir results/uce_embeddings/ \
  --species mouse \
  --nlayers 33 \
  --model_loc data/33l_8ep_1024t_1280.torch \
  --batch_size 10
```

**模型选择**：
- **33-layer** (主模型)：650M 参数，精度更高，但 slower
- **4-layer** (轻量版)：更快，默认选项，适合快速迭代
- 建议：先用 4-layer 跑通，再用 33-layer 做最终分析

**显存估算**：
- 33-layer + batch_size 10：~20GB（我们的 RTX 3090 24GB 刚好够）
- 4-layer + batch_size 50：~10GB

**输出**：每个细胞的 `obsm["X_uce"]` — 1280-dim 向量

---

### Phase 2: Aging Classifier (0.5-1 天)

在固定的 UCE embedding 上训练 aging classifier。

#### 2.1 任务定义

**任务 A: 二分类**（简单，baseline）
- Young: 1m + 3m
- Old: 18m + 21m + 24m + 30m

**任务 B: 连续回归**（更有科学价值）
- 预测精确月龄（1, 3, 18, 21, 24, 30）

#### 2.2 模型选择

| 模型 | 优势 | 劣势 |
|------|------|------|
| Logistic Regression | 可解释性强，快 | 线性，可能 underfit |
| ElasticNet | 特征选择，鲁棒 | 线性 |
| **MLP (PyTorch)** | **非线性，可 fine-tune** | 需要调参 |
| XGBoost | 非线性，特征重要性 | 可能 overfit |

**推荐**：用多种模型尝试，以 ElasticNet 为主（和 Zakar-Polyák 2024 直接对标），MLP 为辅。

#### 2.3 交叉验证策略

```
CV-1: Leave-one-tissue-out
  → 在脾+骨髓+胸腺上训练，在肺上测试
  → 验证跨组织泛化性

CV-2: Leave-one-age-out
  → 在 1,3,18,21,24m 上训练，在 30m 上测试
  → 验证对未见年龄的预测能力

CV-3: Cell-type-stratified
  → 确保每种 cell type 在 train/test 中都有代表
```

#### 2.4 评估指标

- 二分类：Accuracy, AUC-ROC, F1
- 回归：MAE (月), Pearson r, Spearman ρ
- 跨组织：组织-specific 的 AUC 差异

---

### Phase 3: Aging Axis 分析 (0.5 天)

1. 计算 young_mean = mean(UCE_embedding[age ≤ 3m])
2. 计算 old_mean = mean(UCE_embedding[age ≥ 18m])
3. Aging axis = old_mean - young_mean
4. 对每个细胞：aging_score = dot(UCE_embedding, aging_axis)
5. 分析 aging_score 在以下维度的分布：
   - 组织间差异（哪个组织 aging 最快？）
   - Cell type 间差异（哪种细胞 aging 最快？）
   - 性别差异（如果有数据）

---

### Phase 4: Interpretation & Perturbation (1 天)

#### 4.1 特征重要性（解释 UCE embedding 中的 aging signal）

由于 UCE 是 black box，我们用 post-hoc 解释：

```python
# SHAP 分析每个 UCE dimension 对 aging 预测的贡献
import shap
explainer = shap.TreeExplainer(xgb_model)
shap_values = explainer.shap_values(X_uce)
```

然后反向映射：哪些基因对这些重要的 UCE dimensions 贡献最大？

#### 4.2 Perturbation 分析（探索性，受 issue #57 限制）

对于 top candidate genes：

```python
# 1. 获取 baseline embedding
baseline_emb = uce_embed(adata)

# 2. 删除基因 X 的表达
adata_pert = adata.copy()
adata_pert[:, "GeneX"].X = 0

# 3. 重新 embedding
perturbed_emb = uce_embed(adata_pert, batch_size=1, no_filter=True)

# 4. 比较 aging score 变化
delta = aging_classifier(perturbed_emb) - aging_classifier(baseline_emb)
```

**重要警告**：由于 issue #57，perturbed_emb 可能不是完全独立的。结果需保守解读。

---

### Phase 5: 与已知 Aging Signatures 比较 (0.5 天)

1. **Tabula Muris Senis 官方 aging signatures**
   - MSigDB 中有 100+ 个 `TABULA_MURIS_SENIS_*_AGEING` gene sets
   - GitHub 上有完整的 DGE 结果

2. **我们的 aging classifier 预测 vs 官方 signatures**
   - 重叠分析（Fisher exact test）
   - 哪些组织/cell type 的一致性最高？

3. **通路富集**
   - 对 aging-associated UCE dimensions 对应的基因做 GO/KEGG
   - 特别关注：核糖体、线粒体、免疫炎症通路

---

## 4. 技术风险与缓解

| 风险 | 可能性 | 影响 | 缓解 |
|------|--------|------|------|
| UCE 安装失败 | 低 | 高 | 有活跃 GitHub repo，有 requirements.txt |
| 33-layer 显存 OOM | 中 | 中 | Fallback 到 4-layer；或降低 batch_size |
| Protein embedding 下载慢/失败 | 中 | 中 | UCE 会自动下载；可手动缓存 |
| **Aging signal 太弱** | **中** | **高** | **先用 4-layer 快速验证；如果信号弱，考虑聚焦特定组织/cell type** |
| Issue #57 影响 perturbation | 中 | 中 | batch_size=1；结果保守解读 |
| 数据下载慢 | 中 | 低 | Figshare 3.5GB，可用 wget/curl 断点续传 |

---

## 5. 时间线

| 阶段 | 预计时间 | 产出 |
|------|---------|------|
| **Day 1** 数据下载 + UCE 安装 + 4-layer 快速验证 | 1 天 | 确认 aging signal 存在 |
| **Day 2** 33-layer embedding + 数据整合 | 1 天 | 所有组织的 UCE embeddings |
| **Day 3** Aging classifier 训练 + 评估 | 1 天 | 跨组织/跨年龄的 aging model |
| **Day 4** Aging axis + Interpretation + Perturbation | 1 天 | Aging score 分布、特征重要性、perturbation 结果 |
| **Day 5** 与已知 signatures 比较 + 文档 | 0.5-1 天 | 最终报告、图表、结论 |

**总计: 3.5-5 天**

---

## 6. 预期产出

### 6.1 科学产出

1. **第一个 UCE-based 小鼠跨组织衰老时钟**
2. **组织 aging rate 排名**：哪些组织的细胞在 UCE space 中 aging 最快？
3. **Cell-type-specific aging trajectories**
4. **与已知 aging signatures 的验证/对比**

### 6.2 方法学产出

1. **UCE 在 aging 任务上的性能基准**：zero-shot embedding + 简单 classifier 能做到什么程度？
2. **4-layer vs 33-layer 的 trade-off 分析**
3. **UCE 的 interpretability 方法论**：如何解释 1280-dim black box 中的 aging signal？

### 6.3 潜在发现

- 如果在特定 cell type（如 T cell、巨噬细胞）中发现 strong aging signal，可以深挖其 marker genes
- 如果核糖体/线粒体通路再次出现（像我们人类 PBMC 结果一样），这将是跨物种的验证

---

## 7. 关键决策点

在启动前，请确认以下决策：

### 7.1 模型选择

- [ ] **先用 4-layer 做快速验证**（半天出结果，确认 aging signal 存在）
- [ ] **再用 33-layer 做最终分析**（精度更高，但慢 3-5 倍）

### 7.2 组织选择

- [ ] **全组织**（18+ 组织，计算量大但全面）
- [ ] **免疫聚焦**（脾、骨髓、胸腺、肺、肝，计算量适中，aging signal 可能更强）

### 7.3 Aging classifier 类型

- [ ] **二分类**（young vs old，简单稳健）
- [ ] **连续回归**（预测月龄，更有科学价值但更难）
- [ ] **两者都做**

### 7.4 Perturbation 深度

- [ ] **浅尝辄止**（只做 top 10-20 genes 的 perturbation，受 issue #57 限制）
- [ ] **系统筛选**（100-200 genes，需要更多时间，需谨慎解读）

---

## 8. 备选方案

如果 UCE 路线遇到不可克服的障碍（如 aging signal 太弱、安装失败、issue #57 严重影响结果），备选方案是：

**方案 B：scGPT + 小鼠数据**
- scGPT vocab 支持小写基因名（小鼠基因转小写后直接匹配）
- 支持 fine-tuning 和 ISP
- 但预训练主要是人类数据

**方案 C：传统 ML（ElasticNet/LightGBM）**
- 直接在基因表达上训练 aging clock
- 稳健、可解释，但不是 foundation model

---

## 9. 下一步

**请确认以上方案后，我将：**

1. 立即开始下载 Tabula Muris Senis 数据
2. 克隆 UCE repo 并安装
3. 先用 4-layer 模型在 1-2 个组织上做快速验证
4. 如果验证通过，扩展到全部组织和 33-layer 模型

**你的决定是？**
