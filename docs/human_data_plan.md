# 人类数据路线实施计划

> **日期**: 2026-06-09
> **目标**: 全部使用人类数据，驱动 Geneformer + scGPT + tAge 整合分析

---

## 一、调研结论：人类衰老 scRNA-seq 数据可用性

### 候选数据集对比

| 数据集 | 组织 | 细胞数 | 供体数 | 年龄范围 | 获取方式 | 推荐指数 |
|--------|------|--------|--------|----------|----------|----------|
| **AIDA v1** | PBMC (血液) | **1,058,909** | **508** | **19-75岁** | cellxgene | ⭐⭐⭐⭐⭐ |
| **AIDA v2** | PBMC (血液) | **1,265,624** | **~600** | **19-75岁** | cellxgene | ⭐⭐⭐⭐⭐ |
| **Kedlian Muscle** | 骨骼肌 | 183,161 | 17 | 15-75岁 | cellxgene | ⭐⭐⭐⭐ |
| **HLCA (full)** | 肺 | 2,282,447 | 多人 | 部分有年龄组 | cellxgene | ⭐⭐⭐ |
| **Tabula Sapiens** | 多组织 | 483,152 | 15 | 有限年龄信息 | cellxgene | ⭐⭐⭐ |
| **Sikkema Lung** | 肺 | 333,468 | 66 | 15-76岁 | cellxgene | ⭐⭐⭐ |

### 关键发现

1. **AIDA 是唯一兼具"大样本量 + 系统年龄标注 + 健康人"的人类单细胞数据集**
   - 508 个健康供体，年龄连续分布 19-75 岁
   - Zakar-Polyák et al. (2024, *Commun Biol*) 已证明可构建 33 个细胞类型特异性衰老时钟
   - MAE = 5.97 岁 (pseudobulk 水平)

2. **AIDA 的局限：只有 PBMC** — 缺少实体组织
   - 解决方案：用 **Kedlian Muscle Atlas** 作为跨组织验证

3. **Geneformer + scGPT 与人类数据完美兼容**
   - Geneformer 词汇表 = 人类 ENSG ID ✅
   - scGPT 词汇表 = 人类 gene symbol ✅
   - 无需任何物种映射

4. **tAge 模型与人类数据的衔接**
   - tAge 使用 **Entrez Gene ID** (10,487 features)
   - 人类 scRNA-seq 数据通常用 **gene symbol** 或 **Ensembl ID**
   - 需要一层 **Ensembl → Entrez** 映射（标准 bioinformatics 流程，mygene.info / biomart）
   - 预计映射率 ~85-90%

---

## 二、推荐实施方案

### 核心策略：AIDA 为主 + Kedlian 为辅

```
Phase 1 (第 1-2 周): 数据获取与验证
  ├── 下载 AIDA v1 (1M PBMCs)
  ├── 下载 Kedlian Muscle Atlas (183K cells)
  ├── 基因 ID 映射: ENSG / symbol → Entrez (用于 tAge)
  ├── tAge 验证: AIDA pseudobulk 预测年龄 vs 真实年龄是否相关？
  └── 输出: 确认 tAge 在人类数据上有效

Phase 2 (第 2-3 周): Geneformer 微调
  ├── AIDA 数据 tokenization (ENSG ID)
  ├── 微调年龄分类器 (young / middle / old)
  ├── In Silico Perturbation: 遍历 ~20K 基因
  └── 输出: 50-100 个虚拟抗衰老靶点 (人类基因)

Phase 3 (第 3-4 周): scGPT 生成式验证
  ├── scGPT 微调: 年龄回归 + 扰动预测
  ├── 对 Geneformer top 靶点生成预测转录组
  ├── tAge 评分 → 量化 ΔtAge
  └── 交叉验证: Geneformer 方向 vs scGPT 定量

Phase 4 (第 4 周): 跨组织验证 (Kedlian Muscle)
  ├── 在肌肉数据上重复 Geneformer ISP
  ├── 比较: 血液 vs 肌肉中 top 抗衰老靶点是否一致？
  └── 输出: 组织共享 vs 组织特异性抗衰老靶点
```

---

## 三、关键操作细节

### 3.1 基因 ID 映射管线

```python
# 输入: 人类 scRNA-seq (gene symbol 或 ENSG)
# 输出: Entrez Gene ID (用于 tAge 评分)

# 方法: mygene.info API (批量查询)
# ENSG → Entrez: ~90% 映射率
# symbol → Entrez: ~85% 映射率

# 预期损耗: ~10-15% 基因无法映射
# 应对: tAge 模型对这些缺失基因用均值 imputation (与现有 pipeline 一致)
```

### 3.2 tAge 在人类数据上的验证逻辑

```
AIDA 单细胞数据
    ↓ pseudobulk (按供体聚合)
供体级表达谱 (508 个样本)
    ↓ 基因映射 (ENSG → Entrez)
    ↓ tAge 预测
预测 biological age (508 个值)
    ↓ 与真实 chronological age 比较
验证指标:
  - Pearson r (预期 > 0.6)
  - MAE (预期 ~10-15 年)
```

如果 tAge 预测与真实年龄显著相关 → 证明模型跨物种有效，pipeline 可信。
如果不相关 → 需要分析原因（批次效应、组织差异、基因映射问题等）。

### 3.3 In Silico Perturbation 靶点筛选标准

```
Geneformer ISP 输出:
  - 每个基因的"年轻化分数" (rejuvenation score)

筛选标准 (三层过滤):
  Layer 1: rejuvenation score > top 10% (方向性)
  Layer 2: 该基因在 tAge 模型中有非零权重 (可评分)
  Layer 3: scGPT 预测该基因扰动后 ΔtAge < -5 年 (定量验证)

输出: 20-50 个高置信度人类抗衰老靶点
```

---

## 四、风险与应对

| 风险 | 概率 | 影响 | 应对 |
|------|------|------|------|
| AIDA 数据下载慢 (1M cells) | 中 | 中 | 已确认 cellxgene 可下载；如太慢，可先 subsample 100K cells 跑通 pipeline |
| tAge 在人类数据上不工作 | 中 | **高** | 这是核心验证步骤。如不相关，分析原因并记录；可能改用 AIDA 自身训练年龄模型替代 tAge |
| Geneformer 微调爆显存 (24GB) | 中 | 中 | 冻结前 6 层，只微调上层；或用 fp16；或换 104M 模型 |
| 基因映射损耗过高 (>20%) | 低 | 中 | 检查映射工具；必要时手动补充关键基因的映射 |
| AIDA 只有 PBMC，缺少组织多样性 | 确定 | 中 | 用 Kedlian Muscle 补充；明确声明结论限于免疫细胞/肌肉 |

---

## 五、最终产出

1. **tAge 人类验证报告**: 证明 tAge 模型在人类 PBMC 上的预测能力
2. **Geneformer ISP 靶点列表**: 50-100 个人类虚拟抗衰老基因靶点
3. **scGPT 定量验证**: top 靶点的预测 ΔtAge
4. **跨组织一致性分析**: 血液 vs 肌肉中靶点的重叠度
5. **与已知衰老通路的交叉验证**: top 靶点是否富集于 mTOR/AMPK/sirtuin 等通路

---

## 六、需要确认的事项

1. **AIDA 作为主力数据集是否接受？** 它是 PBMC（血液），不是实体组织，但样本量最大
2. **是否加入 Kedlian Muscle 作为跨组织验证？** 增加 ~1 周工作量
3. **tAge 验证作为必做步骤还是可选？** 如果 tAge 在人类上不工作，备选方案是用 AIDA 自身数据训练年龄回归模型
