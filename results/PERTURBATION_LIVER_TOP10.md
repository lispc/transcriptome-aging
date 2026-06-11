# Liver Top 10 Pro-Aging Genes — In-Silico Perturbation (4-Layer UCE)

## Overview

In-silico KO perturbation on the top 10 pro-aging genes identified from Liver aging axis
correlation (4-layer UCE embeddings, n=7,294 cells).

**Method**: Zero-out expression of target gene in all cells, re-run UCE 4-layer inference,
compute aging score shift along pre-computed aging axis.

**Expected**: True pro-aging genes should show negative mean_delta (aging score decreases
after knock-out = rejuvenation effect).

---

## Results

| Rank | Gene | mean_delta | young_delta | old_delta | Effect |
|------|------|-----------:|------------:|----------:|--------|
| 1 | **Cd52** | **-0.000667** | -0.000787 | -0.000535 | Pro-aging |
| 2 | **H2-D1** | **-0.000581** | -0.001102 | -0.000010 | Pro-aging |
| 3 | **Cd53** | **-0.000572** | -0.001818 | +0.000797 | Pro-aging |
| 4 | **Lgmn** | **-0.000486** | -0.000622 | -0.000337 | Pro-aging |
| 5 | **Tyrobp** | **-0.000303** | -0.000713 | +0.000147 | Pro-aging |
| 6 | **Ptpn18** | **-0.000201** | -0.000708 | +0.000356 | Pro-aging |
| 7 | **Arpc1b** | **-0.000019** | -0.000727 | +0.000758 | Weak pro-aging |
| 8 | Unc93b1 | +0.000283 | +0.000300 | +0.000265 | **Anti-aging** |
| 9 | Laptm5 | +0.000414 | +0.000165 | +0.000686 | **Anti-aging** |
| 10 | C1qb | +0.000670 | -0.000384 | +0.001828 | **Anti-aging** |

### Summary
- **Pro-aging (validated)**: 7/10 genes show negative delta
- **Anti-aging (surprising)**: 3/10 genes show positive delta
- **Strongest rejuvenation**: Cd52, H2-D1, Cd53
- **Most ambiguous**: Arpc1b (delta ≈ 0)

---

## Key Insights

### 1. Correlation ≠ Causation
3 genes (Unc93b1, Laptm5, C1qb) that were highly **correlated** with aging score
showed **anti-aging** effects upon perturbation. This demonstrates that:
- High expression in old cells does not imply causal role in aging
- Some genes may be **compensatory** (upregulated to counteract aging)
- Aging axis correlation alone cannot distinguish driver vs. passenger genes

### 2. Effect Sizes Are Small
All |delta| < 0.001, reflecting that:
- Aging is a polygenic trait; single-gene perturbation has limited effect
- UCE 4-layer embeddings may be relatively robust to single-gene changes
- Batch effects and biological noise may dilute signal

### 3. Young vs. Old Cell Responses Differ
- **Cd53**: Young cells show much stronger rejuvenation (-0.0018) than old cells (+0.0008)
  → Suggests Cd53 plays a larger role in young liver aging dynamics
- **C1qb**: Old cells show aging acceleration (+0.0018) while young cells show slight rejuvenation
  → May have context-dependent (age-specific) functions

---

## Technical Notes

- **Model**: UCE 4-layer (3.4GB, 1280-dim embeddings)
- **Runtime**: ~8 minutes for 10 genes (4 GPUs parallel, ~2.5 min/gene)
- **Batch size**: 100 (inference), preprocessing on CPU
- **Work dirs**: `results/perturbation_liver_4l/work_gpu{0-3}/`
- **Output**: `results/perturbation_liver_4l/liver_top10_perturbation.csv`

### Bug Fixes Applied to `perturbation.py`
1. `script_dir` changed from relative `Path("data/uce-repo").resolve()` to absolute path via `Path(__file__)`
2. Added missing `aging_score` computation in `main()` (original adata lacked this column)
3. Added trailing `os.sep` to `--dir` argument (UCE concatenates dir + filename without separator)

---

## Next Steps

1. **Run 33-layer perturbation** for comparison (may show larger effect sizes)
2. **Lung tissue perturbation** on top 10 Lung pro-aging genes
3. **Multi-gene combinatorial KO** (e.g., Cd52 + H2-D1 + Cd53) to test synergy
4. **Compare with known aging signatures** from literature
