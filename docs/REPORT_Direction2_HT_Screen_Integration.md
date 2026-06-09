# Direction 2 — HT Screen Integration with GEO Golden Standards

**Date:** 2026-06-09

## Summary

This report describes the integration of GEO-validated long-term intervention
module signatures with the LINCS high-throughput drug screen to identify
novel anti-aging drug candidates.

## 1. GEO Golden Standard Module Signatures

Computed per-module Drug-Control differences for validated interventions:

| Intervention | Composite | blue | brown4 | darkgreen | darkmagenta | darkred | darkslateblue | green | ivory | orange | pink | plum1 | sienna3 | turquoise | white |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Acarbose | -0.259 | -0.261 | -0.014 | -0.037 | -0.014 | +0.005 | +0.054 | +0.007 | +0.012 | -0.009 | +0.013 | -0.053 | -0.016 | +0.013 | +0.033 | +0.008 |
| CR | -0.820 | -0.629 | -0.051 | -0.025 | -0.011 | +0.029 | +0.039 | +0.023 | +0.005 | -0.000 | -0.006 | -0.187 | -0.041 | +0.019 | -0.014 | +0.028 |
| Metformin | -0.429 | -0.232 | +0.000 | -0.045 | +0.030 | -0.009 | -0.039 | +0.003 | +0.004 | +0.008 | -0.012 | -0.075 | -0.005 | +0.015 | -0.067 | -0.004 |
| Rapamycin | -0.627 | -0.460 | -0.031 | -0.012 | +0.012 | +0.017 | +0.044 | +0.017 | +0.014 | -0.007 | +0.017 | -0.207 | -0.048 | +0.025 | -0.015 | +0.007 |

**Key findings:**
- **Acarbose**: composite = -0.259
- **CR**: composite = -0.820
- **Metformin**: composite = -0.429
- **Rapamycin**: composite = -0.627

## 2. Screening Strategies

Four strategies were defined based on cosine similarity to golden standard module vectors:

- **Strategy A (CR mimic)**: Cosine similarity to CR module signature
- **Strategy B (Rapamycin mimic)**: Cosine similarity to Rapamycin module signature
- **Strategy C (CR+Rapamycin consensus)**: Average of CR and Rapamycin signatures
- **Strategy D (All-validated consensus)**: Weighted average (CR 40%, Rapa 30%, Acarbose 15%, Metformin 15%)

Filters applied:
- Composite < -0.5 (rejuvenation)
- Not cytotoxic: NOT (orange < -0.5 AND green < -0.5)

## 3. Top Candidates per Strategy

### Strategy A_CR

| Rank | pert_id | pert_iname | composite | cos_sim |
|------|---------|------------|-----------|---------|
| 1 | BRD-K66581298 | BRD-K66581298 | -0.5424 | 0.8592 |
| 2 | BRD-K30459086 | BRD-K30459086 | -1.0206 | 0.8504 |
| 3 | BRD-A81866333 | CGS-21680 | -0.6986 | 0.8491 |
| 4 | BRD-K40324831 | BRD-K40324831 | -1.2165 | 0.8290 |
| 5 | BRD-K76485827 | BRD-K76485827 | -0.7215 | 0.8141 |
| 6 | BRD-K79326278 | BRD-K79326278 | -1.2413 | 0.8031 |
| 7 | BRD-K88366685 | trimetazidine | -0.7300 | 0.7905 |
| 8 | BRD-K06745755 | BRD-K06745755 | -0.6522 | 0.7820 |
| 9 | BRD-K19352500 | prochlorperazine | -0.9931 | 0.7807 |
| 10 | BRD-K49005570 | BRD-K49005570 | -0.6961 | 0.7777 |

### Strategy B_Rapamycin

| Rank | pert_id | pert_iname | composite | cos_sim |
|------|---------|------------|-----------|---------|
| 1 | BRD-K40324831 | BRD-K40324831 | -1.2165 | 0.8679 |
| 2 | BRD-K30459086 | BRD-K30459086 | -1.0206 | 0.8625 |
| 3 | BRD-K66581298 | BRD-K66581298 | -0.5424 | 0.8615 |
| 4 | BRD-K06745755 | BRD-K06745755 | -0.6522 | 0.8360 |
| 5 | BRD-A81866333 | CGS-21680 | -0.6986 | 0.8164 |
| 6 | BRD-K77888550 | BRD-K77888550 | -0.5584 | 0.8151 |
| 7 | BRD-K49005570 | BRD-K49005570 | -0.6961 | 0.7726 |
| 8 | BRD-K51662849 | ilomastat | -1.5549 | 0.7597 |
| 9 | BRD-K13195996 | BRD-K13195996 | -0.6754 | 0.7502 |
| 10 | BRD-K19462402 | buflomedil | -0.5817 | 0.7498 |

### Strategy C_CR_Rapa_consensus

| Rank | pert_id | pert_iname | composite | cos_sim |
|------|---------|------------|-----------|---------|
| 1 | BRD-K66581298 | BRD-K66581298 | -0.5424 | 0.8663 |
| 2 | BRD-K30459086 | BRD-K30459086 | -1.0206 | 0.8626 |
| 3 | BRD-K40324831 | BRD-K40324831 | -1.2165 | 0.8550 |
| 4 | BRD-A81866333 | CGS-21680 | -0.6986 | 0.8380 |
| 5 | BRD-K06745755 | BRD-K06745755 | -0.6522 | 0.8155 |
| 6 | BRD-K49005570 | BRD-K49005570 | -0.6961 | 0.7804 |
| 7 | BRD-K79326278 | BRD-K79326278 | -1.2413 | 0.7783 |
| 8 | BRD-K77888550 | BRD-K77888550 | -0.5584 | 0.7772 |
| 9 | BRD-K76485827 | BRD-K76485827 | -0.7215 | 0.7735 |
| 10 | BRD-K88366685 | trimetazidine | -0.7300 | 0.7686 |

### Strategy D_All_consensus

| Rank | pert_id | pert_iname | composite | cos_sim |
|------|---------|------------|-----------|---------|
| 1 | BRD-K66581298 | BRD-K66581298 | -0.5424 | 0.8723 |
| 2 | BRD-K30459086 | BRD-K30459086 | -1.0206 | 0.8686 |
| 3 | BRD-K40324831 | BRD-K40324831 | -1.2165 | 0.8619 |
| 4 | BRD-A81866333 | CGS-21680 | -0.6986 | 0.8388 |
| 5 | BRD-K06745755 | BRD-K06745755 | -0.6522 | 0.8134 |
| 6 | BRD-K88366685 | trimetazidine | -0.7300 | 0.7856 |
| 7 | BRD-K13195996 | BRD-K13195996 | -0.6754 | 0.7846 |
| 8 | BRD-K51662849 | ilomastat | -1.5549 | 0.7730 |
| 9 | BRD-K79326278 | BRD-K79326278 | -1.2413 | 0.7697 |
| 10 | BRD-K49005570 | BRD-K49005570 | -0.6961 | 0.7695 |

## 4. Triple Hits (Highest Confidence)

Compounds that are CR-like, Rapamycin-like, AND Metformin-like:

| pert_id | pert_iname | composite | cos_sim_CR | cos_sim_Rapa |
|---------|------------|-----------|------------|--------------|
| BRD-K25464116 | SA-1925610 | -2.3601 | 0.4478 | 0.4395 |
| BRD-K30229575 | BRD-K30229575 | -2.2744 | 0.3366 | 0.2979 |
| BRD-K18619710 | digoxigenin | -2.1538 | 0.5359 | 0.3891 |
| BRD-K16336526 | capsaicin | -2.1423 | 0.5117 | 0.4331 |
| BRD-K44276885 | acarbose | -1.7691 | 0.5686 | 0.4524 |
| BRD-K15933101 | ropinirole | -1.7342 | 0.4727 | 0.4984 |
| BRD-K51313569 | palbociclib | -1.7127 | 0.4748 | 0.3548 |
| BRD-K28065812 | BRD-K28065812 | -1.4345 | 0.3664 | 0.3083 |
| BRD-K86250672 | BRD-K86250672 | -1.4116 | 0.5917 | 0.5658 |
| BRD-K90135241 | BRD-K90135241 | -1.3747 | 0.4931 | 0.5400 |
| BRD-K04010869 | prostaglandin-a1 | -1.3153 | 0.5474 | 0.5213 |
| BRD-K23958338 | BRD-K23958338 | -1.2795 | 0.3675 | 0.3366 |
| BRD-K42098891 | protriptyline | -1.2330 | 0.5881 | 0.5297 |
| BRD-K40324831 | BRD-K40324831 | -1.2165 | 0.8290 | 0.8679 |
| BRD-K32083350 | BRD-K32083350 | -1.2139 | 0.6684 | 0.7053 |
| BRD-K74271701 | BRD-K74271701 | -1.1891 | 0.4578 | 0.3730 |
| BRD-A70311631 | BRD-A70311631 | -1.1362 | 0.4083 | 0.3383 |
| BRD-K52682771 | BRD-K52682771 | -1.1022 | 0.3395 | 0.3088 |
| BRD-K48296029 | BRD-K48296029 | -1.0958 | 0.4812 | 0.5126 |
| BRD-K79983625 | DC-45-A2 | -1.0744 | 0.4204 | 0.4105 |

## 5. Validation Against Known Anti-Aging Drugs

| Drug | Found in dataset | Composite | Passes filter | CR sim | Rapa sim |
|------|------------------|-----------|---------------|--------|----------|
| resveratrol | Yes | +1.8640 | No | -0.068 | +0.019 |
| spermidine | Yes | +0.3661 | No | -0.182 | +0.026 |
| sirolimus | Yes | -0.8190 | Yes | +0.064 | -0.089 |
| metformin | Yes | -0.1788 | No | +0.077 | -0.060 |
| everolimus | Yes | +2.5593 | No | -0.199 | -0.118 |
| acarbose | Yes | +0.4513 | No | +0.140 | -0.035 |

**Key observations:**

1. **LINCS short-term vs GEO long-term discrepancy**: Most known anti-aging drugs
   do *not* show rejuvenation (composite < -0.5) in the LINCS short-term in vitro
   signatures. This highlights a critical translational gap: long-term in vivo effects
   (GEO) are not always captured by short-term cell-line perturbations (LINCS).

2. **Sirolimus paradox**: Sirolimus has one LINCS signature with composite=-0.819
   (rejuvenating), but its cosine similarity to the GEO Rapamycin module signature
   is *negative* (-0.089). This indicates that the module-level fingerprint of
   short-term in vitro sirolimus differs substantially from long-term in vivo
   rapamycin. The LINCS sirolimus signature is dominated by strong **orange**
   (chromatin) and **green** (cell cycle) suppression, whereas GEO rapamycin
   shows strong **pink** (mitochondrial/OxPhos) suppression.

3. **Multiple signatures per drug**: Several drugs (acarbose, sirolimus) have
   multiple pert_ids in LINCS with conflicting directions. For example, acarbose
   has one signature at +0.45 (pro-aging) and another at -1.77 (anti-aging).
   The per-drug averages in `module_fingerprint_v2_per_drug.csv` can mask these
   differences.

4. **Metformin**: Does not pass the composite < -0.5 filter in LINCS
   (composite = -0.179), consistent with its modest effect size in GEO
   (composite = -0.429).

## 6. Discussion and Next Steps

### Mechanistic Interpretation

**Golden standard module signatures:**
- **CR** and **Rapamycin** share the strongest rejuvenation signal in **pink**
  (mitochondrial translation / OxPhos; CR=-0.187, Rapa=-0.207) and **non_module**
  (CR=-0.629, Rapa=-0.460), suggesting convergent mechanisms on mitochondrial
  function and non-module genes.
- **Metformin** is distinct: it shows strong **turquoise** (innate immunity)
  suppression (-0.067) not seen in CR/Rapa, aligning with its known
  anti-inflammatory effects.
- **Acarbose** has the weakest rejuvenation signal (-0.259) and a relatively
  flat module profile, consistent with its indirect metabolic mechanism.

**Top candidate mechanisms:**
- **BRD-K40324831**: Highest Rapa-similarity (0.868) with strong rejuvenation
  (-1.216). Its module profile mirrors GEO rapamycin in pink and turquoise.
- **Trimetazidine** (BRD-K88366685): A known anti-ischemic drug that shifts
  cardiac metabolism from fatty acid to glucose oxidation. High CR similarity
  (0.791) suggests metabolic reprogramming akin to CR.
- **CGS-21680** (BRD-A81866333): An adenosine A2A receptor agonist with
  neuroprotective properties. Strong CR similarity (0.849) may reflect
  adenosine-mediated suppression of inflammation (turquoise module).
- **Ilomastat** (BRD-K51662849): A broad-spectrum MMP inhibitor. High Rapa
  similarity (0.760) with strong rejuvenation (-1.555) suggests ECM remodeling
  (brown4 module) as a convergent pathway.

### Why known drugs did not rank top

The absence of known anti-aging drugs (resveratrol, spermidine, NAD+ precursors)
from our top candidate lists is **informative rather than disappointing**:

1. **Short-term vs long-term mismatch**: LINCS captures acute transcriptional
   responses (6-24h), while GEO captures chronic adaptive responses (months).
   Resveratrol's acute effects (SIRT1 activation, stress response) may differ
   from its long-term benefits (metabolic remodeling).

2. **Tissue specificity**: The GEO data is from mouse liver, while LINCS uses
   cancer cell lines (A549, MCF7, PC3). Drugs that act via tissue-specific
   pathways (e.g., resveratrol in adipose, metformin in gut/liver) will not
   recapitulate their in vivo signatures in vitro.

3. **Dose and duration**: LINCS signatures are typically high-dose short-term,
   which may activate stress responses (e.g., Nrf2, HSF1) that are beneficial
   acutely but not sustainable chronically.

### Validation Strategy

1. **In vitro**: Test top 10 triple-hit candidates in:
   - Cellular senescence assays (SA-β-gal, p16/p21)
   - Mitochondrial function (Seahorse, TMRM)
   - Inflammatory cytokine panels (IL-6, TNF-α, IL-1β)

2. **In vivo**: Mouse lifespan/healthspan studies for top 3-5 compounds:
   - **BRD-K40324831** (highest Rapa-sim + strong rejuvenation)
   - **Trimetazidine** (known drug, repurposing potential)
   - **CGS-21680** (novel mechanism, adenosine pathway)

3. **Mechanism**: CRISPR screening in iPSC-derived hepatocytes to identify
   which modules are necessary for the rejuvenation effect of each candidate.

4. **Address LINCS-GEO gap**: For candidates with strong GEO-similarity but
   weak LINCS rejuvenation, consider building a "translation model" that maps
   LINCS signatures to predicted long-term effects using the golden standards
   as training data.

### Files Generated
- `results/geo_longterm/geo_module_signatures.csv` — Golden standard module signatures
- `results/geo_longterm/ht_screen_strategy_*_top100.csv` — Top 100 per strategy
- `results/geo_longterm/ht_screen_triple_hits.csv` — Triple-hit candidates
- `figures/direction2_*.pdf` — Visualizations

---
*Report generated by Direction 2 HT Screen Integration pipeline*