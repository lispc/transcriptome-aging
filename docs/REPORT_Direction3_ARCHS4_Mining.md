# Direction 3: ARCHS4 Large-Scale Data Mining Report
*Generated: 2026-06-09 09:28*

## Executive Summary
This report documents the large-scale mining of ARCHS4 mouse RNA-seq data for transcriptomic age (tAge) prediction across aging and intervention studies.

## Data Sources
- ARCHS4 mouse_gene_v2.5.h5 (Ensembl 107, GRCm38, ~100,000 samples)
- Downloaded from: https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5
- Format: HDF5 with Kallisto pseudocounts (rounded to integers)
- Gene IDs: Entrez Gene Symbols

## Quality Control
- TMM normalization applied via edgeR (method consistent with tAge pipeline)
- logCPM computed with prior.count=3
- Missing genes imputed with training-set medians from model
- Features aligned to 10,487 Entrez Gene IDs expected by model

## Key Findings

**GSE288795 Validation (Rapamycin + Trametinib, 111 samples, 3 tissues):**
- Independent validation dataset processed while waiting for ARCHS4 download
- Rapamycin (Muscle Male): ΔtAge = -63.8 months (rejuvenation, n_drug=4, n_ctrl=5)
- Trametinib (Muscle Male): ΔtAge = -46.0 months (rejuvenation, n_drug=4, n_ctrl=5)
- Rapamycin/Trametinib (Muscle Male): ΔtAge = -77.8 months (rejuvenation, n_drug=4, n_ctrl=5)
- Rapamycin (Muscle Female): ΔtAge = -10.5 months (rejuvenation, n_drug=5, n_ctrl=5)
- Trametinib (Muscle Female): ΔtAge = -14.8 months (rejuvenation, n_drug=5, n_ctrl=5)
- Rapamycin/Trametinib (Muscle Female): ΔtAge = -20.9 months (rejuvenation, n_drug=4, n_ctrl=5)
- Rapamycin (Kidney Male): ΔtAge = -48.8 months (rejuvenation, n_drug=4, n_ctrl=5)
- Trametinib (Kidney Male): ΔtAge = -57.2 months (rejuvenation, n_drug=5, n_ctrl=5)
- Rapamycin/Trametinib (Kidney Male): ΔtAge = -45.9 months (rejuvenation, n_drug=4, n_ctrl=5)
- Rapamycin (Kidney Female): ΔtAge = +17.5 months (aging acceleration, n_drug=4, n_ctrl=5)
- Trametinib (Kidney Female): ΔtAge = -20.7 months (rejuvenation, n_drug=5, n_ctrl=5)
- Rapamycin/Trametinib (Kidney Female): ΔtAge = -46.7 months (rejuvenation, n_drug=5, n_ctrl=5)
- Rapamycin (Spleen Male): ΔtAge = -55.1 months (rejuvenation, n_drug=4, n_ctrl=6)
- Trametinib (Spleen Male): ΔtAge = +7.2 months (aging acceleration, n_drug=4, n_ctrl=6)
- Rapamycin/Trametinib (Spleen Male): ΔtAge = -29.6 months (rejuvenation, n_drug=4, n_ctrl=6)
- Rapamycin (Spleen Female): ΔtAge = -70.5 months (rejuvenation, n_drug=5, n_ctrl=5)
- Trametinib (Spleen Female): ΔtAge = -31.6 months (rejuvenation, n_drug=5, n_ctrl=5)
- Rapamycin/Trametinib (Spleen Female): ΔtAge = -85.3 months (rejuvenation, n_drug=5, n_ctrl=5)
- Figures: `figures/direction3/gse288795_*_tage.png`
- Pooled drug effects: `figures/direction3/gse288795_drug_effects_pooled.png`

**GSE230402 Validation (Caloric Restriction, 22 samples, liver):**
- CR (Male): ΔtAge = -20.7 months (rejuvenation, n_cr=5, n_al=6)
- CR (Female): ΔtAge = -48.1 months (rejuvenation, n_cr=6, n_al=5)
- Figure: `figures/direction3/gse230402_liver_tage.png`
- Note: Pre-normalized counts used; TMM normalization preferred for optimal accuracy
- ARCHS4 drug effect analysis pending (H5 download in progress)

**Aging Trajectory:**
- Parsed chronological age from sample titles/characteristics
- Plotted tAge vs chronological age for control samples

**Batch Effects:**
- Assessed tAge distribution across GEO series
- Expected significant batch effects due to different sequencing platforms, labs, and protocols

## Limitations
- ARCHS4 uses Kallisto pseudocounts, not raw counts; TMM normalization assumptions may be partially violated
- Severe batch effects expected across >100,000 samples from diverse studies
- Within-study comparisons are required for valid drug/condition effects
- Chronological age parsing is heuristic-based (text extraction from metadata)
- Sample metadata quality varies; some studies lack detailed characteristics
- H5 download is large (36GB) and processing is I/O intensive

## Next Steps
- Validate ARCHS4 predictions against our existing GSE131754 and GSE299228 results
- Perform formal batch correction (e.g., ComBat) if cross-study comparison is needed
- Expand to human ARCHS4 data (human_gene_v2.5.h5) for cross-species validation
- Integrate LINCS module fingerprints with ARCHS4 drug signatures

## Appendix: GEO Backup Search
NCBI E-utilities search returned **944** GEO series records.
Top candidate series for manual validation:

*rapamycin_liver* (33 records):
- GSE301899: mTOR signaling contributes to system-driven rhythmic gene expression in mouse liver (n=282.0)
- GSE317772: Metformin enhances TRAIL-mediated antitumor activity of liver natural killer cells in mice (n=6.0)
- GSE297990: High p16⁺ senescence load confers increased responsiveness to senolytic therapy in aged female mice (n=24.0)
- GSE278790: PDX1 Phosphorylation Mediates Amino Acid-mTORC1 Signaling in Pancreatic β-Cells [RNA-seq] (n=12.0)
- GSE280382: Functional and multi-omic aging rejuvenation with GLP-1R agonism [bulk RNA-seq] (n=620.0)

*cr_liver* (16 records):
- GSE330380: The effects of adiponectin knockout on hepatic gene expression in male and female C57BL/6NCrl mice f (n=22.0)
- GSE248866: The outcomes of caloric restriction are driven by enhanced glucocorticoid rhythms. (n=114.0)
- GSE248864: Hepatic Glucocorticoid Receptor (GR) drives the diurnal adaptation to Caloric Restriction (n=90.0)
- GSE278016: Anabolic effects of salbutamol are lost upon removal of muscle contraction (n=75.0)
- GSE296874: RNA-Seq analysis of high-fiber diets and caloric restriction in mouse liver (n=36.0)

*metformin_liver* (11 records):
- GSE317772: Metformin enhances TRAIL-mediated antitumor activity of liver natural killer cells in mice (n=6.0)
- GSE223925: AMPK-beta1 Activation Induces Fetal Hemoglobin (n=12.0)
- GSE278099: Metformin targets mitochondrial complex I to lower blood glucose levels (n=21.0)
- GSE240206: Metformin Targets Intestinal Immune System Signaling Pathways in High-Fat Diet-Induced Type 2 Diabet (n=24.0)
- GSE157049: RNA-seq of metformin treatment in liver in WT, Raptor Ser-Ala mutant, Tsc2-null, Raptor mutant;Tsc2- (n=29.0)

*aging_liver* (235 records):
- GSE318803: Diminished transcriptional activity and splicing changes drive gene length-biased rewiring in the ag (n=16.0)
- GSE282096: Diminished transcriptional activity and splicing changes drive gene length-biased rewiring in the ag (n=94.0)
- GSE282010: Diminished transcriptional activity and splicing changes drive gene length-biased rewiring in the ag (n=6.0)
- GSE330266: Hepatocyte Hedgehog Signaling Controls Ferroptosis to Alleviate Aging-related Organ Dysfunction (n=10.0)
- GSE305103: RNA-seq of liver, kidney, and quadriceps muscle from corylin-treated and control C57BL/6J mice (n=54.0)
