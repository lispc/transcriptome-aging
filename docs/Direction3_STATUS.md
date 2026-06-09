# Direction 3: ARCHS4 Mining — Current Status

*Last updated: 2026-06-09 09:33 UTC*

## Background Tasks Running

### 1. ARCHS4 H5 Download
- **File**: `data/archs4/mouse_gene_v2.5.h5` (38.96 GB)
- **Progress**: ~10.0 GB / 38.96 GB (**25.7%**)
- **Speed**: ~8-11 MB/s (variable)
- **ETA**: ~45-60 minutes
- **Process**: Running in `screen` session `archs4_download` (PID via screen)
- **Log**: `logs/direction3/wget_download.log`
- **Resume**: Supported (wget -c)

### 2. Auto-Run Watcher
- **Process**: `screen` session `archs4_watcher`
- **Purpose**: Polls H5 file size every 5 minutes; auto-triggers full pipeline when download completes
- **Pipeline**: `src/python/direction3/run_pipeline.sh`
- **Log**: `logs/direction3/auto_run.log`

## Completed Work

### Phase 1: API Test & Download
- ARCHS4 API endpoint (`/archs4/search/search`) returned **404** — confirmed broken
- Found correct H5 download URL: `https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5`
- Verified resume support and disk space (~256 GB available)
- Download started and running stably

### Phase 2: GEO Backup Search
- NCBI E-utilities search completed: **944 GEO series** records found
- Queries: rapamycin, CR, metformin, aging, dexamethasone, doxorubicin, young/old liver
- Results saved: `data/geo_expansion/geo_eutils_search_results.csv`
- Top candidates identified:
  - GSE288795: Rapamycin + Trametinib (111 samples)
  - GSE248866: CR liver (114 samples)
  - GSE282096: Aging liver (94 samples)

### Phase 3: Immediate Validation Datasets

#### GSE288795 (Rapamycin + Trametinib)
- **Samples**: 111 mouse samples across Muscle, Kidney, Spleen
- **Treatments**: Control, Rapamycin, Trametinib, Rapamycin/Trametinib combo
- **Result**: Consistent rejuvenation across tissues
  - Rapamycin pooled: **-38.2 months**
  - Trametinib pooled: **-27.3 months**
  - Combo pooled: **-50.7 months**
  - Strongest in Spleen Female Combo: **-85.3 months**
- **Figures**: `figures/direction3/gse288795_*`
- **Data**: `results/direction3/gse288795_predictions.csv`

#### GSE230402 (Caloric Restriction)
- **Samples**: 22 mouse liver samples (11 male, 11 female)
- **Design**: Ad libitum vs 30% CR
- **Result**: CR shows rejuvenation
  - Male: **-20.7 months**
  - Female: **-48.1 months**
- **Figures**: `figures/direction3/gse230402_liver_tage.png`
- **Data**: `results/direction3/gse230402_predictions.csv`

### Phase 4: Pipeline Infrastructure
All scripts written and tested in `src/python/direction3/`:

| Script | Purpose | Status |
|--------|---------|--------|
| `download_archs4.py` | Python-based H5 downloader with resume | Ready |
| `explore_archs4_h5.py` | Extract metadata & keyword search H5 | Ready |
| `extract_and_predict_archs4.py` | Extract counts, TMM normalize, predict tAge | Ready |
| `analyze_archs4_predictions.py` | Aging trajectory, drug effects, batch effects | Ready |
| `generate_report.py` | Auto-generate markdown report | Ready |
| `run_pipeline.sh` | Master pipeline (runs all above) | Ready |
| `auto_run_when_ready.sh` | Watcher that triggers pipeline on completion | Running |
| `geo_eutils_search.py` | GEO backup search via NCBI E-utilities | Completed |
| `process_gse288795.py` | Validation dataset processing | Completed |
| `process_gse230402.py` | Validation dataset processing | Completed |

## What Happens Next (Automatic)

When the H5 download completes (~45 min):

1. **Watcher detects completion** and runs `run_pipeline.sh`
2. **Metadata extraction**: Reads H5 metadata in batches of 10K samples
   - Output: `results/direction3/archs4_sample_metadata.csv.gz`
3. **Keyword searches**: Filters samples for target queries
   - rapamycin, caloric restriction, metformin, aging, dexamethasone, doxorubicin
   - Output: `results/direction3/archs4_search_*.csv`
4. **Count extraction**: Pulls expression matrix for selected samples (max 500 per set)
   - Chunked gene reading (5K genes at a time) to control memory
5. **TMM normalization**: edgeR via Rscript subprocess
6. **tAge prediction**: Aligns 10,487 Entrez IDs, imputes missing, centers, predicts
7. **Analysis**: Aging trajectory, drug-control differences, batch effect plots
8. **Report update**: `docs/REPORT_Direction3_ARCHS4_Mining.md` regenerated with ARCHS4 findings

## Manual Checkpoints (if needed)

```bash
# Check download progress
python3 src/python/direction3/check_download_status.py

# View active screen sessions
screen -ls

# Attach to download session
screen -r archs4_download

# Attach to watcher session
screen -r archs4_watcher

# Manually trigger pipeline (when H5 is complete)
bash src/python/direction3/run_pipeline.sh
```

## Known Limitations

1. **API failure**: ARCHS4 REST API endpoint is non-functional (404)
2. **Pre-normalized GEO counts**: GSE288795 and GSE230402 used pre-normalized counts (not raw + TMM)
3. **Batch effects**: ARCHS4 cross-study comparisons will require within-study normalization
4. **Download time**: 36GB file requires ~1 hour download at current speeds
5. **H5 memory**: Scripts use chunked reading; full matrix never loaded into RAM

## Files & Deliverables

- **Report**: `docs/REPORT_Direction3_ARCHS4_Mining.md`
- **Status**: `docs/Direction3_STATUS.md` (this file)
- **Results**: `results/direction3/`
- **Figures**: `figures/direction3/`
- **Scripts**: `src/python/direction3/`
- **GEO search**: `data/geo_expansion/geo_eutils_search_results.csv`
