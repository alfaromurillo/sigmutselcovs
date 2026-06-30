# sigmutselcovs

## Project overview

Companion package to `sigmutsel` for building gene-level covariate
matrices used in mutation rate modeling. Wraps genomic data sources:
GTEx gene expression, TCGA gene expression, Repli-seq replication
timing, and BigWig histone/chromatin tracks.

## Setup

```bash
pip install -e .
pip install -e ../sigmutsel   # sigmutsel.locations is a dependency
```

External data files (not bundled; paths set per-project):
- GTEx GCT file (v10): `GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct`
- TCGA STAR count TSVs: downloaded from GDC
- Repli-seq `.mat` file (HCT116): from GEO GSE137764
- BigWig files: Roadmap Epigenomics or TCGA ATAC-seq

## Key modules

| Module | Role |
|--------|------|
| `covariates_gene_expression.py` | GTEx and TCGA expression loaders |
| `covariates_replication_timing.py` | Repli-seq MRT per gene |
| `covariates_chromatin.py` | BigWig signal over gene bodies/promoters |
| `covariates_utilities.py` | GTF parsing, BED reading, genomic annotation, PCA |
| `covariates_locations.py` | Default paths for GTEx GCT and mapping JSON |

## Usage pattern

Each covariate module has a `load_or_generate_*` wrapper that reads
from a cached CSV on subsequent calls:

```python
from sigmutselcovs.covariates_replication_timing import load_or_generate_mrt
from sigmutsel.locations import location_gencode38_annotation

mrt = load_or_generate_mrt(
    location_mrt_csv,       # cache path
    location_repli_seq_hct, # raw data path
    location_gencode38_annotation,
    force_generation=False)
```

See `coad_analysis/code/covariates.py` for a full worked example.

## New loaders (added 2026-06)

- `load_or_generate_rt_fractions` — per-fraction Repli-seq; apply
  `clr_transform(...).add_prefix('clr_')` before use (CLR removes
  compositional constraint; prefix makes transform explicit downstream)
- `load_or_generate_tcga_gexp_per_sample` — wide gene×(barcode_metric)
  DataFrame; cached as Parquet (not CSV) due to ~3,084 columns
- `import_tcga_gene_expression` takes `tissue_type` kwarg; default None
  (all samples); `load_or_generate_mean_tcga_gexp` defaults to
  `tissue_type="Tumor"` — cached CSVs generated before this change need
  `force_generation=True`
- `clr_transform` is in `covariates_utilities.py`

## TCGA ATAC-seq is already per-sample

`load_or_generate_chromatin_covariates` with `average_by_assay=False`
(the default) produces one column per BigWig file. TCGA ATAC has 81
files × 2 regions = 162 columns — do NOT add per-sample code for it.

## Notes

- `run_pca_on_covariates` lives in `sigmutsel.utils`, not here —
  import from there to avoid duplication
- Chromatin loading requires `pyBigWig` (Linux/Mac only)
- GTF loading handles both gzip and plain text automatically
