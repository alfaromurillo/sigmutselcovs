# Developing sigmutselcovs

Internals reference for anyone modifying this codebase.
Full data-source citations: `SOURCES.md`.

## Setup

```bash
pip install -e ".[dev]"
pip install -e ../sigmutsel   # until sigmutsel is on PyPI
pip install -e ../gdcfetch    # until gdcfetch is on PyPI
pytest                        # network-marked tests are skipped
```

## The generalized workflow

Everything is keyed on a TCGA study code registered in
`src/sigmutselcovs/data/projects.json` (currently COAD and BRCA):

```bash
sigmutselcovs projects
sigmutselcovs download BRCA --data-dir brca_data
sigmutselcovs build    BRCA --data-dir brca_data
sigmutselcovs validate BRCA --data-dir brca_data
sigmutselcovs fetch    BRCA          # pre-built from Zenodo (once published)
sigmutselcovs check-updates          # run every ~6 months
sigmutselcovs download-gtex
```

or from Python:

```python
from sigmutselcovs import download_covariates, build_covariate_matrix
download_covariates("BRCA", "brca_data")
full, simple, tcga = build_covariate_matrix("BRCA", "brca_data",
                                            cache_matrices=True)
```

`build_covariate_matrix` reproduces the historical COAD pipeline
exactly (same loaders, caches, and concatenation order), so running
it over an existing `coad_data/` tree is bit-for-bit identical to
the old `coad_analysis/code/covariates.py` `build()`.

## Module map

| Module | Role |
|--------|------|
| `registry.py` + `data/projects.json` | Per-cancer-type data sources |
| `paths.py` | `ProjectPaths`: data-dir layout (mirrors `coad_data/`) |
| `download.py` | Idempotent per-source downloaders + orchestrator |
| `encode.py` | ENCODE accession → S3 URL resolution |
| `builder.py` | `build_covariate_matrix` + column dictionary |
| `validate.py` | Data-sanity + biological-plausibility checks |
| `fetch.py` + `data/zenodo.json` | Pre-built matrices from Zenodo |
| `updates.py` + `data/sources.json` | Source update checks |
| `cli.py` | The `sigmutselcovs` console command |
| `covariates_gene_expression.py` | GTEx and TCGA expression loaders |
| `covariates_replication_timing.py` | Repli-seq: mat / fraction bigWigs / wavelet |
| `covariates_chromatin.py` | BigWig signal over gene bodies/promoters |
| `covariates_utilities.py` | GTF parsing, BED reading, annotation, PCA |
| `covariates_checks.py` | fix_all / fix_variance / fix_skewness |
| `covariates_locations.py` | Packaged-data paths (mapping JSON, GCT) |

GDC search/download/manifest functions (`search_files`,
`download_files`, `get_data_size`, `write_manifest`,
`write_sample_sheet`) come from the external
[gdcfetch](https://github.com/alfaromurillo/gdcfetch) package, not
a module here — it was split out since it's useful independently of
sigmutselcovs. See its `DEVELOPMENT.md` for the "two kinds of GDC
UUID" distinction (indexed files vs. publication-pinned blobs like
the ATAC tarballs) if you're touching anything GDC-related.

## Adding a new cancer type

1. Add a row to `data/projects.json`: GTEx mapping key (must exist
   in `gtex_tcga_mapping.json`), the TCGA project id, the ATAC
   tarball UUID from
   https://gdc.cancer.gov/about-data/publications/ATACseq-AWG,
   Roadmap EIDs, and a repliseq spec.
2. Replication timing source types:
   - `mat` — 16-fraction transposed table (GEO GSE137764: HCT116,
     H1, H9 only)
   - `fraction_bigwigs` — N per-fraction ENCODE bigWigs early→late
     (UW 6-fraction series: MCF-7, HeLa-S3, K562, GM12878, IMR-90,
     HepG2, keratinocyte, SK-N-SH, BJ, HUVEC, BG02)
   - `wavelet` — single smoothed early/late track (Gilbert series:
     LNCaP, Caki2, A549, NCI-H460, G401, SK-N-MC, T47D)
3. `pytest tests/test_registry.py` validates the row; then
   download, build, and `validate`.

## External data notes

- GTEx GCT (v10) resolves explicit → packaged → user cache and is
  fetched by `download-gtex`; a new GTEx major version renames
  tissue columns, so update `gtex_tcga_mapping.json` together with
  the GCT (and `sources.json`).
- GDC sample sheets and manifests are generated from the API — no
  more checked-in copies. A changed `check-updates` file count
  means a GDC data release touched the project: regenerate and
  rebuild the expression caches with `force_generation`.
- MCF-7 Repli-seq (and Roadmap) are hg19 → those blocks use the
  GENCODE v19 GTF; TCGA ATAC and the HCT116 MAT are hg38 → v38.
  `build_covariate_matrix` picks the GTF per source by the
  registry's `assembly`.

## Loader-specific gotchas

- `load_or_generate_rt_fractions` — per-fraction Repli-seq; the
  builder applies `clr_transform(...).add_prefix('clr_')` (CLR
  removes the compositional constraint) and gates CLR off below
  three fractions.
- `load_or_generate_tcga_gexp_per_sample` — wide
  gene×(barcode_metric) DataFrame; cached as Parquet due to
  thousands of columns (COAD: 3,084; BRCA: ~7,386).
- `load_or_generate_mean_tcga_gexp` defaults to
  `tissue_type="Tumor"` — cached CSVs generated before this default
  changed need `force_generation=True` to pick it up.
- **TCGA ATAC-seq is already per-sample**:
  `load_or_generate_chromatin_covariates` with
  `average_by_assay=False` (the default) produces one column per
  BigWig file (body + promoter each) — do not add separate
  per-sample handling for it.
- `run_pca_on_covariates` currently exists in BOTH
  `sigmutsel.utils` and `covariates_utilities.py` — a known
  duplication not yet resolved; until then change both or neither.
- Chromatin loading requires `pyBigWig` (Linux/Mac only). Never
  open bigWigs from URLs — download first (`download.py` does).
- GTF loading handles both gzip and plain text automatically.
