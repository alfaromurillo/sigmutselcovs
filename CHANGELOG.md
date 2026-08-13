# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-08-12

### Added
- Generalization from COAD-only to any registered TCGA cancer
  type, driven by a new per-project registry
  (`data/projects.json`) with entries for COAD and BRCA.
- `build_covariate_matrix(project, data_dir, ...)`: registry-driven
  assembly of the full/simple/tcga covariate matrices, reproducing
  the historical COAD pipeline bit-for-bit. Writes a column
  dictionary (`cov_matrix_columns.csv`) and a build manifest
  alongside the cached matrices.
- `download_covariates(project, data_dir, ...)`: idempotent,
  resumable downloaders for GDC gene expression (live API query,
  replacing checked-in manifests), replication timing (ENCODE
  Repli-seq or GEO MAT files), Roadmap Epigenomics chromatin, and
  TCGA ATAC-seq bigWig tarballs.
- Three replication-timing source types: `mat` (16-fraction
  transposed table), `fraction_bigwigs` (N ENCODE fraction tracks,
  early to late), and `wavelet` (single smoothed early/late
  signal), unlocking Repli-seq coverage for many more cell lines.
- `validate_covariates(project, ...)`: data-sanity and
  biological-plausibility checks (expression vs. replication
  timing, active/repressive chromatin marks, housekeeping-gene
  panel), run automatically at the end of every build.
- `fetch_covariate_matrix(project, which=...)`: download pre-built
  matrices from OSF instead of rebuilding from raw tracks (the OSF
  project itself is not yet populated).
- `check_updates()`: compare external covariate sources (GTEx, GDC,
  ENCODE, Roadmap, GEO) against their last-known state.
- `sigmutselcovs` console command: `projects`, `download`, `build`,
  `validate`, `fetch`, `check-updates`, `download-gtex`.
- `gdc.py` and `encode.py`: GDC files-API and ENCODE-portal clients
  used by the download layer.

### Changed
- `import_gtex` now defaults `mapping_path` to the packaged
  `gtex_tcga_mapping.json` when `columns` is a TCGA study code,
  matching its docstring's original promise.
- Replication-timing loaders were renamed from `repli_seq_hct` to
  `repliseq_source` (the old keyword still works for one release,
  with a `DeprecationWarning`), and gained `source_type`,
  `bin_size`, and `mrt_fraction_cols` parameters.
- The GTEx GCT now resolves from an explicit path, then the
  packaged data directory, then a user cache
  (`$XDG_CACHE_HOME/sigmutselcovs/gtex/`) — downloads never write
  into the installed package anymore.
- `coad_analysis/code/covariates.py`'s `build()` is now a thin
  wrapper over `build_covariate_matrix("COAD", ...)`; verified
  bit-for-bit identical to the pre-migration output.
- `requires-python` raised to `>=3.12`; added `pyarrow` and
  `requests` dependencies; adopted black/ruff at 70-character line
  length (matching sigmutsel).

### Fixed
- `load_or_generate_tcga_gexp_per_sample` was silently dropping its
  `tissue_type` argument.
- `annotate_rt_twidth` was ignoring its `bin_size` argument.

## [0.1.0] - 2026 (initial COAD-only version)

### Added
- Initial covariate-processing modules: GTEx/TCGA gene expression,
  replication timing (HCT116 MAT), chromatin (Roadmap and TCGA
  ATAC-seq bigWigs), data-quality checks (`fix_all` and friends),
  and shared genomic utilities. Downloading and per-project paths
  lived in the downstream `coad_analysis` repository.
