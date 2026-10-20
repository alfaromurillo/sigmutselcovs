# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `build_covariate_matrix` records the source block of every column
  of the full matrix in `full.attrs["column_sources"]` (the same
  sources as the column dictionary), and `combine_with_generic`
  carries both matrices' records over.
- `covariate_column_blocks(matrix)` groups those sources into kinds
  of measurement (`SOURCE_BLOCKS`: gtex, expression, replication
  timing, atac, chromatin), for callers that treat a gene's missing
  covariates by block.
- `fetch_covariate_matrix(project, "full", columns=..., blocks=...)`
  reads only the named columns, or whole blocks such as
  `["chromatin", "atac"]`, of a published full matrix by HTTP range
  requests (`HTTPRangeFile`), instead of downloading the whole file:
  3 columns of a 538 MB matrix transfer about 3.4 MB, mostly the
  parquet footer. Names and blocks are checked against the record's
  column dictionary before any large transfer. Not cached, and not
  checksum-verified, since only parts of the file are read. The
  `fetch` command takes the same selection as `--columns` and
  `--blocks`, and `--out FILE` (.parquet or .csv) to keep the result.
  README leads with downloading pre-built matrices.
- Opt-in chromatin summaries. Per-gene windows `cds_span` (first to
  last coding base, from the GTF's CDS features merged across
  transcripts; NaN without a CDS) and `cds_span2k` (that span
  +-2 kb) next to `body` and `promoter`; `transform="log1p"`, the
  per-base mean of log1p(signal) over the bases a track covers
  (columns `<label>_<window>_log`); and peak coverage, the fraction
  of a window inside the union of a track's narrowPeak/broadPeak
  calls (columns `<label>_<peaktype>_<window>`). Registry fields
  `windows`, `transform` and `peaks` (Roadmap) / `peak_tracks`
  (ENCODE) on `roadmap` and `encode_chromatin`, validated; peak
  downloaders `download_roadmap_peaks` and
  `download_encode_chromatin_peaks`; functions
  `summarize_bigwig_windows`, `summarize_peak_coverage`,
  `gene_windows`, `load_or_generate_peak_covariates` and
  `chromatin_cache_path`; `load_gene_bodies_from_gtf(...,
  with_cds_span=True)`. `transform` also takes a list, e.g.
  `[null, "log1p"]`, so one source emits both linear and log
  columns. Each option set is cached in its own file.
  The defaults are unchanged: default builds are bit for bit the
  same.

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
  matrices from Zenodo instead of rebuilding from raw tracks (no
  records are published yet).
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
- The original per-cohort pipeline's `covariates.py` `build()` is
  now a thin wrapper over `build_covariate_matrix("COAD", ...)`;
  verified bit-for-bit identical to the pre-migration output.
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
  lived in a downstream, cohort-specific pipeline repository.
