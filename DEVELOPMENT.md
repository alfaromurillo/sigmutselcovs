# Developing sigmutselcovs

Internals reference for anyone modifying this codebase.
Full data-source citations: `SOURCES.md`.

## Setup

```bash
pip install -e ".[dev]"       # gdcfetch resolves from PyPI
pytest                        # network-marked tests are skipped
```

No sigmutsel install needed: this package has no dependency on it
(GENCODE GTFs are fetched/cached directly, see
`download.ensure_gencode_gtf`). Only install sigmutsel separately
if you're also working on the PCA/Riemannian modeling path (see
"Loader-specific gotchas" below).

## The generalized workflow

Everything is keyed on a TCGA study code registered in
`src/sigmutselcovs/data/projects.json` (13 TCGA cohorts as of
2026-09, plus the tissue-agnostic `GENERIC` pseudo-project -- see
"Generic (pan-tissue) covariates" below):

```bash
sigmutselcovs projects
sigmutselcovs download BRCA --data-dir data/BRCA
sigmutselcovs build    BRCA --data-dir data/BRCA
sigmutselcovs validate BRCA --data-dir data/BRCA
sigmutselcovs fetch    BRCA          # pre-built from Zenodo (once published)
sigmutselcovs check-updates          # run every ~6 months
sigmutselcovs download-gtex
```

or from Python:

```python
from sigmutselcovs import download_covariates, build_covariate_matrix
download_covariates("BRCA", "data/BRCA")
full, simple, tcga = build_covariate_matrix("BRCA", "data/BRCA",
                                            cache_matrices=True)
```

`--data-dir` accepts any path; `data/<CODE>` is `tcga_analysis`'s
convention (see its `CLAUDE.md`), not something this package
enforces. `build_covariate_matrix` reproduces the historical COAD
pipeline exactly (same loaders, caches, and concatenation order), so
running it over an existing `coad_data/`-layout tree is bit-for-bit
identical to the old `coad_analysis/code/covariates.py` `build()`.

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

## Adding a new cohort

1. Add a row to `data/projects.json`: GTEx mapping key (must exist
   in `gtex_tcga_mapping.json`), the TCGA project id, the ATAC
   tarball UUID from
   https://gdc.cancer.gov/about-data/publications/ATACseq-AWG,
   Roadmap EIDs, a repliseq spec, and (when no matched Roadmap
   epigenome exists for the tissue, or to add DNase-seq or a mark
   Roadmap's default panel doesn't cover) `encode_chromatin` tracks
   -- see "Generic (pan-tissue) covariates" below for the ENCODE
   portal search pattern used to find these.
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

## Adding a new covariate source

Different from "adding a new cohort" above, which wires a new
*cohort* into the six existing source kinds (gtex, gexp, repliseq,
roadmap, atac, encode_chromatin). This is for a genuinely new *kind*
of data — DNA methylation, Hi-C compartments, a new track catalog
entirely, etc. Touches more of the codebase; budget a full pass
through it. (A ChIP mark or DNase track Roadmap doesn't cover is
*not* a new kind — that's `encode_chromatin`, an existing source; see
"Adding a new cohort" above and "Generic (pan-tissue) covariates"
below.)

1. **Registry** (`registry.py`): add a frozen dataclass for the
   source's parameters (mirror `RoadmapSpec`/`AtacSpec`), add it as
   an `Optional[YourSpec] = None` field on `ProjectSpec`, and wire
   it through `load_registry()`'s `source(...)` helper. If it needs
   required sub-fields, validate them in `validate_registry()`
   (mirror the `repliseq.type` checks).
2. **Paths** (`paths.py`): add the source's directory/cache-file
   paths to `ProjectPaths` (mirror `roadmap_covs_csv` /
   `roadmap_bigwig_files()`).
3. **Download** (`download.py`): add a downloader function, and add
   the source's key to `_WHICH` so `sigmutselcovs download PROJECT
   --which yoursource` and the `download_covariates()` orchestrator
   both pick it up.
4. **Per-gene loader**: add a new `covariates_yoursource.py` module
   with a `load_or_generate_...()` function following the
   cache-first pattern in `covariates_chromatin.py` /
   `covariates_replication_timing.py` — check for a cached
   CSV/parquet first, else compute from the raw downloaded files and
   write the cache.
5. **Builder wiring** (`builder.py`): add the source's key to
   `_SOURCES`; add a block in `build_covariate_matrix()` following
   the `if "roadmap" in selected: ...` pattern (skip-with-warning
   when the spec is `None` or files are missing, append to `frames`
   and to `blocks` for the column dictionary); add a case to
   `_describe_block()` so the column dictionary gets a real
   description/units instead of the generic fallback.
6. **Validation** (`validate.py`): if the source has a natural value
   range, add a check to `_sanity_tier`. More importantly, add a
   `_check_direction(...)` call in `_biology_tier()` stating which
   way the new source should correlate with expression (or `mrt`, or
   another block) — this is what turns a broken download or a wrong
   sign into a validation failure instead of a silently wrong
   number.
7. **Docs**: cite the source in `SOURCES.md` (new `##` section,
   following the existing pattern: what it is, exact accession/
   version, license, how to cite) and add it to `README.md`'s data
   sources table.
8. **Tests**: the registry dataclass + `validate_registry` case
   (`test_registry.py`), the per-gene loader against a tiny
   synthetic file, and the builder's skip-with-warning + column-
   dictionary behavior for the new source (`test_builder.py`).

## Generic (pan-tissue) covariates

`GENERIC` is a pseudo-project code (not a real TCGA study) registered
the same way as any cohort, with three differences:

- `gtex.mapping_key` resolves (via `gtex_tcga_mapping.json`, same
  mechanism as any code) to *every* GTEx tissue column rather than
  one representative tissue, and `gtex.reduce: "median"` collapses
  them into a single `gtex_pantissue_median` column — see
  `GtexSpec.reduce` and `import_gtex`'s `reduce` parameter. `None`
  (the default for every real cohort) keeps today's one-column
  behavior.
- `gexp`, `atac`, `repliseq` are `null` — these are TCGA-tumor-
  specific data types with no tissue-agnostic equivalent.
- `roadmap`/`encode_chromatin` are populated with a large pool of
  epigenomes/tracks spanning many tissues (rather than one matched
  tissue), pooled through the same PCA/concatenation machinery as
  any cohort's matrix.

`combine_with_generic(cov_matrix_full, generic_matrix)` (in
`builder.py`, exported from the package root) concatenates a
cohort's own `cov_matrix_full` with `GENERIC`'s (outer join on
`ensembl_gene_id`, same as any other source block) — this is the
supported way to get a matrix combining tissue-specific and
pan-tissue signal; it is not automatic, since not every downstream
use wants both.

**Finding tracks for `encode_chromatin`** (used for both a cohort's
own registry row and `GENERIC`'s pool): search the ENCODE portal for
released, GRCh38, bigWig files —
`https://www.encodeproject.org/search/?type=File&status=released&assembly=GRCh38&file_format=bigWig&assay_title=Histone+ChIP-seq&target.label=<MARK>`
(swap `target.label` per mark; DNase-seq uses
`assay_title=DNase-seq` with no `target.label`). Cross-check against
EpiMap's metadata table
(https://personal.broadinstitute.org/cboix/epimap/metadata/main_metadata_table.tsv)
for biosample coverage before committing to a track. Not every
mark has broad tissue coverage on native GRCh38 — some are cell-line/
differentiated-cell only, in which case they belong in `GENERIC`'s
pool (any coverage helps a pooled estimate) but not in a specific
cohort's matched-tissue row (a poor tissue match is worse than no
data for that source).

**Version/refit policy.** `GENERIC`'s precision comes from the
number of pooled epigenomes/tracks, so it should be rebuilt whenever
a meaningful number of new tracks become available (e.g. after
bringing up several new cohorts' `encode_chromatin` rows) — but a
rebuild changes `GENERIC`'s columns for every cohort that concatenates
it, which is a consequential enough change to need its own
before/after check, not silent regeneration. Every
`build_covariate_matrix(..., cache_matrices=True)` call already
stamps `registry_sha256_16` (a hash of the whole `projects.json` at
build time) into `build_manifest.json` — for `GENERIC` specifically,
this hash *is* its version identifier: it changes exactly when
`GENERIC`'s registry row (or the shared `defaults` block) changes, so
two `GENERIC` builds sharing a `registry_sha256_16` are guaranteed to
carry the same pooled tracks. Anyone concatenating `GENERIC` into a
cohort's matrix (`combine_with_generic`) should keep the
`registry_sha256_16` from the `GENERIC` build they used alongside the
result, so a later `GENERIC` rebuild is a visible, deliberate
version change rather than a silent difference. Treat a `GENERIC`
rebuild the same as any other source-data version bump affecting a
published/production matrix: revalidate downstream fits against the
previous version before adopting the new one, rather than assuming a
richer pool is automatically better.

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
- PCA over a built covariate matrix: `from sigmutsel.utils import
  run_pca_on_covariates` (sigmutsel needs to be installed
  separately for this -- as of 2026-08-13 this package no longer
  depends on sigmutsel for anything; an unused duplicate of this
  function in `covariates_utilities.py` was also removed the same
  day).
- **Resolved (2026-08-27): sizing the published PCA artifact does
  live here now**, via `pca_artifact.build_pca_artifact` /
  `save_pca_artifact`. `sigmutsel` stays an optional dependency
  (the `pca` extra, `pip install sigmutselcovs[pca]`), imported
  lazily inside `build_pca_artifact` only -- the core
  covariate-building path still doesn't need it. The sizing
  criterion is deliberately *not* tied to any downstream fitting
  task: cumulative explained variance (99% default) computed from
  the covariate matrix alone, since PCA never sees a fitting
  outcome and a task-tuned component count would understate the
  resource for any other downstream use. 95% was considered and
  rejected -- across every TCGA cohort checked, it fell at or below
  the component count a specific downstream regression already
  found useful on its own cross-validation, which would have
  understated the resource for exactly that use; 99% clears every
  checked cohort's task-specific optimum with substantial headroom
  instead. See `build_pca_artifact`'s docstring for the mechanism.
- Chromatin loading requires `pyBigWig` (Linux/Mac only). Never
  open bigWigs from URLs — download first (`download.py` does).
- GTF loading handles both gzip and plain text automatically.
