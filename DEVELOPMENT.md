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

`--data-dir` accepts any path; a downstream pipeline's own layout
convention (e.g. `data/<CODE>`) is not enforced here. `build_covariate_matrix`
reproduces the historical COAD pipeline exactly (same loaders, caches,
and concatenation order), so running it over an existing
`coad_data/`-layout tree reproduces the original pre-migration
pipeline's `build()` output bit-for-bit.

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
| `covariates_chromatin.py` | BigWig signal and peak coverage over gene windows |
| `chromatin_stream.py` | Parallel, resumable chromatin caches with no bigWig kept |
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
the same way as any cohort. It is the complete tissue-agnostic matrix:
pooled chromatin, pan-tissue expression and pooled replication timing.
It differs from a cohort's row in four ways:

- `gtex.mapping_key` resolves (via `gtex_tcga_mapping.json`, same
  mechanism as any code) to *every* GTEx tissue column rather than
  one representative tissue, and `gtex.reduce: "median"` collapses
  them into a single `gtex_pantissue_median` column — see
  `GtexSpec.reduce` and `import_gtex`'s `reduce` parameter. `None`
  (the default for every real cohort) keeps today's one-column
  behavior.
- `gexp` and `atac` are `null` — these are TCGA-tumor-specific data
  types with no tissue-agnostic equivalent.
- `repliseq` is a `pooled_mrt` source: one profile per replication-
  timing experiment (`mat`, `fraction_bigwigs` or `wavelet`, each on
  its own assembly), giving a single `rt_pool_z` column. Each
  profile's per-gene timing is oriented larger = later (a wavelet's
  log2 early/late is negated), z-scored across genes, averaged within
  a biosample (`cell_line`) and then across biosamples, so each
  biosample has one vote (`pool_rt_profiles`). Timing is largely
  conserved across cell types, which is what makes a pool a sensible
  tissue-agnostic covariate. Per-profile results are cached under
  `replication_timing/pool/`. Bigwig sources are binned on chr1-22
  and chrX and every chromosome's genes are timed (2026-10-21; before,
  autosomes only), so chrX genes get a pooled time from the profiles
  that cover chrX; an autosome-only table (the HCT116 MAT) leaves them
  out of its own profile. Per-gene caches written before keep their
  autosome-only genes until regenerated.
- `roadmap`/`encode_chromatin` are populated with a large pool of
  epigenomes/tracks spanning many tissues (rather than one matched
  tissue), pooled through the same PCA/concatenation machinery as
  any cohort's matrix. Since 2026-10-21 they are summarised exactly as
  a cohort's rows (the registry defaults, "Chromatin summaries"
  below), one column per track (`average_by_assay` false), over all
  32 Roadmap consolidated marks. The pool is too large for disk:
  build its caches with `sigmutselcovs stream-chromatin GENERIC
  --source roadmap` (and `encode_chromatin`), which downloads,
  summarises and deletes a few tracks at a time, resumably
  (`chromatin_stream.py`); `build` then reads the caches.

`combine_with_generic(cov_matrix_full, generic_matrix)` (in
`builder.py`, exported from the package root) concatenates a
cohort's own `cov_matrix_full` with `GENERIC`'s (outer join on
`ensembl_gene_id`, same as any other source block) — this is the
supported way to get a matrix combining tissue-specific and
pan-tissue signal; it is not automatic, since not every downstream
use wants both. Pass `drop_own_tracks=True` to leave out of the pool
every track the cohort's matrix already summarises, matched by track
identity (Roadmap epigenome and mark, ENCODE file accession;
`chromatin_track_identity`), whatever the windows and transforms.

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

## Chromatin summaries

Every chromatin track (a Roadmap or ENCODE bigWig, a TCGA ATAC
bigWig) is reduced to a few numbers per gene. The spec defaults are
two: the mean signal over the gene body and over the promoter (TSS
-2 kb / +200 bp, strand-aware), in `summarize_bigwig_over_genes` and
`summarize_promoter_signal`, columns `<label>_body` and
`<label>_promoter`. **The shipped registry's `defaults` block puts
every row, a cohort's and `GENERIC`'s, on one method instead
(2026-10-21): `windows: ["cds_span"]`, `transform: "log1p_mean"`,
`min_completeness: 0.9`, for `roadmap`, `encode_chromatin` and
`atac` alike.** log1p of the mean signal over the CDS span is the
summary that best reproduces the chromatin principal components a
widely used selection method ships (cross-validated R^2 0.855 against
0.754 for the gene body), and it was the better covariate on held-out
mutation counts. `build_covariate_matrix(chromatin_windows=[...])`
(CLI `--chromatin-windows`) overrides the windows of every chromatin
source at once, e.g. `cds_span promoter`.

- **`windows`** (default `["body", "promoter"]`): any of `body`,
  `promoter`, `cds_span` (first to last coding base: the gene
  without its UTRs, from the GTF's `CDS` features merged across
  transcripts, via `load_gene_bodies_from_gtf(...,
  with_cds_span=True)`; NaN for a gene with no CDS) and
  `cds_span2k` (the CDS span +-2 kb). All 0-based half-open
  (`gene_windows`), starts clipped at 0.
- **`transform`** (default `null`): `"log1p_mean"` replaces each mean
  by log1p of the mean (clipped at 0), columns ending `_logmean`; the
  skewness fix (`fix_all(skip_skewness=...)`) leaves those columns
  alone, as they are already on their scale. `"log1p"` replaces each mean by
  the mean of log1p(signal) *per base* (not the log of the mean)
  over the bases the track covers, the bases
  `pyBigWig.stats(..., exact=True)` counts. `pyBigWig.stats` cannot
  compute it: `_log1p_interval_means` reads `bw.values(...,
  numpy=True)` in chunks of at most 8 Mb holding whole windows and
  sums by cumulative sums. Columns gain `_log`. A list, e.g.
  `[null, "log1p"]`, emits each track's columns once per transform
  (linear then log, in the order given); a scalar or a one-element
  list behaves as before.
- **`min_completeness`** (default `null`): drop a summary column with
  values for fewer than that fraction of the source's genes (a track
  missing whole chromosomes; the hg19/hg38 gene-set gaps that per-mark
  averaging used to hide).
- **Peaks**: `roadmap.peaks` (a list of `narrowPeak`/`broadPeak`)
  or `encode_chromatin.peak_tracks` (`{label, accession,
  peak_type}`) add `summarize_peak_coverage` columns
  `<label>_<peaktype>_<window>` (e.g.
  `e075_h3k4me3_narrowpeak_body`): the fraction of each window's
  bases inside the union of the track's peaks (overlapping peaks
  merged first). The transform does not apply to them. Peak files
  live under `chromatin/<source>/peaks/`
  (`download_roadmap_peaks`, `download_encode_chromatin_peaks`),
  and their column source is `roadmap_peaks`/
  `encode_chromatin_peaks`, both in the `chromatin` block.

Edge cases: without a transform, `body`/`promoter` are the original
functions, which give NaN for a window running past the chromosome
end; the CDS windows and the log transform clip it to the
chromosome instead. Peak windows are clipped the same way, with the
packaged UCSC chromosome sizes of the source's assembly
(`load_chrom_sizes("hg19")`, `data/chrom_sizes/`); a direct call to
`summarize_peak_coverage` clips only when given `chrom_sizes`. A BED
file does not list its chromosomes, so a chromosome with no peaks
has coverage 0 (with a warning).

**Caches.** `chromatin_cache_path` names one cache per option set:
the defaults keep `chromatin_covs_<source>.csv`, others get e.g.
`chromatin_covs_roadmap_cds_span-cds_span2k_log1p.csv`,
`chromatin_covs_roadmap_body-promoter_linear-log1p.csv` (a
transform list, `null` spelled `linear`) or
`chromatin_covs_roadmap_peaks_narrowpeak_body-promoter.csv`, so a
build never reads a cache built under other options. Each cache
records the files it was built from (`<cache>.tracks.json`) and is
rebuilt when they change.

**Downstream.** `_collapse_by_assay` (`average_by_assay`) keys on
the mark plus the whole summary suffix (peak type, window, `_log`),
so `h3k4me3_cds_span_log` and `h3k4me3_narrowpeak_body` are
averaged separately. `cov_matrix_simple`'s patterns skip `_log`
columns unless the pattern itself ends in `_log`. `validate.py`
sizes its column-count checks by the windows and peaks, and its
direction checks fall back to the log-means when a build has no
linear means.

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
