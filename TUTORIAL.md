# Tutorial: building covariates for BRCA

This walks through the full workflow for a cancer type that isn't
COAD, using BRCA (breast invasive carcinoma) as the worked example
— the first cancer type added after COAD.

## 1. List what's registered

```bash
sigmutselcovs projects
```

```
BRCA: Breast invasive carcinoma [gtex, gexp, repliseq, roadmap, atac]
COAD: Colon adenocarcinoma [gtex, gexp, repliseq, roadmap, atac]
```

Each bracketed name is a covariate source with an entry in
`data/projects.json` for that project. A project can have any
subset of the five (see `DEVELOPMENT.md` for the registry schema).

## 2. Download the source data

```bash
sigmutselcovs download BRCA --data-dir brca_data
```

This fetches, into `brca_data/` (laid out exactly like the
historical `coad_data/` tree — see `ProjectPaths`):

- **GTEx**: nothing to download here; the GTEx GCT resolves
  separately (see step 2b).
- **gexp**: queries the GDC files API live for
  `TCGA-BRCA` STAR-Counts gene expression, writes a manifest and
  sample sheet, and downloads ~1,231 files into
  `gene_expression/tcga/star_gene_counts/`.
- **repliseq**: for BRCA, six ENCODE bigWigs (the UW 6-fraction
  Repli-seq for MCF-7: G1b, S1–S4, G2) into
  `replication_timing/encode/`.
- **roadmap**: Roadmap Epigenomics fold-change signal bigWigs for
  the breast epigenomes E027, E028, E119, into
  `chromatin/roadmap/`. A few (EID, mark) combinations
  legitimately don't exist (e.g. E027 has no H3K27ac track) — those
  are skipped with a warning, not an error.
- **atac**: the TCGA ATAC-seq bigWig tarball for BRCA, extracted
  into `chromatin/tcga/`.

Downloads are idempotent — rerunning only fetches what's missing —
and resumable (`.part` files survive interruption).

To fetch only some sources, or preview without downloading:

```bash
sigmutselcovs download BRCA --data-dir brca_data --which repliseq roadmap
sigmutselcovs download BRCA --data-dir brca_data --dry-run
```

## 2b. The GTEx GCT

The GTEx gene-median-TPM matrix is shared across all projects and
resolved separately:

```bash
sigmutselcovs download-gtex
```

This checks (in order) an explicit path, the packaged data
directory, and `$XDG_CACHE_HOME/sigmutselcovs/gtex/`, downloading
into the cache if nothing is found. `build_covariate_matrix` calls
this resolution automatically, so you usually don't need to run it
by hand.

## 3. Build the covariate matrices

```bash
sigmutselcovs build BRCA --data-dir brca_data
```

This produces three matrices, cached under
`brca_data/covariate_matrices/`:

- `cov_matrix_full.parquet` — every covariate column
- `cov_matrix_simple.csv` — one representative column per source
  (GTEx expression, MRT, active/repressive chromatin marks)
- `cov_matrix_tcga.csv` — TCGA-only columns (expression, ATAC)
- `cov_matrix_columns.csv` — a dictionary, one row per column of
  `full`, describing its source, cell line/epigenome, assembly,
  units, and any `fix_all` transform applied
- `build_manifest.json` — registry hash, package version, per-source
  inventory

From Python:

```python
from sigmutselcovs import build_covariate_matrix

full, simple, tcga = build_covariate_matrix(
    "BRCA", "brca_data", cache_matrices=True)
```

Rebuilding is fast: each per-source block is cached as CSV/parquet
under `brca_data/`, so a rerun without `force_generation=True` just
reloads them and re-concatenates.

### Assemblies

BRCA's chromatin (Roadmap) and replication-timing (MCF-7) sources
are hg19; its ATAC-seq and gene-expression sources are hg38. The
builder picks the matching GENCODE GTF per source automatically
(from `sigmutsel.locations` by default) — this is why the RT and
Roadmap blocks of a BRCA matrix have a different (and typically
higher) NaN rate than the ATAC/expression blocks: genes only in the
hg19 GENCODE annotation but not hg38 (or vice versa) get NaN in the
other blocks.

## 4. Validate

```bash
sigmutselcovs validate BRCA --data-dir brca_data
```

Runs two tiers of checks against the raw (pre-`fix_all`) full
matrix:

- **Data sanity**: NaN fractions, constant columns, value ranges
  (MRT in [0, 1], expression/signal ≥ 0), column counts against the
  registry, index integrity.
- **Biological plausibility**: expression should anticorrelate with
  late replication timing; active histone marks (H3K4me3, H3K27ac,
  H3K9ac) and ATAC accessibility should correlate with expression;
  repressive marks (H3K9me3, H3K27me3) should anticorrelate; the
  same mark across Roadmap epigenomes should agree; a housekeeping
  panel (GAPDH, ACTB, RPL13A, B2M) should be highly expressed and
  early-replicating.

A single failed direction is a warning (a tissue can be genuinely
unusual); three or more together fail the run, since that usually
means an assembly or column-mapping bug. `build_covariate_matrix`
runs this automatically at the end of every build (log-only by
default; pass `strict_validation=True` to raise on failure).

## 5. Or skip straight to a pre-built matrix

If BRCA's matrices have already been published to OSF, you don't
need to download or build anything:

```python
from sigmutselcovs import fetch_covariate_matrix

pca = fetch_covariate_matrix("BRCA")          # small, default
full = fetch_covariate_matrix("BRCA", "full") # larger
```

Until the OSF project is populated, this raises
`CovariateArtifactsUnavailable` (a `FileNotFoundError` subclass)
with the download/build commands above in the error message.
