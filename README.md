# sigmutselcovs

**Tissue-specific covariate matrices for mutation rate modeling**

## Overview

`sigmutselcovs` is a companion package to
[sigmutsel](https://github.com/alfaromurillo/sigmutsel) that builds
per-gene covariate matrices — replication timing, gene expression,
and chromatin accessibility — for any registered TCGA cancer type,
and downloads the source data needed to build them.

- **Registry-driven**: a per-cancer-type registry
  (`data/projects.json`) maps a TCGA study code (e.g. `COAD`,
  `BRCA`) to its GTEx tissue, TCGA gene expression, TCGA ATAC-seq,
  Roadmap Epigenomics chromatin tracks, and replication-timing
  sources.
- **Downloads**: idempotent, resumable downloaders for GDC (via its
  files API), ENCODE, Roadmap Epigenomics, and GEO — no more
  manually maintained manifests or wget scripts.
- **Builds**: `build_covariate_matrix` assembles the full, simple,
  and TCGA-only covariate matrices, with a column dictionary
  documenting every column's source, assembly, and units.
- **Validates**: automatic data-sanity and biological-plausibility
  checks on every build.
- **Fetches**: pre-built matrices can be pulled from OSF instead of
  rebuilding from ~10s of GB of raw tracks.

## Installation

### From source (development)

```bash
git clone https://github.com/alfaromurillo/sigmutselcovs.git
cd sigmutselcovs
pip install "sigmutsel @ git+https://github.com/alfaromurillo/sigmutsel.git"
pip install -e ".[dev]"
```

(`sigmutsel` is not yet on PyPI; install it from GitHub until it
is. `gdcfetch` is on PyPI, so it installs normally as part of
`pip install -e ".[dev]"` above.)

## Quick start

```bash
sigmutselcovs projects
sigmutselcovs download BRCA --data-dir brca_data
sigmutselcovs build    BRCA --data-dir brca_data
sigmutselcovs validate BRCA --data-dir brca_data
```

or from Python:

```python
from sigmutselcovs import download_covariates, build_covariate_matrix

download_covariates("BRCA", "brca_data")
full, simple, tcga = build_covariate_matrix(
    "BRCA", "brca_data", cache_matrices=True)
```

See `TUTORIAL.md` for a full walkthrough and `DEVELOPMENT.md` for
the module map and internals.

## Main functions

The five functions below cover the whole obtain-a-covariate-matrix
workflow; each is also available as a `sigmutselcovs <name>`
console command (`download`, `build`, `validate`, `fetch`,
`check-updates`).

| Function | Purpose |
|---|---|
| `download_covariates(project, data_dir, *, which=(...))` | Fetch source data (GDC gene expression, replication timing, Roadmap chromatin, TCGA ATAC-seq) for a registered project into `data_dir`. Idempotent and resumable — safe to rerun after a partial failure. `which` selects a subset of sources. |
| `build_covariate_matrix(project, data_dir, *, cache_matrices=False)` | Assemble the `full`, `simple`, and `tcga` covariate matrices from data already in `data_dir` (returns a `CovariateMatrices` NamedTuple; unpacks as a 3-tuple). Also writes a per-column dictionary (`cov_matrix_columns.csv`) and a build manifest when `cache_matrices=True`. Runs `validate_covariates` automatically at the end. |
| `validate_covariates(project, data_dir)` | Run data-sanity (NaN fractions, value ranges, column counts) and biological-plausibility (expression vs. replication timing, histone-mark direction, a housekeeping-gene panel) checks; returns a DataFrame of pass/warn/fail results. |
| `fetch_covariate_matrix(project, which="pca")` | Download a pre-built matrix from OSF instead of running download + build locally. Raises `CovariateArtifactsUnavailable` (a `FileNotFoundError`) with the equivalent download/build commands when nothing is published yet. |
| `check_updates()` | Compare each external source (GTEx, GDC, ENCODE, Roadmap, GEO) against its last-known state; run every few months to catch upstream data releases. |

`load_registry()` / `available_projects()` / `get_project(code)`
inspect the per-project source registry (`data/projects.json`)
without touching the network.

## Supported cancer types

Currently registered: `COAD`, `BRCA`. Adding a new type is a
registry entry — see `DEVELOPMENT.md` § "Adding a new cancer type".

## Data sources and citations

Every covariate comes from a public resource with its own citation
requirements. Short form, with the paper this package's methods
section should cite for each:

| Covariate | Source | Cite |
|---|---|---|
| Gene expression (baseline) | GTEx v10 | The GTEx Consortium (2020), *Science* 369:1318–1330 |
| Gene expression (per-sample) | GDC / TCGA | Grossman et al. (2016), *NEJM* 375:1109–1112 |
| Chromatin accessibility | TCGA ATAC-seq | Corces et al. (2018), *Science* 362:eaav1898 |
| Histone marks | Roadmap Epigenomics | Kundaje et al. (2015), *Nature* 518:317–330 |
| Replication timing | ENCODE Repli-seq | ENCODE Project Consortium (2012), *Nature* 489:57–74; Hansen et al. (2010), *PNAS* 107:139–144 |
| Replication timing (COAD) | GEO GSE137764 | Zhao et al. (2020), *Genome Biology* 21:76 |
| Gene coordinates | GENCODE | Frankish et al. (2019), *NAR* 47:D766–D773 |

**See `SOURCES.md` for the full citations (with DOIs), the exact
dataset/version/accession used per source, license terms, and
ENCODE's stricter three-part citation policy.**

## Related packages

- [sigmutsel](https://github.com/alfaromurillo/sigmutsel): mutation
  rate estimation and selection inference, the consumer of these
  covariate matrices.
- [gdcfetch](https://github.com/alfaromurillo/gdcfetch): the GDC
  search/download client this package's gene-expression and
  update-check code is built on. Useful on its own for downloading
  any GDC data — mutation calls, copy number, structural variants
  — not just what sigmutselcovs happens to need.

## License

MIT License — see LICENSE file for details.

## Support

For questions and issues, please use the [GitHub issue tracker](https://github.com/alfaromurillo/sigmutselcovs/issues).
