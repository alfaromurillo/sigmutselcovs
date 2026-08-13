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
is.)

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

## Supported cancer types

Currently registered: `COAD`, `BRCA`. Adding a new type is a
registry entry — see `DEVELOPMENT.md` § "Adding a new cancer type".

## Related packages

- [sigmutsel](https://github.com/alfaromurillo/sigmutsel): mutation
  rate estimation and selection inference, the consumer of these
  covariate matrices.

## License

MIT License — see LICENSE file for details.

## Support

For questions and issues, please use the [GitHub issue tracker](https://github.com/alfaromurillo/sigmutselcovs/issues).
