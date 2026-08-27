"""Unsupervised sizing for the PCA-reduced covariate matrix this
package publishes (e.g. to Zenodo).

PCA is unsupervised: it is fit on a covariate matrix alone and never
sees any downstream fitting outcome. That means "how many components
does this covariate matrix need" is a different question from "how
many does one particular downstream regression benefit from" -- the
second question is task-specific and belongs to whichever downstream
analysis does that fitting (e.g. via its own cross-validated
component-count selection), not to this package. This module answers
only the first, task-agnostic question, sizing the artifact generously
enough to remain useful for downstream tasks this package was not
built with in mind.

The criterion is cumulative explained variance: the smallest number
of components whose cumulative explained-variance ratio reaches a
threshold (99% by default). A lower threshold such as 95% was
considered and rejected -- across every TCGA cohort checked so far,
95% cumulative variance falls at or below the component count a
downstream regression already found useful on its own held-out
cross-validation, which would understate the published resource for
exactly the analysis this package was built to support. 99% clears
every checked cohort's own task-specific optimum with substantial
headroom.

Needs ``sigmutsel`` installed (imported lazily, only inside
:func:`build_pca_artifact`) -- this package's core covariate-building
does not otherwise depend on it, matching the decoupling decision in
``sigmutselcovs/DEVELOPMENT.md``.
"""

import json
import logging

import numpy as np
import pandas as pd

from .paths import ProjectPaths

logger = logging.getLogger(__name__)

DEFAULT_VARIANCE_THRESHOLD = 0.99


def build_pca_artifact(
    cov_matrix_full: pd.DataFrame,
    *,
    variance_threshold: float = DEFAULT_VARIANCE_THRESHOLD,
    **pca_kwargs,
) -> pd.DataFrame:
    """Reduce a covariate matrix to the components needed for
    ``variance_threshold`` cumulative explained variance.

    Fits PCA at full rank once (``sigmutsel.utils.
    run_pca_on_covariates`` with ``n_components=None``) and slices
    the leading columns -- PCA components are nested, so the first
    *k* columns of a full-rank fit equal a direct *k*-component fit,
    up to sign. No reference to any fitting task's outcome is made
    anywhere in this function.

    Parameters
    ----------
    cov_matrix_full : pandas.DataFrame
        Gene-indexed covariate matrix, e.g. ``sigmutselcovs.
        build_covariate_matrix(...).full``.
    variance_threshold : float, default 0.99
        Cumulative explained-variance ratio to reach. Must be in
        (0, 1).
    **pca_kwargs
        Forwarded to ``run_pca_on_covariates`` (e.g. ``standardize``,
        ``dropna``, ``random_state``).

    Returns
    -------
    pandas.DataFrame
        Gene-indexed PC scores, truncated to the selected number of
        components. ``result.attrs`` carries ``n_components_selected``,
        ``variance_threshold``, ``cumulative_variance_at_selection``,
        ``raw_columns`` (the input's column count), and
        ``explained_variance_ratio`` (truncated to match).

    Raises
    ------
    ValueError
        If ``variance_threshold`` is not in (0, 1).
    """
    if not 0 < variance_threshold < 1:
        raise ValueError(
            f"variance_threshold must be in (0, 1), got "
            f"{variance_threshold!r}"
        )

    from sigmutsel.utils import run_pca_on_covariates

    scores_full = run_pca_on_covariates(
        cov_matrix_full, n_components=None, **pca_kwargs
    )
    evr = scores_full.attrs["explained_variance_ratio"]
    cumvar = np.cumsum(evr)

    n_components = (
        int(np.searchsorted(cumvar, variance_threshold)) + 1
    )
    n_components = min(n_components, scores_full.shape[1])

    scores = scores_full.iloc[:, :n_components].copy()
    scores.attrs["explained_variance_ratio"] = evr[:n_components]
    scores.attrs["n_components_selected"] = n_components
    scores.attrs["variance_threshold"] = variance_threshold
    scores.attrs["cumulative_variance_at_selection"] = float(
        cumvar[n_components - 1]
    )
    scores.attrs["raw_columns"] = int(cov_matrix_full.shape[1])

    logger.info(
        "PCA artifact: %d components reach %.1f%% variance "
        "(threshold %.1f%%), from %d raw columns",
        n_components,
        100 * scores.attrs["cumulative_variance_at_selection"],
        100 * variance_threshold,
        scores.attrs["raw_columns"],
    )
    return scores


def _package_version(name: str) -> str:
    try:
        from importlib.metadata import version

        return version(name)
    except Exception:  # noqa: BLE001 - version is informational only
        return "unknown"


def save_pca_artifact(
    scores: pd.DataFrame, paths: ProjectPaths, *, project: str
) -> None:
    """Write a :func:`build_pca_artifact` result and its metadata to
    ``<data_dir>/covariate_matrices/`` (``cov_matrix_pca.parquet`` +
    ``pca_manifest.json``), matching the file-naming convention used
    by ``build_covariate_matrix``'s own cached outputs.
    """
    paths.matrices_dir.mkdir(parents=True, exist_ok=True)
    # .attrs (e.g. explained_variance_ratio, an ndarray) isn't
    # parquet-serializable and pandas' own best-effort attrs
    # round-trip chokes on it rather than silently dropping it (at
    # least on this pandas version) -- write a clean copy instead;
    # the scalar metadata that matters is already in pca_manifest.json.
    to_write = scores.copy()
    to_write.attrs = {}
    to_write.to_parquet(paths.matrix_pca_parquet)

    manifest = {
        "project": project,
        "sigmutselcovs_version": _package_version("sigmutselcovs"),
        "sigmutsel_version": _package_version("sigmutsel"),
        "shape": list(scores.shape),
        "raw_columns": scores.attrs["raw_columns"],
        "variance_threshold": scores.attrs["variance_threshold"],
        "n_components_selected": scores.attrs[
            "n_components_selected"
        ],
        "cumulative_variance_at_selection": scores.attrs[
            "cumulative_variance_at_selection"
        ],
        "built": pd.Timestamp.now(tz="UTC").isoformat(),
    }
    paths.pca_manifest_json.write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    logger.info(
        "Wrote PCA artifact for %s (%d components) to %s",
        project,
        scores.attrs["n_components_selected"],
        paths.matrix_pca_parquet,
    )
