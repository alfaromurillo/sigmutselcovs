"""Tests for the unsupervised PCA-artifact sizing module."""

import json

import numpy as np
import pandas as pd
import pytest

from sigmutselcovs.paths import project_paths
from sigmutselcovs.pca_artifact import (
    DEFAULT_VARIANCE_THRESHOLD,
    build_pca_artifact,
    save_pca_artifact,
)


def _synthetic_cov_matrix(n_genes=200, n_cols=12, seed=0):
    rng = np.random.default_rng(seed)
    index = [f"ENSG{i:08d}" for i in range(n_genes)]
    data = rng.normal(size=(n_genes, n_cols))
    return pd.DataFrame(
        data, index=index, columns=[f"cov{i}" for i in range(n_cols)]
    )


def test_build_pca_artifact_reaches_threshold_minimally():
    cov = _synthetic_cov_matrix()
    scores = build_pca_artifact(
        cov, variance_threshold=0.90, random_state=0
    )

    n = scores.attrs["n_components_selected"]
    cumvar_at_n = scores.attrs["cumulative_variance_at_selection"]
    assert cumvar_at_n >= 0.90
    if n > 1:
        evr = scores.attrs["explained_variance_ratio"]
        # one fewer component must not already clear the threshold --
        # confirms n is the *smallest* count reaching it, not just any
        assert np.cumsum(evr)[n - 2] < 0.90
    assert scores.shape == (cov.shape[0], n)


def test_build_pca_artifact_higher_threshold_needs_more_components():
    cov = _synthetic_cov_matrix()
    low = build_pca_artifact(
        cov, variance_threshold=0.80, random_state=0
    )
    high = build_pca_artifact(
        cov, variance_threshold=0.99, random_state=0
    )
    assert (
        high.attrs["n_components_selected"]
        >= low.attrs["n_components_selected"]
    )


def test_build_pca_artifact_attrs_present():
    cov = _synthetic_cov_matrix()
    scores = build_pca_artifact(cov, random_state=0)
    assert (
        scores.attrs["variance_threshold"]
        == DEFAULT_VARIANCE_THRESHOLD
    )
    assert scores.attrs["raw_columns"] == cov.shape[1]
    assert (
        len(scores.attrs["explained_variance_ratio"])
        == scores.attrs["n_components_selected"]
    )


def test_build_pca_artifact_nested_components_match_full_fit():
    """The first k columns of a low-threshold fit must equal the
    first k columns of a high-threshold fit (up to per-component
    sign) -- both are truncations of the same full-rank decomposition,
    the nested-components property this module's design relies on."""
    cov = _synthetic_cov_matrix()
    small = build_pca_artifact(
        cov, variance_threshold=0.80, random_state=0
    )
    big = build_pca_artifact(
        cov, variance_threshold=0.99, random_state=0
    )
    for col in small.columns:
        a = small[col].to_numpy()
        b = big[col].to_numpy()
        assert np.allclose(a, b, atol=1e-8) or np.allclose(
            a, -b, atol=1e-8
        )


def test_build_pca_artifact_invalid_threshold_raises():
    cov = _synthetic_cov_matrix()
    with pytest.raises(ValueError, match="variance_threshold"):
        build_pca_artifact(cov, variance_threshold=1.5)
    with pytest.raises(ValueError, match="variance_threshold"):
        build_pca_artifact(cov, variance_threshold=0.0)


def test_save_pca_artifact_writes_parquet_and_manifest(tmp_path):
    cov = _synthetic_cov_matrix()
    scores = build_pca_artifact(
        cov, variance_threshold=0.90, random_state=0
    )
    paths = project_paths(tmp_path)

    save_pca_artifact(scores, paths, project="TEST")

    assert paths.matrix_pca_parquet.exists()
    assert paths.pca_manifest_json.exists()

    loaded = pd.read_parquet(paths.matrix_pca_parquet)
    pd.testing.assert_frame_equal(loaded, pd.DataFrame(scores))

    manifest = json.loads(paths.pca_manifest_json.read_text())
    assert manifest["project"] == "TEST"
    assert (
        manifest["n_components_selected"]
        == scores.attrs["n_components_selected"]
    )
    assert manifest["shape"] == list(scores.shape)
    assert manifest["raw_columns"] == cov.shape[1]
    assert manifest["variance_threshold"] == 0.90
