"""The ``log1p_mean`` transform: log1p of a window's linear mean."""

import numpy as np
import pandas as pd
import pytest
from test_chromatin_windows import WINDOWS, _assert_close, _reference

from sigmutselcovs import covariates_chromatin as cc
from sigmutselcovs.covariates_checks import fix_all


@pytest.fixture
def genes(synthetic_gtf):
    return cc._load_gene_bodies(
        synthetic_gtf, ("protein_coding",), with_cds_span=True
    )


def test_is_log1p_of_the_clipped_linear_mean(synthetic_bigwig, genes):
    path, dense = synthetic_bigwig
    out = cc.summarize_bigwig_windows(
        path,
        genes,
        windows=WINDOWS,
        transform="log1p_mean",
        label="trk",
    )
    assert list(out.columns) == [f"trk_{w}_logmean" for w in WINDOWS]
    for window in WINDOWS:
        coords = cc.gene_windows(genes, window)
        expected = [
            np.log1p(
                max(
                    _reference(
                        dense,
                        row["Chromosome"],
                        row["start"],
                        row["end"],
                        None,
                    ),
                    0,
                )
            )
            for _, row in coords.iterrows()
        ]
        _assert_close(out[f"trk_{window}_logmean"], expected)


def test_differs_from_the_per_base_log(synthetic_bigwig, genes):
    path, _ = synthetic_bigwig
    kw = {"windows": ["cds_span"], "label": "trk"}
    per_base = cc.summarize_bigwig_windows(
        path, genes, transform="log1p", **kw
    )["trk_cds_span_log"]
    of_mean = cc.summarize_bigwig_windows(
        path, genes, transform="log1p_mean", **kw
    )["trk_cds_span_logmean"]
    both = per_base.notna() & of_mean.notna()
    # Jensen: the log of the mean is at least the mean of the logs
    assert (of_mean[both] >= per_base[both] - 1e-12).all()
    assert (of_mean[both] > per_base[both]).any()


def test_suffixes_round_trip():
    for transform in (None, "log1p", "log1p_mean"):
        name = cc.window_column("E075-H3K4me3", "cds_span", transform)
        assert cc.split_transform_suffix(name) == (
            "e075_h3k4me3_cds_span",
            transform,
        )
    # a label that merely ends in "log" is not a transform suffix
    assert cc.split_transform_suffix("dialog") == ("dialog", None)


def test_collapse_keeps_the_two_logs_apart():
    df = pd.DataFrame(
        {
            "e1_h3k4me3_cds_span_log": [1.0],
            "e2_h3k4me3_cds_span_log": [3.0],
            "e1_h3k4me3_cds_span_logmean": [10.0],
            "e2_h3k4me3_cds_span_logmean": [30.0],
        }
    )
    out = cc._collapse_by_assay(df)
    assert out.loc[0, "h3k4me3_cds_span_log"] == 2.0
    assert out.loc[0, "h3k4me3_cds_span_logmean"] == 20.0


def test_cache_path_names_the_transform(tmp_path):
    base = tmp_path / "chromatin_covs_roadmap.csv"
    path = cc.chromatin_cache_path(
        base, windows=["cds_span"], transform="log1p_mean"
    )
    assert (
        path.name == "chromatin_covs_roadmap_cds_span_log1p_mean.csv"
    )


def test_fix_all_leaves_skipped_columns_alone():
    rng = np.random.default_rng(0)
    skewed = rng.lognormal(0, 2, 500)
    df = pd.DataFrame({"a_cds_span_logmean": skewed, "b": skewed})
    fixed, reports = fix_all(df, skip_skewness=["a_cds_span_logmean"])
    assert np.array_equal(fixed["a_cds_span_logmean"], skewed)
    assert not np.array_equal(fixed["b"], skewed)
    assert "a_cds_span_logmean" not in reports["skewness"].index


def test_build_does_not_log_the_columns_again(tmp_path):
    import json

    from conftest import legacy_registry_raw, write_bigwig, write_gtf

    from sigmutselcovs.builder import build_covariate_matrix
    from sigmutselcovs.paths import project_paths

    paths = project_paths(tmp_path / "data")
    paths.roadmap_dir.mkdir(parents=True)
    for seed, track in enumerate(("E075-H3K4me3", "E106-H3K4me3")):
        write_bigwig(
            paths.roadmap_dir / f"{track}.fc.signal.bigwig", seed
        )
    gtf = write_gtf(tmp_path / "g.gtf")
    raw = legacy_registry_raw()
    raw["projects"]["COAD"]["roadmap"].update(
        windows=["cds_span"], transform="log1p_mean"
    )
    reg = tmp_path / "projects.json"
    reg.write_text(json.dumps(raw))
    kwargs = {
        "registry_path": reg,
        "gencode_gtfs": {"hg19": gtf, "hg38": gtf},
        "include": ["roadmap"],
        "validate_output": False,
    }
    fixed = build_covariate_matrix("COAD", paths.root, **kwargs).full
    raw_full = build_covariate_matrix(
        "COAD", paths.root, apply_fixes=False, **kwargs
    ).full
    assert list(fixed.columns) == [
        "e075_h3k4me3_fc_signal_cds_span_logmean",
        "e106_h3k4me3_fc_signal_cds_span_logmean",
    ]
    pd.testing.assert_frame_equal(fixed, raw_full)
