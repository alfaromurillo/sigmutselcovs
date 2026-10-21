"""validate.py on chromatin columns built with the new options."""

import json

import numpy as np
import pandas as pd
from conftest import legacy_registry_raw

from sigmutselcovs.validate import _mark_mean, validate_covariates

GENES = [f"ENSG{i:011d}" for i in range(1, 51)]


def _frame(columns, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {c: rng.random(len(GENES)) for c in columns}, index=GENES
    )


def test_mark_mean_prefers_linear_then_log():
    both = _frame(
        [
            "e075_h3k4me3_fc_signal_body",
            "e075_h3k4me3_fc_signal_body_log",
        ]
    )
    pd.testing.assert_series_equal(
        _mark_mean(both, "h3k4me3", "body"),
        both[["e075_h3k4me3_fc_signal_body"]].mean(axis=1),
    )
    log_only = both[["e075_h3k4me3_fc_signal_body_log"]]
    pd.testing.assert_series_equal(
        _mark_mean(log_only, "h3k4me3", "body"),
        log_only.mean(axis=1),
    )
    peaks_only = _frame(["e075_h3k4me3_narrowpeak_body"])
    assert _mark_mean(peaks_only, "h3k4me3", "body") is None


def _check(report, name):
    return report.set_index("check").loc[name]


def test_collapsed_column_count_counts_new_summaries(tmp_path):
    columns = [
        f"{mark}{summary}"
        for mark in ("h3k4me3", "h3k9me3")
        for summary in (
            "_body_log",
            "_cds_span_log",
            "_narrowpeak_body",
            "_narrowpeak_cds_span",
        )
    ]
    raw = legacy_registry_raw()
    raw["projects"]["GENERIC"]["roadmap"].update(
        windows=["body", "cds_span"],
        transform="log1p",
        peaks=["narrowPeak"],
    )
    registry = tmp_path / "projects.json"
    registry.write_text(json.dumps(raw))
    report = validate_covariates(
        "GENERIC",
        cov_matrix_raw=_frame(columns),
        registry_path=registry,
    )
    row = _check(report, "chromatin_collapsed_column_count")
    assert row["status"] == "pass"
    assert row["value"] == len(columns)


def test_roadmap_column_count_scales_with_windows(tmp_path):
    columns = [
        f"e075_h3k4me3_fc_signal_{w}"
        for w in ("body", "promoter", "cds_span")
    ]
    raw = legacy_registry_raw()
    raw["projects"]["COAD"]["roadmap"].update(
        eids=["E075"], marks=["H3K4me3"]
    )
    default = tmp_path / "default.json"
    default.write_text(json.dumps(raw))
    raw["projects"]["COAD"]["roadmap"]["windows"] = [
        "body",
        "promoter",
        "cds_span",
    ]
    three = tmp_path / "three.json"
    three.write_text(json.dumps(raw))
    frame = _frame(columns)
    for registry, status in ((default, "fail"), (three, "pass")):
        report = validate_covariates(
            "COAD", cov_matrix_raw=frame, registry_path=registry
        )
        assert (
            _check(report, "roadmap_column_count")["status"] == status
        )
