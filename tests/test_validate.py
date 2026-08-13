"""Validation tests: synthetic matrices with planted violations."""

import numpy as np
import pandas as pd
import pytest

from sigmutselcovs.validate import (
    HOUSEKEEPING_GENES,
    validate_covariates,
)

RNG = np.random.default_rng(7)
N = 400


def _index():
    ids = [f"ENSG{i:011d}" for i in range(1, N + 1)]
    # plant the housekeeping genes into the index
    for k, gid in enumerate(HOUSEKEEPING_GENES.values()):
        ids[k] = gid
    return pd.Index(ids, name="ensembl_gene_id")


def _healthy_matrix() -> pd.DataFrame:
    """A biologically sensible synthetic COAD-shaped matrix."""
    index = _index()
    # latent "activity" drives everything with the right signs
    activity = RNG.normal(size=N)
    noise = lambda: RNG.normal(scale=0.4, size=N)  # noqa: E731

    def expr_like(shift=0.0):
        return np.exp(activity + shift + noise())

    df = pd.DataFrame(index=index)
    df["gtex_colon_sigmoid"] = expr_like()
    df["gtex_colon_transverse_mucosa"] = expr_like()
    df["tpm_unstranded"] = expr_like()
    df["tpm_unstranded_normal"] = expr_like()
    df["TCGA-AA-0001-01A_tpm_unstranded"] = expr_like()
    mrt = 1 / (1 + np.exp(activity + noise()))  # active -> early
    df["mrt"] = mrt
    for i in range(1, 17):
        df[f"clr_rt_s{i}"] = RNG.normal(size=N)
    df["coad_s1_insertions_body"] = expr_like(0.5)
    df["coad_s1_insertions_promoter"] = expr_like(0.2)
    for eid in ("e075", "e106"):
        for mark, sign in (
            ("h3k4me3", 1),
            ("h3k27ac", 1),
            ("h3k9ac", 1),
            ("h3k9me3", -1),
            ("h3k27me3", -1),
            ("h3k36me3", 1),
        ):
            df[f"{eid}_{mark}_fc_signal_body"] = np.exp(
                sign * activity + noise()
            )
            df[f"{eid}_{mark}_fc_signal_promoter"] = np.exp(
                sign * activity + noise()
            )
    # make housekeeping genes highly expressed and early
    for gid in HOUSEKEEPING_GENES.values():
        for col in (
            "gtex_colon_sigmoid",
            "gtex_colon_transverse_mucosa",
            "tpm_unstranded",
        ):
            df.loc[gid, col] = df[col].max() * 2
        df.loc[gid, "mrt"] = 0.01
    return df


def _statuses(frame: pd.DataFrame) -> dict[str, str]:
    return dict(zip(frame["check"], frame["status"]))


def test_healthy_matrix_passes():
    frame = validate_covariates(
        "COAD", cov_matrix_raw=_healthy_matrix()
    )
    assert not (frame["status"] == "fail").any()
    statuses = _statuses(frame)
    assert statuses["expression_vs_mrt"] == "pass"
    assert statuses["h3k4me3_vs_expression"] == "pass"
    assert statuses["h3k9me3_vs_expression"] == "pass"
    assert statuses["atac_vs_expression"] == "pass"
    assert statuses["h3k4me3_eid_coherence"] == "pass"
    assert statuses["housekeeping_panel"] == "pass"
    assert statuses["mrt_range"] == "pass"
    assert statuses["index_unique"] == "pass"


def test_almost_all_nan_column_fails():
    df = _healthy_matrix()
    df["e075_h3k4me3_fc_signal_body"] = np.nan
    frame = validate_covariates("COAD", cov_matrix_raw=df)
    assert _statuses(frame)["nan_fraction"] == "fail"
    row = frame.set_index("check").loc["nan_fraction"]
    assert "e075_h3k4me3_fc_signal_body" in row["note"]
    with pytest.raises(ValueError, match="nan_fraction"):
        validate_covariates("COAD", cov_matrix_raw=df, strict=True)


def test_constant_column_warns():
    df = _healthy_matrix()
    df["gtex_colon_sigmoid"] = 1.0
    frame = validate_covariates("COAD", cov_matrix_raw=df)
    assert _statuses(frame)["constant_columns"] == "warn"


def test_mrt_out_of_range_fails():
    df = _healthy_matrix()
    df.loc[df.index[10], "mrt"] = 1.7
    frame = validate_covariates("COAD", cov_matrix_raw=df)
    assert _statuses(frame)["mrt_range"] == "fail"


def test_inverted_biology_warns_and_batches():
    df = _healthy_matrix()
    # invert expression so every cross-block direction flips
    for col in (
        "gtex_colon_sigmoid",
        "gtex_colon_transverse_mucosa",
        "tpm_unstranded",
    ):
        df[col] = df[col].max() - df[col]
    frame = validate_covariates("COAD", cov_matrix_raw=df)
    statuses = _statuses(frame)
    assert statuses["expression_vs_mrt"] == "warn"
    assert statuses["biology_batch"] == "fail"


def test_duplicate_index_fails():
    df = _healthy_matrix()
    df.index = [df.index[0]] + list(df.index[1:])
    df = pd.concat([df, df.iloc[[0]]])
    frame = validate_covariates("COAD", cov_matrix_raw=df)
    assert _statuses(frame)["index_unique"] == "fail"


def test_clr_count_mismatch_fails():
    df = _healthy_matrix().drop(columns=["clr_rt_s16"])
    frame = validate_covariates("COAD", cov_matrix_raw=df)
    assert _statuses(frame)["clr_fraction_count"] == "fail"


def test_missing_blocks_not_applicable():
    df = _healthy_matrix()[
        ["gtex_colon_sigmoid", "gtex_colon_transverse_mucosa"]
    ]
    frame = validate_covariates("COAD", cov_matrix_raw=df)
    row = frame.set_index("check").loc["expression_vs_mrt"]
    assert row["status"] == "pass"
    assert "not applicable" in row["note"]


def test_requires_matrix_or_data_dir():
    with pytest.raises(
        ValueError, match="data_dir or cov_matrix_raw"
    ):
        validate_covariates("COAD")
