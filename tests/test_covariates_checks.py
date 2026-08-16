"""Tests for covariates_checks' variance screening.

Regression coverage for a real finding: an all-NaN column (e.g.
TGCT's tpm_unstranded_normal -- TCGA-TGCT has zero matched-normal
RNA-seq samples) was explicitly marked as *not* droppable by
check_variance's n==0 branch, contradicting its own docstring ("at
most one unique non-NaN value" trivially includes zero). Left in
the matrix, that single column poisoned every downstream
dropna(how="any") call into discarding every row, so PCA over a
975-column, 1946-gene matrix crashed with a 0-sample array.
"""

import numpy as np
import pandas as pd

from sigmutselcovs.covariates_checks import (
    check_variance,
    fix_variance,
)


def test_all_nan_column_is_flagged_constant_and_dropped():
    df = pd.DataFrame(
        {
            "all_nan": [np.nan, np.nan, np.nan, np.nan],
            "varies": [1.0, 2.0, 3.0, 4.0],
        }
    )
    report = check_variance(df)
    assert report.loc["all_nan", "is_constant"]
    assert report.loc["all_nan", "drop"]
    assert not report.loc["varies", "drop"]

    fixed, _ = fix_variance(df)
    assert "all_nan" not in fixed.columns
    assert "varies" in fixed.columns


def test_all_nan_column_does_not_poison_dropna_any():
    """The scenario that actually broke TGCT: an all-NaN column left
    in the matrix means every row has at least one NaN, so
    df.dropna(how="any") -- PCA's default missing-data handling --
    discards every single row."""
    df = pd.DataFrame(
        {
            "all_nan": [np.nan] * 5,
            "signal_a": [1.0, 2.0, 3.0, 4.0, 5.0],
            "signal_b": [5.0, 4.0, 3.0, 2.0, 1.0],
        }
    )
    fixed, _ = fix_variance(df)
    assert len(fixed.dropna(how="any")) == 5


def test_single_value_column_is_constant():
    df = pd.DataFrame({"constant": [3.0, 3.0, 3.0]})
    report = check_variance(df)
    assert report.loc["constant", "is_constant"]
    assert report.loc["constant", "drop"]


def test_varying_column_is_not_dropped():
    df = pd.DataFrame({"varies": [1.0, 5.0, 100.0, 2.0]})
    report = check_variance(df)
    assert not report.loc["varies", "is_constant"]
    assert not report.loc["varies", "drop"]
