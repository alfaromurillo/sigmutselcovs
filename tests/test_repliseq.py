"""Replication-timing loader tests (synthetic data, no network)."""

import numpy as np
import pandas as pd
import pytest

from sigmutselcovs.covariates_replication_timing import (
    load_repliseq_fractions_bins,
    load_repliseq_mrt_bins,
)


def _write_mat(path, chroms, starts, ends, fractions):
    """Write a transposed Repli-seq table (features as rows)."""
    rows = [chroms, starts, ends, *fractions]
    lines = ["\t".join(str(x) for x in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def mat_file(tmp_path):
    """3 bins x 4 fractions: uniform, empty, all-early."""
    path = tmp_path / "synthetic.mat"
    _write_mat(
        path,
        chroms=["chr1", "chr1", "chr2"],
        starts=[0, 50000, 0],
        ends=[50000, 100000, 50000],
        fractions=[
            [1.0, 0.0, 2.0],   # s1 per bin
            [1.0, 0.0, 0.0],   # s2
            [1.0, 0.0, 0.0],   # s3
            [1.0, 0.0, 0.0],   # s4
        ])
    return path


def test_fractions_bins_normalization_pin(mat_file):
    """Pins the pre-refactor behavior of the normalization."""
    out = load_repliseq_fractions_bins(mat_file)
    assert list(out.columns) == ["Chromosome", "region_start",
                                 "region_end", "rt_s1", "rt_s2",
                                 "rt_s3", "rt_s4"]
    rt = out[[f"rt_s{i}" for i in range(1, 5)]].to_numpy()
    np.testing.assert_allclose(rt[0], [0.25, 0.25, 0.25, 0.25])
    assert np.isnan(rt[1]).all()          # zero-signal bin -> NaN
    np.testing.assert_allclose(rt[2], [1.0, 0.0, 0.0, 0.0])
    assert out["region_start"].tolist() == [0, 50000, 0]
    assert out["region_end"].tolist() == [50000, 100000, 50000]


def test_mrt_bins_midpoints(mat_file):
    out = load_repliseq_mrt_bins(mat_file)
    # uniform bin: mean of midpoints (i+0.5)/4 = 0.5
    assert out["mrt"].iloc[0] == pytest.approx(0.5)
    assert pd.isna(out["mrt"].iloc[1])
    # all-early bin: first midpoint 0.125
    assert out["mrt"].iloc[2] == pytest.approx(0.125)
