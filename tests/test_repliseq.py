"""Replication-timing loader tests (synthetic data, no network)."""

import numpy as np
import pandas as pd
import pytest

from sigmutselcovs.covariates_replication_timing import (
    load_repliseq_fractions_bins,
    load_repliseq_fractions_bins_from_bigwigs,
    load_repliseq_mrt_bins,
)

BIN = 50_000


def _write_bigwig(path, values, *, chrom="chr1", bin_size=BIN):
    """One constant value per bin on a single chromosome."""
    import pyBigWig

    length = bin_size * len(values)
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, length)])
    starts = list(range(0, length, bin_size))
    ends = [s + bin_size for s in starts]
    bw.addEntries(
        [chrom] * len(starts),
        starts,
        ends=ends,
        values=[float(v) for v in values],
    )
    bw.close()
    return path


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
            [1.0, 0.0, 2.0],  # s1 per bin
            [1.0, 0.0, 0.0],  # s2
            [1.0, 0.0, 0.0],  # s3
            [1.0, 0.0, 0.0],  # s4
        ],
    )
    return path


def test_fractions_bins_normalization_pin(mat_file):
    """Pins the pre-refactor behavior of the normalization."""
    out = load_repliseq_fractions_bins(mat_file)
    assert list(out.columns) == [
        "Chromosome",
        "region_start",
        "region_end",
        "rt_s1",
        "rt_s2",
        "rt_s3",
        "rt_s4",
    ]
    rt = out[[f"rt_s{i}" for i in range(1, 5)]].to_numpy()
    np.testing.assert_allclose(rt[0], [0.25, 0.25, 0.25, 0.25])
    assert np.isnan(rt[1]).all()  # zero-signal bin -> NaN
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


def test_bigwig_adapter_two_fractions(tmp_path):
    early = _write_bigwig(tmp_path / "early.bigWig", [1, 0, 3, 0])
    late = _write_bigwig(tmp_path / "late.bigWig", [1, 2, 1, 0])
    out = load_repliseq_fractions_bins_from_bigwigs([early, late])
    assert list(out.columns) == [
        "Chromosome",
        "region_start",
        "region_end",
        "rt_s1",
        "rt_s2",
    ]
    assert len(out) == 4
    assert out["region_start"].tolist() == [0, BIN, 2 * BIN, 3 * BIN]
    np.testing.assert_allclose(
        out["rt_s1"].to_numpy()[:3], [0.5, 0.0, 0.75]
    )
    np.testing.assert_allclose(
        out["rt_s2"].to_numpy()[:3], [0.5, 1.0, 0.25]
    )
    assert pd.isna(out["rt_s1"].iloc[3])  # all-zero bin


def test_bigwig_adapter_six_fractions_sum_to_one(tmp_path):
    paths = [
        _write_bigwig(tmp_path / f"f{i}.bigWig", [i + 1, 2 * i + 1])
        for i in range(6)
    ]
    out = load_repliseq_fractions_bins_from_bigwigs(paths)
    rt = out[[f"rt_s{i}" for i in range(1, 7)]].to_numpy()
    np.testing.assert_allclose(rt.sum(axis=1), [1.0, 1.0])


def test_bigwig_adapter_chrom_name_reconciliation(tmp_path):
    # one track uses '1', the other 'chr1'
    a = _write_bigwig(tmp_path / "a.bigWig", [2, 2], chrom="1")
    b = _write_bigwig(tmp_path / "b.bigWig", [2, 6], chrom="chr1")
    out = load_repliseq_fractions_bins_from_bigwigs([a, b])
    assert out["Chromosome"].unique().tolist() == ["chr1"]
    np.testing.assert_allclose(out["rt_s1"].to_numpy(), [0.5, 0.25])


def test_bigwig_adapter_rejects_single_track(tmp_path):
    a = _write_bigwig(tmp_path / "a.bigWig", [1])
    with pytest.raises(ValueError, match="at least two"):
        load_repliseq_fractions_bins_from_bigwigs([a])


def test_mrt_two_fractions_closed_form(tmp_path):
    """N=2: mrt = 0.25 + 0.5 * f2 (midpoints 0.25 and 0.75)."""
    early = _write_bigwig(tmp_path / "e.bigWig", [1, 3, 0])
    late = _write_bigwig(tmp_path / "l.bigWig", [1, 1, 2])
    out = load_repliseq_mrt_bins([early, late])
    f2 = np.array([0.5, 0.25, 1.0])
    np.testing.assert_allclose(out["mrt"].to_numpy(), 0.25 + 0.5 * f2)


def test_mrt_six_fractions_midpoints(tmp_path):
    """N=6 midpoints [1,3,5,7,9,11]/12 (UW weighted average)."""
    paths = [
        _write_bigwig(
            tmp_path / f"f{i}.bigWig", [1 if i == 2 else 0, 1]
        )
        for i in range(6)
    ]
    out = load_repliseq_mrt_bins(paths)
    # all signal in fraction 3 -> its midpoint 5/12
    assert out["mrt"].iloc[0] == pytest.approx(5 / 12)
    # uniform -> 0.5
    assert out["mrt"].iloc[1] == pytest.approx(0.5)


def test_mrt_fraction_subset_renormalizes(tmp_path):
    paths = [
        _write_bigwig(tmp_path / f"f{i}.bigWig", [v])
        for i, v in enumerate([5.0, 1.0, 3.0, 2.0])
    ]
    out = load_repliseq_mrt_bins(
        paths, mrt_fraction_cols=["rt_s2", "rt_s3"]
    )
    # subset renormalized: f = [1,3]/4; midpoints [0.25, 0.75]
    expected = (0.25 * 1 + 0.75 * 3) / 4
    assert out["mrt"].iloc[0] == pytest.approx(expected)


def test_mrt_fraction_subset_unknown_column(mat_file):
    with pytest.raises(ValueError, match="mrt_fraction_cols"):
        load_repliseq_mrt_bins(mat_file, mrt_fraction_cols=["rt_s99"])


def test_explicit_source_type_overrides(tmp_path, mat_file):
    out = load_repliseq_fractions_bins(mat_file, source_type="mat")
    assert len(out) == 3
    with pytest.raises(ValueError, match="wavelet"):
        load_repliseq_fractions_bins(mat_file, source_type="wavelet")
    with pytest.raises(ValueError, match="source_type"):
        load_repliseq_fractions_bins(mat_file, source_type="bogus")


def _write_gtf(path, genes):
    """genes: list of (gene_id, chrom, start1, end1) 1-based incl."""
    lines = []
    for gid, chrom, start, end in genes:
        attrs = f'gene_id "{gid}"; gene_type "protein_coding";'
        lines.append(
            "\t".join(
                [
                    chrom,
                    "TEST",
                    "gene",
                    str(start),
                    str(end),
                    ".",
                    "+",
                    ".",
                    attrs,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def test_wavelet_per_gene(tmp_path):
    from sigmutselcovs.covariates_replication_timing import (
        generate_rt_wavelet_per_gene,
        load_or_generate_rt_wavelet,
    )

    track = _write_bigwig(tmp_path / "wavelet.bigWig", [2, 4, 6])
    gtf = _write_gtf(
        tmp_path / "genes.gtf",
        [
            # spans bins 1-2 fully -> mean(2, 4) = 3
            ("ENSG00000000001", "chr1", 1, 2 * BIN),
            # inside bin 3 only -> 6
            ("ENSG00000000002", "chr1", 2 * BIN + 1, 2 * BIN + 1000),
        ],
    )
    ser = generate_rt_wavelet_per_gene(track, gtf)
    assert ser.name == "rt_wavelet"
    assert ser.loc["ENSG00000000001"] == pytest.approx(3.0)
    assert ser.loc["ENSG00000000002"] == pytest.approx(6.0)

    cache = tmp_path / "rt_wavelet_per_gene.csv"
    first = load_or_generate_rt_wavelet(cache, track, gtf)
    assert cache.exists()
    again = load_or_generate_rt_wavelet(
        cache, "nonexistent.bigWig", gtf
    )
    pd.testing.assert_series_equal(first, again)
