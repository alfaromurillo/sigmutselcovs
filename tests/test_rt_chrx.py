"""Replication timing keeps chrX genes, as the chromatin summaries do."""

import numpy as np

from sigmutselcovs import covariates_replication_timing as rt


def _track(path, value):
    import pyBigWig

    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([("chr1", 200_000), ("chrX", 150_000)])
    for chrom, size in (("chr1", 200_000), ("chrX", 150_000)):
        bw.addEntries([chrom], [0], ends=[size], values=[value])
    bw.close()


def test_bins_cover_chrx_by_default(tmp_path):
    paths = []
    for i, value in enumerate((1.0, 3.0)):
        paths.append(tmp_path / f"s{i}.bw")
        _track(paths[-1], value)
    bins = rt._bin_bigwigs(paths, bin_size=50_000)
    assert set(bins["Chromosome"]) == {"chr1", "chrX"}
    on_x = bins[bins["Chromosome"] == "chrX"]
    assert len(on_x) == 3
    assert np.allclose(on_x["track_2"], 3.0)


def test_default_chromosomes_are_autosomes_and_chrx():
    assert rt._RT_CHROMOSOMES[-1] == "chrX"
    assert len(rt._RT_CHROMOSOMES) == 23
