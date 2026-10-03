"""Pooled replication timing (``pooled_mrt``)."""

import numpy as np
import pandas as pd
import pytest

from sigmutselcovs import builder
from sigmutselcovs.covariates_replication_timing import (
    pool_rt_profiles,
)
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import RepliseqSpec, TrackRef

GENES = pd.Index(
    [f"ENSG{i:011d}" for i in range(6)], name="ensembl_gene_id"
)


def test_pool_zscores_and_gives_each_biosample_one_vote():
    base = pd.Series(np.arange(6.0), index=GENES)
    # Two profiles of one biosample on another scale, one of another:
    # after z-scoring all three agree, so the pool is base's z-score.
    pooled = pool_rt_profiles(
        [("A", base * 10 + 3), ("A", base * 10 + 3), ("B", base)]
    )
    expected = (base - base.mean()) / base.std()
    assert pooled.name == "rt_pool_z"
    np.testing.assert_allclose(pooled.to_numpy(), expected.to_numpy())

    # One vote per biosample: A twice does not outweigh B once.
    flipped = pool_rt_profiles(
        [("A", base), ("A", base), ("B", -base)]
    )
    np.testing.assert_allclose(flipped.to_numpy(), 0.0, atol=1e-12)


def test_pool_uses_the_profiles_that_cover_a_gene():
    a = pd.Series([0.0, 1.0, 2.0, np.nan], index=GENES[:4])
    b = pd.Series([0.0, 1.0, 2.0, 3.0], index=GENES[:4])
    pooled = pool_rt_profiles([("A", a), ("B", b)])
    assert pooled.notna().all()
    zb = (b - b.mean()) / b.std()
    assert pooled.iloc[3] == pytest.approx(zb.iloc[3])


def test_pool_needs_a_profile():
    with pytest.raises(ValueError, match="at least one"):
        pool_rt_profiles([])


def _spec():
    return RepliseqSpec(
        type="pooled_mrt",
        assembly="hg19",
        cell_line="pool",
        profiles=(
            RepliseqSpec(
                type="wavelet",
                assembly="hg19",
                cell_line="A549",
                tracks=(TrackRef("wavelet", "ENCFF000AAA"),),
            ),
            RepliseqSpec(
                type="mat",
                assembly="hg38",
                cell_line="HCT116",
                filename="x.mat",
            ),
            RepliseqSpec(
                type="wavelet",
                assembly="hg19",
                cell_line="MISSING",
                tracks=(TrackRef("wavelet", "ENCFF000ZZZ"),),
            ),
        ),
    )


def test_builder_flips_wavelets_skips_missing_and_caches(
    tmp_path, monkeypatch
):
    paths = project_paths(tmp_path)
    paths.rt_encode_dir.mkdir(parents=True)
    (paths.rt_encode_dir / "ENCFF000AAA.bigWig").write_bytes(b"")
    (paths.rt_dir / "x.mat").write_bytes(b"")
    later = pd.Series(np.arange(6.0), index=GENES)
    seen_gtfs = []

    def fake_mrt(source, gtf, **kw):
        seen_gtfs.append(gtf)
        return later.rename("mrt")

    def fake_wavelet(path, gtf, **kw):
        seen_gtfs.append(gtf)
        return (-later).rename("rt_wavelet")  # log2 early/late

    monkeypatch.setattr(builder, "generate_mrt_per_gene", fake_mrt)
    monkeypatch.setattr(
        builder, "generate_rt_wavelet_per_gene", fake_wavelet
    )
    frames = builder._load_pooled_rt(
        _spec(), paths, lambda a: f"gtf-{a}", force_generation=False
    )
    pooled = frames[0]["rt_pool_z"]
    expected = (later - later.mean()) / later.std()
    np.testing.assert_allclose(pooled.to_numpy(), expected.to_numpy())
    assert seen_gtfs == ["gtf-hg19", "gtf-hg38"]
    assert paths.rt_pool_csv.exists()
    assert len(list(paths.rt_pool_dir.glob("*.csv"))) == 2

    # A second call reads the cache rather than recomputing.
    monkeypatch.setattr(builder, "generate_mrt_per_gene", None)
    again = builder._load_pooled_rt(
        _spec(), paths, lambda a: f"gtf-{a}", force_generation=False
    )[0]["rt_pool_z"]
    # The cache is written at 6 significant digits.
    np.testing.assert_allclose(
        again.to_numpy(), pooled.to_numpy(), rtol=1e-5
    )


def test_download_recurses_into_profiles(tmp_path, monkeypatch):
    import sigmutselcovs.download as dl

    seen = []
    real = dl.download_repliseq

    def spy(spec, paths, **kw):
        if spec.type == "pooled_mrt":
            return real(spec, paths, **kw)
        seen.append(spec.cell_line)
        return [spec.cell_line]

    monkeypatch.setattr(dl, "download_repliseq", spy)
    out = spy(_spec(), project_paths(tmp_path))
    assert seen == ["A549", "HCT116", "MISSING"]
    assert out == ["A549", "HCT116", "MISSING"]
