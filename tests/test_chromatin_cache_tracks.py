"""The chromatin cache is rebuilt when its tracks change.

A registry row that gained or swapped a track (a cohort gaining a
released mark, a treated cell-line experiment replaced by its vehicle
control) used to keep loading the cache built from the old tracks,
with no warning: the new columns never appeared.
"""

import pandas as pd

from sigmutselcovs import builder


def test_cache_follows_its_tracks(tmp_path, monkeypatch):
    calls = []

    def fake_load(covs_csv, bigwigs, gtf, force_generation, **kwargs):
        calls.append(force_generation)
        if force_generation or not covs_csv.exists():
            pd.DataFrame(
                {Path_stem(b): [1.0] for b in bigwigs},
                index=["ENSG_A"],
            ).to_csv(covs_csv)
        return pd.read_csv(covs_csv, index_col=0)

    def Path_stem(b):
        return str(b).rsplit("/", 1)[-1]

    monkeypatch.setattr(
        builder, "load_or_generate_chromatin_covariates", fake_load
    )
    cache = tmp_path / "chromatin_covs_encode.csv"
    a, b, c = (
        tmp_path / n for n in ("A.bigWig", "B.bigWig", "C.bigWig")
    )

    def load(tracks):
        return builder._load_chromatin(
            "encode", cache, tracks, "g.gtf", False
        )

    assert list(load([a, b]).columns) == ["A.bigWig", "B.bigWig"]
    assert list(load([b, a]).columns) == ["A.bigWig", "B.bigWig"]
    assert calls == [False, False]  # same tracks: the cache is used

    # A track swapped in the registry: rebuilt, not reused.
    assert list(load([a, c]).columns) == ["A.bigWig", "C.bigWig"]
    assert calls[-1] is True

    # bigWigs moved to an archive after the build: the cache stands.
    assert list(load([]).columns) == ["A.bigWig", "C.bigWig"]
