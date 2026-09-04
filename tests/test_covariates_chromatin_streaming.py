"""Tests for the streaming chromatin-covariate path (added 2026-09-04
for building large track pools, e.g. GENERIC, without needing every
bigwig on disk simultaneously -- see covariates_chromatin.py's
summarize_tracks_to_genes_streaming/
load_or_generate_chromatin_covariates_streaming and download.py's
stream_roadmap_tracks/stream_encode_chromatin_tracks)."""

import pandas as pd
import pytest

from sigmutselcovs import covariates_chromatin as cc
from sigmutselcovs import download as dl
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import (
    EncodeChromatinSpec,
    RoadmapSpec,
    TrackRef,
)

GENES = pd.DataFrame(
    {
        "Chromosome": ["chr1", "chr1"],
        "start": [0, 1000],
        "end": [500, 1500],
        "strand": ["+", "+"],
    },
    index=pd.Index(["ENSG_A", "ENSG_B"], name="ensembl_gene_id"),
)


def _fake_summaries(monkeypatch, calls):
    """Stub the two low-level bigwig readers so these tests never
    touch a real bigwig file -- each call is recorded in `calls` so
    tests can assert on ordering/disk state at the moment it ran."""

    def fake_body(path, genes, *, label, statistic="mean"):
        calls.append(("body", str(path)))
        return pd.Series([1.0, 2.0], index=genes.index, name=label)

    def fake_promoter(
        path, genes, *, upstream, downstream, label, statistic="mean"
    ):
        calls.append(("promoter", str(path)))
        return pd.Series([3.0, 4.0], index=genes.index, name=label)

    monkeypatch.setattr(cc, "summarize_bigwig_over_genes", fake_body)
    monkeypatch.setattr(
        cc, "summarize_promoter_signal", fake_promoter
    )


# --- covariates_chromatin.summarize_tracks_to_genes_streaming ---


def test_streaming_matches_batch_output(monkeypatch, tmp_path):
    calls = []
    _fake_summaries(monkeypatch, calls)

    paths = [tmp_path / "a.bigwig", tmp_path / "b.bigwig"]
    for p in paths:
        p.touch()

    batch = cc.summarise_tracks_to_genes(
        [str(p) for p in paths], GENES
    )
    streaming = cc.summarize_tracks_to_genes_streaming(
        ((p.stem, p) for p in paths), GENES
    )

    pd.testing.assert_frame_equal(
        batch.sort_index(axis=1), streaming.sort_index(axis=1)
    )


def test_streaming_consumes_one_track_at_a_time(
    monkeypatch, tmp_path
):
    """The whole point of streaming: a generator that deletes its
    file after yielding must see that file gone by the time the next
    track's summary call happens -- i.e. summarize_tracks_to_genes_
    streaming must not pre-fetch/hold onto more than one track."""
    calls = []
    _fake_summaries(monkeypatch, calls)

    seen_disk_state = []

    def track_source():
        for i in range(3):
            p = tmp_path / f"track{i}.bigwig"
            p.write_bytes(b"x")
            yield f"t{i}", p
            seen_disk_state.append(
                sorted(f.name for f in tmp_path.glob("*.bigwig"))
            )
            p.unlink()

    cc.summarize_tracks_to_genes_streaming(track_source(), GENES)

    # After each track is summarized (both body+promoter calls made)
    # and control returns to the generator, exactly that one file
    # was on disk -- never more than one at a time.
    assert seen_disk_state == [
        ["track0.bigwig"],
        ["track1.bigwig"],
        ["track2.bigwig"],
    ]
    assert not list(tmp_path.glob("*.bigwig"))  # all cleaned up


def test_streaming_no_promoter(monkeypatch, tmp_path):
    calls = []
    _fake_summaries(monkeypatch, calls)
    p = tmp_path / "a.bigwig"
    p.touch()

    out = cc.summarize_tracks_to_genes_streaming(
        [("a", p)], GENES, include_promoter=False
    )
    assert list(out.columns) == ["a_body"]
    assert all(kind == "body" for kind, _ in calls)


# --- covariates_chromatin.load_or_generate_chromatin_covariates_streaming ---


def test_streaming_cache_hit_never_consumes_track_source(tmp_path):
    location_df = tmp_path / "cached.csv"
    pd.DataFrame({"x_body": [1.0, 2.0]}, index=GENES.index).to_csv(
        location_df
    )

    def poison():
        raise AssertionError(
            "track_source must not be consumed on a cache hit"
        )
        yield  # pragma: no cover - unreachable, makes this a generator

    out = cc.load_or_generate_chromatin_covariates_streaming(
        location_df, poison(), "unused.gtf"
    )
    assert list(out.columns) == ["x_body"]


def test_streaming_average_by_assay(monkeypatch, tmp_path):
    monkeypatch.setattr(
        cc, "_load_gene_bodies", lambda gtf_path, biotypes: GENES
    )
    calls = []
    _fake_summaries(monkeypatch, calls)
    # Labels here are already in sanitize_feature_label's form
    # (lowercase, underscore-separated) since the stubbed summarizers
    # in _fake_summaries bypass real sanitization -- real callers get
    # this for free from summarize_bigwig_over_genes/
    # summarize_promoter_signal.
    paths = [
        tmp_path / "e1_h3k27ac.bigwig",
        tmp_path / "e2_h3k27ac.bigwig",
    ]
    for p in paths:
        p.touch()

    out = cc.load_or_generate_chromatin_covariates_streaming(
        tmp_path / "out.csv",
        ((p.stem, p) for p in paths),
        "unused.gtf",
        average_by_assay=True,
    )
    assert list(out.columns) == ["h3k27ac_body", "h3k27ac_promoter"]


# --- download.stream_roadmap_tracks ---


class _FakeResponse:
    def __init__(self, content: bytes, status_code: int = 200):
        self.content = content
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise dl.requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size):
        yield self.content

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class _RoadmapSession:
    def __init__(self, missing=()):
        self.missing = set(missing)
        self.requests = []

    def get(self, url, stream=True, headers=None, timeout=None):
        self.requests.append(url)
        if any(m in url for m in self.missing):
            return _FakeResponse(b"", status_code=404)
        return _FakeResponse(b"bigwigbytes")


def test_stream_roadmap_tracks_deletes_after_each_yield(tmp_path):
    spec = RoadmapSpec(eids=("E001", "E002"), marks=("H3K4me1",))
    paths = project_paths(tmp_path)
    session = _RoadmapSession()

    seen = []
    for label, path in dl.stream_roadmap_tracks(
        spec, paths, session=session
    ):
        assert path.exists()
        seen.append((label, path.exists()))
    # after the generator loop fully exits, nothing is left on disk
    assert not list(paths.roadmap_dir.glob("*.bigwig"))
    assert [label for label, _ in seen] == [
        "E001-H3K4me1",
        "E002-H3K4me1",
    ]


def test_stream_roadmap_tracks_keep_files(tmp_path):
    spec = RoadmapSpec(eids=("E001",), marks=("H3K4me1",))
    paths = project_paths(tmp_path)
    session = _RoadmapSession()

    list(
        dl.stream_roadmap_tracks(
            spec, paths, session=session, keep_files=True
        )
    )
    assert (
        paths.roadmap_dir / "E001-H3K4me1.fc.signal.bigwig"
    ).exists()


def test_stream_roadmap_tracks_tolerates_missing_marks(tmp_path):
    spec = RoadmapSpec(eids=("E001",), marks=("H3K4me1", "H3K27ac"))
    paths = project_paths(tmp_path)
    session = _RoadmapSession(missing=("H3K27ac",))

    out = list(dl.stream_roadmap_tracks(spec, paths, session=session))
    assert [label for label, _ in out] == ["E001-H3K4me1"]


def test_stream_roadmap_tracks_raises_on_required_mark(tmp_path):
    spec = RoadmapSpec(
        eids=("E001",),
        marks=("H3K27ac",),
        required_marks=("H3K27ac",),
    )
    paths = project_paths(tmp_path)
    session = _RoadmapSession(missing=("H3K27ac",))

    with pytest.raises(dl.requests.HTTPError):
        list(dl.stream_roadmap_tracks(spec, paths, session=session))


# --- download.stream_encode_chromatin_tracks ---


def test_stream_encode_chromatin_tracks_deletes_after_each_yield(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        dl,
        "resolve_encode_file",
        lambda acc, session=None: {
            "accession": acc,
            "url": f"http://s3/{acc}",
            "md5sum": None,
            "file_size": None,
            "assembly": "GRCh38",
        },
    )
    spec = EncodeChromatinSpec(
        tracks=(
            TrackRef(label="DNase", accession="ENCFF000AAA"),
            TrackRef(label="H3K4me1", accession="ENCFF000BBB"),
        )
    )
    paths = project_paths(tmp_path)
    session = _RoadmapSession()  # generic 200-OK-everything session

    out = list(
        dl.stream_encode_chromatin_tracks(
            spec, paths, session=session
        )
    )
    assert [label for label, _ in out] == [
        "DNase_fc_signal_ENCFF000AAA",
        "H3K4me1_fc_signal_ENCFF000BBB",
    ]
    assert not list(paths.encode_chromatin_dir.glob("*.bigWig"))


def test_stream_encode_chromatin_tracks_unique_labels_for_duplicate_marks(
    tmp_path, monkeypatch
):
    """The actual bug this guards: GENERIC pools many tracks that
    all share the same TrackRef.label (e.g. 83 tracks all labeled
    "DNase") -- the yielded label must still be unique per track
    (accession-qualified), or summarize_tracks_to_genes_streaming's
    dict-keyed column assembly silently collapses them to one
    column."""
    monkeypatch.setattr(
        dl,
        "resolve_encode_file",
        lambda acc, session=None: {
            "accession": acc,
            "url": f"http://s3/{acc}",
            "md5sum": None,
            "file_size": None,
            "assembly": "GRCh38",
        },
    )
    spec = EncodeChromatinSpec(
        tracks=(
            TrackRef(label="DNase", accession="ENCFF000AAA"),
            TrackRef(label="DNase", accession="ENCFF000BBB"),
        )
    )
    paths = project_paths(tmp_path)
    session = _RoadmapSession()

    out = list(
        dl.stream_encode_chromatin_tracks(
            spec, paths, session=session
        )
    )
    labels = [label for label, _ in out]
    assert len(labels) == len(set(labels)) == 2
