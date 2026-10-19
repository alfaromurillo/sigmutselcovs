"""Peak-file downloaders (Roadmap and ENCODE), with a fake session."""

import gzip
import hashlib

import pytest
import requests

import sigmutselcovs.download as dl
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import (
    EncodeChromatinSpec,
    PeakRef,
    RoadmapSpec,
    TrackRef,
)

PEAKS = gzip.compress(b"chr1\t100\t200\tp\t0\t.\t1\t1\t1\t50\n")


class _Response:
    def __init__(self, payload, status_code=200):
        self.payload = payload
        self.status_code = status_code

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size):
        yield self.payload


class _Session:
    """Serves a url -> bytes mapping; other urls are 404."""

    def __init__(self, mapping):
        self.mapping = mapping
        self.fetched = []

    def get(self, url, stream=True, headers=None, timeout=None):
        self.fetched.append(url)
        if url not in self.mapping:
            return _Response(b"", status_code=404)
        return _Response(self.mapping[url])


def _roadmap(**kwargs):
    return RoadmapSpec(
        eids=("E075", "E106"),
        marks=("H3K4me3", "H3K27ac"),
        peaks=("narrowPeak", "broadPeak"),
        **kwargs,
    )


def test_roadmap_peaks_url_layout():
    spec = _roadmap()
    assert spec.peaks_url_template.format(
        eid="E106", mark="H3K27ac", peak_type="broadPeak"
    ) == (
        "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/"
        "consolidated/broadPeak/E106-H3K27ac.broadPeak.gz"
    )


def test_download_roadmap_peaks(tmp_path):
    spec = _roadmap(
        peaks_url_template="http://r/{peak_type}/{eid}-{mark}"
    )
    # E106-H3K27ac has no peaks: skipped, not fatal
    mapping = {
        f"http://r/{t}/{e}-{m}": PEAKS
        for e in spec.eids
        for m in spec.marks
        for t in spec.peaks
        if (e, m) != ("E106", "H3K27ac")
    }
    paths = project_paths(tmp_path)
    session = _Session(mapping)
    got = dl.download_roadmap_peaks(spec, paths, session=session)
    assert [p.name for p in got] == [
        "E075-H3K4me3.narrowPeak.gz",
        "E075-H3K4me3.broadPeak.gz",
        "E075-H3K27ac.narrowPeak.gz",
        "E075-H3K27ac.broadPeak.gz",
        "E106-H3K4me3.narrowPeak.gz",
        "E106-H3K4me3.broadPeak.gz",
    ]
    assert all(p.parent == paths.roadmap_peaks_dir for p in got)
    assert all(p.read_bytes() == PEAKS for p in got)  # kept gzipped
    assert not list(paths.roadmap_peaks_dir.glob("*.part"))
    # listed by peak type, then by name
    assert paths.roadmap_peak_files(["narrowPeak"]) == sorted(
        p for p in got if "narrowPeak" in p.name
    )

    # idempotent: present files are not fetched again
    session.fetched.clear()
    again = dl.download_roadmap_peaks(spec, paths, session=session)
    assert again == got
    assert session.fetched == [
        "http://r/narrowPeak/E106-H3K27ac",
        "http://r/broadPeak/E106-H3K27ac",
    ]


def test_download_roadmap_peaks_required_mark(tmp_path):
    spec = _roadmap(
        peaks_url_template="http://r/{peak_type}/{eid}-{mark}",
        required_marks=("H3K27ac",),
    )
    with pytest.raises(requests.HTTPError, match="404"):
        dl.download_roadmap_peaks(
            spec, project_paths(tmp_path), session=_Session({})
        )


def test_download_encode_chromatin_peaks(tmp_path, monkeypatch):
    spec = EncodeChromatinSpec(
        tracks=(TrackRef(label="H3K4me3", accession="ENCFF000SIG"),),
        peak_tracks=(
            PeakRef(label="H3K4me3", accession="ENCFF000PKA"),
            PeakRef(
                label="H3K27ac",
                accession="ENCFF000PKB",
                peak_type="broadPeak",
                url="http://direct/pkb",
            ),
        ),
    )

    def resolve(accession, session=None):
        return {
            "url": f"http://s3/{accession}",
            "md5sum": hashlib.md5(PEAKS).hexdigest(),
            "file_size": len(PEAKS),
            "file_format": "bed",
            "assembly": "GRCh38",
        }

    monkeypatch.setattr(dl, "resolve_encode_file", resolve)
    paths = project_paths(tmp_path)
    session = _Session(
        {"http://s3/ENCFF000PKA": PEAKS, "http://direct/pkb": PEAKS}
    )
    got = dl.download_encode_chromatin_peaks(
        spec, paths, session=session
    )
    assert [p.name for p in got] == [
        "H3K4me3_ENCFF000PKA.narrowPeak.gz",
        "H3K27ac_ENCFF000PKB.broadPeak.gz",
    ]
    assert all(
        p.parent == paths.encode_chromatin_peaks_dir for p in got
    )
    assert (
        paths.encode_chromatin_peak_files(["narrowPeak", "broadPeak"])
        == got
    )

    session.fetched.clear()
    assert (
        dl.download_encode_chromatin_peaks(
            spec, paths, session=session
        )
        == got
    )
    assert session.fetched == []


def test_orchestration_reports_peaks(tmp_path, monkeypatch):
    monkeypatch.setattr(
        dl, "get_project", lambda code, path=None: _project()
    )
    monkeypatch.setattr(
        dl,
        "download_roadmap_tracks",
        lambda spec, paths, **kw: ["a"] * 4,
    )
    monkeypatch.setattr(
        dl,
        "download_roadmap_peaks",
        lambda spec, paths, **kw: ["p"] * 7,
    )
    monkeypatch.setattr(
        dl,
        "download_encode_chromatin_tracks",
        lambda spec, paths, **kw: ["a"],
    )
    monkeypatch.setattr(
        dl,
        "download_encode_chromatin_peaks",
        lambda spec, paths, **kw: ["p"],
    )
    report = dl.download_covariates(
        "X", tmp_path, which=("roadmap", "encode_chromatin")
    )
    assert report.sources["roadmap"] == {
        "status": "ok",
        "n_files": 4,
        "expected": 4,
        "n_peak_files": 7,
        "expected_peaks": 8,
    }
    assert report.sources["encode_chromatin"]["n_peak_files"] == 1
    assert report.sources["encode_chromatin"]["expected_peaks"] == 1


def _project():
    from types import SimpleNamespace

    return SimpleNamespace(
        code="X",
        roadmap=_roadmap(),
        encode_chromatin=EncodeChromatinSpec(
            tracks=(TrackRef(label="DNase", accession="ENCFF1"),),
            peak_tracks=(PeakRef(label="DNase", accession="ENCFF2"),),
        ),
    )
