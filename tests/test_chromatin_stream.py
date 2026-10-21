"""A streamed cache equals a build from bigWigs on disk."""

import json
import shutil

import pandas as pd
from conftest import write_bigwig, write_gtf

from sigmutselcovs import chromatin_stream as cs
from sigmutselcovs.builder import build_covariate_matrix
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import location_projects_registry

TRACKS = ("E075-H3K4me3", "E106-H3K4me3")


def _registry(tmp_path):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["roadmap"].update(
        eids=["E075", "E106"],
        marks=["H3K4me3"],
        windows=["cds_span", "body"],
        transform=["log1p_mean", None],
    )
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    return path


def test_streamed_cache_equals_the_batch_build(tmp_path, monkeypatch):
    gtf = write_gtf(tmp_path / "g.gtf")
    gtfs = {"hg19": gtf, "hg38": gtf}
    reg = _registry(tmp_path)
    store = tmp_path / "remote"
    store.mkdir()
    for seed, track in enumerate(TRACKS):
        write_bigwig(store / f"{track}.fc.signal.bigwig", seed)

    # the batch build, from bigWigs on disk
    batch = project_paths(tmp_path / "batch")
    batch.roadmap_dir.mkdir(parents=True)
    for f in store.iterdir():
        shutil.copy(f, batch.roadmap_dir / f.name)
    kwargs = {
        "registry_path": reg,
        "gencode_gtfs": gtfs,
        "include": ["roadmap"],
        "apply_fixes": False,
        "validate_output": False,
    }
    expected = build_covariate_matrix(
        "COAD", batch.root, **kwargs
    ).full

    # the streamed one: "downloads" copy from the store, and one track
    # of a missing (eid, mark) pair is skipped
    def fake_download(url, dest, **_):
        src = store / url.rsplit("/", 1)[1]
        if not src.exists():
            raise FileNotFoundError(url)
        shutil.copy(src, dest)
        return dest

    monkeypatch.setattr(cs, "download_file", fake_download)
    streamed = project_paths(tmp_path / "streamed")
    cache = cs.stream_chromatin_cache(
        "COAD",
        streamed.root,
        "roadmap",
        registry_path=reg,
        gencode_gtfs=gtfs,
        workers=2,
    )
    assert not list(streamed.roadmap_dir.glob("*.bigwig"))
    got = build_covariate_matrix("COAD", streamed.root, **kwargs).full
    pd.testing.assert_frame_equal(got, expected)
    assert cache.name == (
        "chromatin_covs_roadmap_cds_span-body_log1p_mean-linear.csv"
    )

    # resumable: a second run summarises nothing again
    calls = []
    monkeypatch.setattr(
        cs, "download_file", lambda *a, **k: calls.append(a) or None
    )
    cs.stream_chromatin_cache(
        "COAD",
        streamed.root,
        "roadmap",
        registry_path=reg,
        gencode_gtfs=gtfs,
    )
    assert calls == []


def test_missing_tracks_are_skipped(tmp_path, monkeypatch):
    gtf = write_gtf(tmp_path / "g.gtf")
    store = tmp_path / "remote"
    store.mkdir()
    write_bigwig(store / "E075-H3K4me3.fc.signal.bigwig")

    def fake_download(url, dest, **_):
        src = store / url.rsplit("/", 1)[1]
        if not src.exists():
            raise FileNotFoundError(url)
        shutil.copy(src, dest)
        return dest

    monkeypatch.setattr(cs, "download_file", fake_download)
    cache = cs.stream_chromatin_cache(
        "COAD",
        project_paths(tmp_path / "d").root,
        "roadmap",
        registry_path=_registry(tmp_path),
        gencode_gtfs={"hg19": gtf, "hg38": gtf},
    )
    columns = pd.read_csv(cache, index_col=0).columns
    assert all(c.startswith("e075_") for c in columns)


def test_cli_stream_chromatin(tmp_path, monkeypatch, capsys):
    from sigmutselcovs import cli

    seen = {}

    def fake(project, data_dir, source, **kwargs):
        seen.update(project=project, source=source, **kwargs)
        return tmp_path / "cache.csv"

    monkeypatch.setattr(cs, "stream_chromatin_cache", fake)
    code = cli.main(
        [
            "stream-chromatin",
            "GENERIC",
            "--data-dir",
            str(tmp_path),
            "--source",
            "roadmap",
            "--workers",
            "3",
            "--chromatin-windows",
            "cds_span",
            "promoter",
        ]
    )
    assert code == 0
    assert seen == {
        "project": "GENERIC",
        "source": "roadmap",
        "workers": 3,
        "chromatin_windows": ["cds_span", "promoter"],
    }
