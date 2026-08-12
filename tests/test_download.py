"""Download layer tests (mocked HTTP, no network)."""

import hashlib

import pytest

from sigmutselcovs.download import download_file


class _FakeResponse:
    def __init__(self, content: bytes, status_code: int = 200):
        self.content = content
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise OSError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size):
        for i in range(0, len(self.content), chunk_size):
            yield self.content[i:i + chunk_size]

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class _FakeSession:
    """Serves `payload`; honours Range when `supports_range`."""

    def __init__(self, payload: bytes, supports_range: bool = True,
                 fail_after: int | None = None):
        self.payload = payload
        self.supports_range = supports_range
        self.fail_after = fail_after
        self.requests: list[dict] = []

    def get(self, url, stream=True, headers=None, timeout=None):
        headers = headers or {}
        self.requests.append(headers)
        body = self.payload
        status = 200
        if "Range" in headers and self.supports_range:
            offset = int(headers["Range"].split("=")[1].rstrip("-"))
            body = body[offset:]
            status = 206
        if self.fail_after is not None:
            body = body[:self.fail_after]
            self.fail_after = None
        return _FakeResponse(body, status)


PAYLOAD = b"x" * 100 + b"y" * 50


def test_download_atomic_and_verified(tmp_path):
    dest = tmp_path / "file.bin"
    session = _FakeSession(PAYLOAD)
    got = download_file("http://example/f", dest,
                        expected_md5=hashlib.md5(PAYLOAD).hexdigest(),
                        expected_size=len(PAYLOAD),
                        session=session)
    assert got == dest
    assert dest.read_bytes() == PAYLOAD
    assert not dest.with_suffix(".bin.part").exists()


def test_skip_if_present(tmp_path):
    dest = tmp_path / "file.bin"
    dest.write_bytes(PAYLOAD)
    session = _FakeSession(PAYLOAD)
    download_file("http://example/f", dest, session=session,
                  expected_size=len(PAYLOAD))
    assert session.requests == []  # no HTTP at all


def test_refetch_on_size_mismatch(tmp_path):
    dest = tmp_path / "file.bin"
    dest.write_bytes(b"truncated")
    session = _FakeSession(PAYLOAD)
    download_file("http://example/f", dest, session=session,
                  expected_size=len(PAYLOAD))
    assert dest.read_bytes() == PAYLOAD


def test_resume_from_part(tmp_path):
    dest = tmp_path / "file.bin"
    part = tmp_path / "file.bin.part"
    part.write_bytes(PAYLOAD[:60])
    session = _FakeSession(PAYLOAD)
    download_file("http://example/f", dest, session=session)
    assert dest.read_bytes() == PAYLOAD
    assert session.requests[0].get("Range") == "bytes=60-"


def test_range_ignored_falls_back_to_full(tmp_path):
    dest = tmp_path / "file.bin"
    part = tmp_path / "file.bin.part"
    part.write_bytes(PAYLOAD[:60])
    session = _FakeSession(PAYLOAD, supports_range=False)
    download_file("http://example/f", dest, session=session)
    assert dest.read_bytes() == PAYLOAD


def test_md5_mismatch_raises_and_removes_part(tmp_path):
    dest = tmp_path / "file.bin"
    session = _FakeSession(PAYLOAD)
    with pytest.raises(OSError, match="md5 mismatch"):
        download_file("http://example/f", dest, session=session,
                      expected_md5="0" * 32)
    assert not dest.exists()
    assert not (tmp_path / "file.bin.part").exists()


class _URLSession:
    """Serves a url -> bytes mapping; unknown urls raise."""

    def __init__(self, mapping):
        self.mapping = mapping
        self.fetched: list[str] = []

    def get(self, url, stream=True, headers=None, timeout=None):
        self.fetched.append(url)
        if url not in self.mapping:
            return _FakeResponse(b"", status_code=404)
        return _FakeResponse(self.mapping[url])


def test_download_repliseq_mat_gunzips(tmp_path):
    import gzip as _gzip

    import sigmutselcovs.download as dl
    from sigmutselcovs.paths import project_paths
    from sigmutselcovs.registry import RepliseqSpec

    table = b"chr1\t0\t50000\n"
    spec = RepliseqSpec(type="mat", assembly="hg38", cell_line="HCT116",
                        filename="rt.mat", url="http://geo/rt.mat.gz")
    paths = project_paths(tmp_path)
    session = _URLSession({"http://geo/rt.mat.gz":
                           _gzip.compress(table)})
    got = dl.download_repliseq(spec, paths, session=session)
    assert got == [paths.rt_dir / "rt.mat"]
    assert got[0].read_bytes() == table
    assert not (paths.rt_dir / "rt.mat.gz").exists()
    # idempotent: second call fetches nothing
    dl.download_repliseq(spec, paths, session=session)
    assert len(session.fetched) == 1


def test_download_repliseq_bigwigs(tmp_path, monkeypatch):
    import sigmutselcovs.download as dl
    from sigmutselcovs.paths import project_paths
    from sigmutselcovs.registry import RepliseqSpec, TrackRef

    payload = b"bigwigbytes"
    monkeypatch.setattr(
        dl, "resolve_encode_file",
        lambda acc, session=None: {
            "accession": acc,
            "url": f"http://s3/{acc}",
            "md5sum": hashlib.md5(payload).hexdigest(),
            "file_size": len(payload),
            "assembly": "hg19",
        })
    spec = RepliseqSpec(
        type="fraction_bigwigs", assembly="hg19", cell_line="MCF-7",
        tracks=(TrackRef(label="s1", accession="ENCFF000AAA"),
                TrackRef(label="s2", accession="ENCFF000BBB")))
    paths = project_paths(tmp_path)
    session = _URLSession({"http://s3/ENCFF000AAA": payload,
                           "http://s3/ENCFF000BBB": payload})
    got = dl.download_repliseq(spec, paths, session=session)
    assert [p.name for p in got] == ["ENCFF000AAA.bigWig",
                                     "ENCFF000BBB.bigWig"]
    assert all(p.read_bytes() == payload for p in got)


def test_download_roadmap_tolerates_missing_marks(tmp_path):
    import sigmutselcovs.download as dl
    from sigmutselcovs.paths import project_paths
    from sigmutselcovs.registry import RoadmapSpec

    template = "http://roadmap/{eid}-{mark}.fc.signal.bigwig"
    spec = RoadmapSpec(eids=("E027",),
                       marks=("H3K4me1", "H3K27ac"),
                       url_template=template)
    paths = project_paths(tmp_path)
    session = _URLSession(
        {template.format(eid="E027", mark="H3K4me1"): b"data"})
    got = dl.download_roadmap_tracks(spec, paths, session=session)
    assert [p.name for p in got] == ["E027-H3K4me1.fc.signal.bigwig"]

    required = RoadmapSpec(eids=("E027",),
                           marks=("H3K27ac",),
                           required_marks=("H3K27ac",),
                           url_template=template)
    with pytest.raises(Exception):
        dl.download_roadmap_tracks(required, paths, session=session)


def _make_atac_tarball(payloads: dict[str, bytes]) -> bytes:
    import io
    import tarfile

    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        for name, data in payloads.items():
            import time
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = int(time.time())
            tar.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


def test_download_tcga_atac_flattens(tmp_path):
    import sigmutselcovs.download as dl
    from sigmutselcovs.paths import project_paths
    from sigmutselcovs.registry import AtacSpec

    tar_bytes = _make_atac_tarball({
        "oak/stanford/deep/path/COAD_A_T1.insertions.bw": b"aaa",
        "oak/stanford/deep/path/COAD_B_T1.insertions.bw": b"bbb",
        "oak/stanford/deep/path/README.txt": b"skip me",
        "../evil.bw": b"nope",
    })
    spec = AtacSpec(gdc_uuid="uuid-1234", column_prefix="coad")
    paths = project_paths(tmp_path)
    session = _URLSession({f"{dl.GDC_DATA_URL}/uuid-1234": tar_bytes})
    got = dl.download_tcga_atac(spec, paths, session=session)
    assert [p.name for p in got] == ["COAD_A_T1.insertions.bw",
                                     "COAD_B_T1.insertions.bw"]
    assert (paths.atac_dir / "COAD_A_T1.insertions.bw"
            ).read_bytes() == b"aaa"
    assert not (paths.atac_dir / "README.txt").exists()
    assert not (tmp_path / "evil.bw").exists()
    assert not list(paths.atac_dir.glob("*.tgz"))
    # idempotent once bigWigs exist
    dl.download_tcga_atac(spec, paths, session=session)
    assert len(session.fetched) == 1


def test_download_tcga_atac_empty_tarball_raises(tmp_path):
    import sigmutselcovs.download as dl
    from sigmutselcovs.paths import project_paths
    from sigmutselcovs.registry import AtacSpec

    tar_bytes = _make_atac_tarball({"just/a/README.txt": b"x"})
    spec = AtacSpec(gdc_uuid="uuid-9", column_prefix="coad")
    paths = project_paths(tmp_path)
    session = _URLSession({f"{dl.GDC_DATA_URL}/uuid-9": tar_bytes})
    with pytest.raises(OSError, match="No bigWig members"):
        dl.download_tcga_atac(spec, paths, session=session)


def test_short_download_keeps_part_for_resume(tmp_path):
    dest = tmp_path / "file.bin"
    session = _FakeSession(PAYLOAD, fail_after=80)
    with pytest.raises(OSError, match="size"):
        download_file("http://example/f", dest, session=session,
                      expected_size=len(PAYLOAD), resume=False)
    assert not dest.exists()
    part = tmp_path / "file.bin.part"
    assert part.exists() and part.stat().st_size == 80
    # second attempt resumes and completes
    download_file("http://example/f", dest, session=session,
                  expected_size=len(PAYLOAD))
    assert dest.read_bytes() == PAYLOAD
