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
