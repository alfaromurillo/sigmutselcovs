"""Zenodo publish tests (mocked HTTP -- never touches the real API)."""

import json

import pytest

from sigmutselcovs.publish import (
    get_deposition,
    publish_deposition,
    upload_artifact_files,
)


class _Response:
    def __init__(self, payload=None, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._payload


class _Session:
    """Records every call made against it; returns scripted responses."""

    def __init__(self, deposition):
        self.deposition = deposition
        self.calls = []
        self.put_bodies = {}

    def get(self, url, headers=None, timeout=None):
        self.calls.append(("GET", url, headers))
        return _Response(payload=self.deposition)

    def put(self, url, data=None, headers=None, timeout=None):
        self.calls.append(("PUT", url, headers))
        filename = url.rsplit("/", 1)[-1]
        content = data.read() if hasattr(data, "read") else data
        self.put_bodies[filename] = content
        return _Response(
            payload={"size": len(content), "checksum": "md5:deadbeef"}
        )

    def post(self, url, headers=None, timeout=None):
        self.calls.append(("POST", url, headers))
        return _Response(payload={"id": 998877, "state": "done"})


def _deposition():
    return {
        "id": 998877,
        "links": {
            "bucket": "https://zenodo.example/api/files/bucket123"
        },
    }


def test_get_deposition_uses_bearer_token():
    session = _Session(_deposition())
    result = get_deposition(998877, "secret-token", session=session)
    assert result["id"] == 998877
    method, url, headers = session.calls[0]
    assert method == "GET"
    assert url.endswith("/998877")
    assert headers["Authorization"] == "Bearer secret-token"


def test_upload_artifact_files_uploads_present_files_only(tmp_path):
    (tmp_path / "cov_matrix_pca.parquet").write_bytes(b"pca-bytes")
    (tmp_path / "build_manifest.json").write_text(
        json.dumps({"project": "TEST"})
    )
    # cov_matrix_full.parquet, cov_matrix_columns.csv, and
    # pca_manifest.json deliberately absent -- must be skipped, not
    # raise.

    session = _Session(_deposition())
    results = upload_artifact_files(
        998877, tmp_path, "secret-token", session=session
    )

    assert set(results) == {
        "cov_matrix_pca.parquet",
        "build_manifest.json",
    }
    assert (
        session.put_bodies["cov_matrix_pca.parquet"] == b"pca-bytes"
    )
    put_calls = [c for c in session.calls if c[0] == "PUT"]
    assert len(put_calls) == 2
    for _, url, headers in put_calls:
        assert url.startswith(
            "https://zenodo.example/api/files/bucket123/"
        )
        assert headers["Authorization"] == "Bearer secret-token"


def test_upload_artifact_files_no_files_present_uploads_nothing(
    tmp_path,
):
    session = _Session(_deposition())
    results = upload_artifact_files(
        998877, tmp_path, "secret-token", session=session
    )
    assert results == {}
    assert not any(c[0] == "PUT" for c in session.calls)


def test_publish_deposition_posts_to_actions_publish():
    session = _Session(_deposition())
    result = publish_deposition(
        998877, "secret-token", session=session
    )
    assert result["state"] == "done"
    method, url, headers = session.calls[0]
    assert method == "POST"
    assert url.endswith("/998877/actions/publish")
    assert headers["Authorization"] == "Bearer secret-token"


def test_upload_artifact_files_raises_on_http_error(tmp_path):
    (tmp_path / "cov_matrix_pca.parquet").write_bytes(b"pca-bytes")

    class _FailingSession(_Session):
        def put(self, url, data=None, headers=None, timeout=None):
            self.calls.append(("PUT", url, headers))
            return _Response(status_code=403)

    session = _FailingSession(_deposition())
    with pytest.raises(RuntimeError, match="403"):
        upload_artifact_files(
            998877, tmp_path, "secret-token", session=session
        )
