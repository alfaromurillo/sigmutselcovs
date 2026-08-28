"""Zenodo publish tests (mocked HTTP -- never touches the real API)."""

import json

import pytest

from sigmutselcovs.publish import (
    create_deposition,
    default_metadata,
    get_deposition,
    publish_deposition,
    update_deposition_metadata,
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
    """Records every call made against it; returns scripted responses.

    Distinguishes a file upload PUT (bucket URL + `data=`, a raw file
    body) from a metadata-update PUT (deposition URL + `json=`, a
    {"metadata": ...} body) by which kwarg is passed, matching how
    the real requests.Session().put(...) calls in publish.py differ.
    """

    def __init__(self, deposition, created=None):
        self.deposition = deposition
        self.created = created or {
            "id": 555111,
            "links": {
                "bucket": "https://zenodo.example/api/files/newbucket"
            },
        }
        self.calls = []
        self.put_bodies = {}
        self.put_metadata = {}
        self.post_metadata = None

    def get(self, url, headers=None, timeout=None):
        self.calls.append(("GET", url, headers))
        return _Response(payload=self.deposition)

    def put(
        self, url, data=None, json=None, headers=None, timeout=None
    ):
        self.calls.append(("PUT", url, headers))
        if json is not None:
            self.put_metadata[url] = json
            updated = dict(self.deposition)
            updated["metadata"] = json.get("metadata")
            return _Response(payload=updated)
        filename = url.rsplit("/", 1)[-1]
        content = data.read() if hasattr(data, "read") else data
        self.put_bodies[filename] = content
        return _Response(
            payload={"size": len(content), "checksum": "md5:deadbeef"}
        )

    def post(self, url, json=None, headers=None, timeout=None):
        self.calls.append(("POST", url, headers))
        if url.endswith("/actions/publish"):
            return _Response(payload={"id": 998877, "state": "done"})
        self.post_metadata = json
        return _Response(payload=self.created)


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


def test_create_deposition_posts_metadata_and_returns_new_deposition():
    session = _Session(_deposition())
    metadata = {"title": "sigmutselcovs covariate matrices: TCGA-ACC"}

    result = create_deposition(
        "secret-token", metadata=metadata, session=session
    )

    assert result["id"] == 555111
    assert result["links"]["bucket"].endswith("newbucket")
    method, url, headers = session.calls[0]
    assert method == "POST"
    assert not url.endswith("/actions/publish")
    assert headers["Authorization"] == "Bearer secret-token"
    assert session.post_metadata == {"metadata": metadata}


def test_create_deposition_without_metadata_sends_empty_body():
    session = _Session(_deposition())
    create_deposition("secret-token", session=session)
    assert session.post_metadata == {}


def test_update_deposition_metadata_puts_to_deposition_url():
    session = _Session(_deposition())
    metadata = {"title": "updated title"}

    result = update_deposition_metadata(
        998877, metadata, "secret-token", session=session
    )

    assert result["metadata"] == metadata
    method, url, headers = session.calls[0]
    assert method == "PUT"
    assert url.endswith("/998877")
    assert headers["Authorization"] == "Bearer secret-token"
    assert session.put_metadata[url] == {"metadata": metadata}


def test_default_metadata_fills_in_project_and_component_count():
    metadata = default_metadata(
        "coad",
        "colon adenocarcinoma",
        n_components=481,
        cumulative_variance=0.9900301457204893,
    )

    assert (
        metadata["title"]
        == "sigmutselcovs covariate matrices: TCGA-COAD"
    )
    assert "TCGA-COAD" in metadata["keywords"]
    assert "colon adenocarcinoma" in metadata["description"]
    assert "481 components" in metadata["description"]
    assert "99%" in metadata["description"]
    assert metadata["license"] == "cc-by-4.0"
    assert metadata["upload_type"] == "dataset"
    # standard GitHub-repo link is always present
    assert any(
        r["relation"] == "isSupplementTo"
        for r in metadata["related_identifiers"]
    )


def test_default_metadata_includes_extra_related_identifiers():
    metadata = default_metadata(
        "ACC",
        "adrenocortical carcinoma",
        n_components=90,
        cumulative_variance=0.99,
        related_identifiers=[
            {
                "identifier": "10.5281/zenodo.21923082",
                "relation": "isPartOf",
                "resource_type": "dataset",
                "scheme": "doi",
            }
        ],
    )
    relations = {
        r["relation"] for r in metadata["related_identifiers"]
    }
    assert relations == {"isSupplementTo", "isPartOf"}
