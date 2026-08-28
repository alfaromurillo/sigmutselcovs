"""Zenodo fetch tests (mocked HTTP)."""

import hashlib

import pandas as pd
import pytest

from sigmutselcovs import fetch
from sigmutselcovs.fetch import (
    CovariateArtifactsUnavailable,
    fetch_covariate_artifacts,
    fetch_covariate_matrix,
    zenodo_available,
)


def test_unconfigured_project_raises_actionable_error():
    with pytest.raises(
        CovariateArtifactsUnavailable, match="published: none"
    ):
        fetch.zenodo_record("COAD", config={"records": {}})
    # subclasses FileNotFoundError so callers can catch broadly
    assert issubclass(
        CovariateArtifactsUnavailable, FileNotFoundError
    )
    assert not zenodo_available("COAD", config={"records": {}})


class _NotFoundResponse:
    status_code = 404

    def raise_for_status(self):
        raise AssertionError(
            "zenodo_record must check status_code == 404 before "
            "calling raise_for_status, not rely on it"
        )


class _NotFoundSession:
    def get(self, url, headers=None, timeout=None):
        return _NotFoundResponse()


def test_configured_but_unpublished_record_raises_actionable_error():
    """A record id can be wired into zenodo.json before the
    deposition is actually published (creating/uploading to a
    deposition and publishing it are separate steps) -- that should
    give the same clean, actionable error as an unconfigured
    project, not a raw HTTPError leaking the Zenodo API's 404."""
    config = {
        "api_url": "http://zenodo/api/records",
        "records": {"COAD": 21923082},
    }
    with pytest.raises(
        CovariateArtifactsUnavailable, match="not published yet"
    ):
        fetch.zenodo_record(
            "COAD", config=config, session=_NotFoundSession()
        )
    assert not zenodo_available(
        "COAD", config=config, session=_NotFoundSession()
    )


class _Response:
    def __init__(self, payload=None, content=b""):
        self._payload = payload
        self.content = content
        self.status_code = 200

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload

    def iter_content(self, chunk_size):
        yield self.content

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class _Session:
    def __init__(self, record, files):
        self.record = record
        self.files = files

    def get(self, url, stream=False, headers=None, timeout=None):
        if url.endswith("/998877"):
            return _Response(payload=self.record)
        return _Response(content=self.files[url])


def _setup(tmp_path):
    csv = b"ensembl_gene_id,PC1\nENSG00000000001,0.5\n"
    pca = b"parquet"
    record = {
        "id": 998877,
        "files": [
            {
                "key": "cov_matrix_pca.parquet",
                "size": len(pca),
                "checksum": "md5:" + hashlib.md5(pca).hexdigest(),
                "links": {"self": "http://zenodo/pca"},
            },
            {
                "key": "cov_matrix_simple.csv",
                "size": len(csv),
                "checksum": "md5:" + hashlib.md5(csv).hexdigest(),
                "links": {"self": "http://zenodo/simple"},
            },
        ],
    }
    session = _Session(
        record,
        {"http://zenodo/pca": pca, "http://zenodo/simple": csv},
    )
    config = {
        "api_url": "http://zenodo/api/records",
        "records": {"COAD": 998877},
        "layout_version": "v1",
    }
    return session, config


def test_fetch_artifacts(tmp_path):
    session, config = _setup(tmp_path)
    got = fetch_covariate_artifacts(
        "coad",
        artifacts=("cov_matrix_pca",),
        dest=tmp_path,
        config=config,
        session=session,
    )
    assert got["cov_matrix_pca"].read_bytes() == b"parquet"


def test_fetch_matrix_csv(tmp_path):
    session, config = _setup(tmp_path)
    frame = fetch_covariate_matrix(
        "COAD",
        "simple",
        dest=tmp_path,
        config=config,
        session=session,
    )
    assert isinstance(frame, pd.DataFrame)
    assert frame.loc["ENSG00000000001", "PC1"] == 0.5


def test_unpublished_project_raises(tmp_path):
    session, config = _setup(tmp_path)
    with pytest.raises(
        CovariateArtifactsUnavailable, match="published: \\['COAD'\\]"
    ):
        fetch_covariate_artifacts(
            "BRCA", dest=tmp_path, config=config, session=session
        )


def test_missing_artifact_raises(tmp_path):
    session, config = _setup(tmp_path)
    with pytest.raises(
        CovariateArtifactsUnavailable, match="cov_matrix_full.parquet"
    ):
        fetch_covariate_artifacts(
            "COAD",
            artifacts=("cov_matrix_full",),
            dest=tmp_path,
            config=config,
            session=session,
        )


def test_unknown_artifact_rejected(tmp_path):
    session, config = _setup(tmp_path)
    with pytest.raises(ValueError, match="Unknown artifacts"):
        fetch_covariate_artifacts(
            "COAD",
            artifacts=("bogus",),
            dest=tmp_path,
            config=config,
            session=session,
        )


def test_md5_mismatch_detected(tmp_path):
    session, config = _setup(tmp_path)
    session.files["http://zenodo/pca"] = b"corrupted!!"
    with pytest.raises(OSError):
        fetch_covariate_artifacts(
            "COAD",
            artifacts=("cov_matrix_pca",),
            dest=tmp_path,
            config=config,
            session=session,
        )
