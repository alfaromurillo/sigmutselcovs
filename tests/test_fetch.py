"""OSF fetch tests (mocked HTTP)."""

import hashlib
import json

import pandas as pd
import pytest

import sigmutselcovs.fetch as fetch
from sigmutselcovs.fetch import (
    CovariateArtifactsUnavailable,
    fetch_covariate_artifacts,
    fetch_covariate_matrix,
    osf_available,
)


def test_unconfigured_osf_raises_actionable_error():
    with pytest.raises(CovariateArtifactsUnavailable,
                       match="not published yet"):
        fetch.osf_index(config={"index_url": None})
    # subclasses FileNotFoundError so callers can catch broadly
    assert issubclass(CovariateArtifactsUnavailable, FileNotFoundError)
    assert not osf_available("COAD", config={"index_url": None})


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
    def __init__(self, index, files):
        self.index = index
        self.files = files

    def get(self, url, stream=False, headers=None, timeout=None):
        if url.endswith("index.json"):
            return _Response(payload=self.index)
        return _Response(content=self.files[url])


def _setup(tmp_path):
    csv = b"ensembl_gene_id,PC1\nENSG00000000001,0.5\n"
    index = {"schema_version": 1, "projects": {"COAD": {
        "layout_version": "v1",
        "files": [{"name": "cov_matrix_pca.parquet",
                   "url": "http://osf/pca",
                   "md5": hashlib.md5(b"parquet").hexdigest(),
                   "size": len(b"parquet")},
                  {"name": "cov_matrix_simple.csv",
                   "url": "http://osf/simple",
                   "md5": hashlib.md5(csv).hexdigest(),
                   "size": len(csv)}]}}}
    session = _Session(index, {"http://osf/pca": b"parquet",
                               "http://osf/simple": csv})
    config = {"index_url": "http://osf/index.json",
              "layout_version": "v1"}
    return session, config


def test_fetch_artifacts(tmp_path):
    session, config = _setup(tmp_path)
    got = fetch_covariate_artifacts(
        "coad", artifacts=("cov_matrix_pca",), dest=tmp_path,
        config=config, session=session)
    assert got["cov_matrix_pca"].read_bytes() == b"parquet"


def test_fetch_matrix_csv(tmp_path):
    session, config = _setup(tmp_path)
    frame = fetch_covariate_matrix("COAD", "simple", dest=tmp_path,
                                   config=config, session=session)
    assert isinstance(frame, pd.DataFrame)
    assert frame.loc["ENSG00000000001", "PC1"] == 0.5


def test_unpublished_project_raises(tmp_path):
    session, config = _setup(tmp_path)
    with pytest.raises(CovariateArtifactsUnavailable,
                       match="published: \\['COAD'\\]"):
        fetch_covariate_artifacts("BRCA", dest=tmp_path,
                                  config=config, session=session)


def test_missing_artifact_raises(tmp_path):
    session, config = _setup(tmp_path)
    with pytest.raises(CovariateArtifactsUnavailable,
                       match="cov_matrix_full.parquet"):
        fetch_covariate_artifacts(
            "COAD", artifacts=("cov_matrix_full",), dest=tmp_path,
            config=config, session=session)


def test_unknown_artifact_rejected(tmp_path):
    session, config = _setup(tmp_path)
    with pytest.raises(ValueError, match="Unknown artifacts"):
        fetch_covariate_artifacts("COAD", artifacts=("bogus",),
                                  dest=tmp_path, config=config,
                                  session=session)


def test_md5_mismatch_detected(tmp_path):
    session, config = _setup(tmp_path)
    session.files["http://osf/pca"] = b"corrupted!!"
    with pytest.raises(OSError):
        fetch_covariate_artifacts(
            "COAD", artifacts=("cov_matrix_pca",), dest=tmp_path,
            config=config, session=session)
