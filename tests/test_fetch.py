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


# --- column-selective reads of the full matrix -----------------------


class _RangeResponse:
    def __init__(self, content, status_code):
        self.content = content
        self.status_code = status_code

    def raise_for_status(self):
        pass


class _RangeSession(_Session):
    """Serves the record, the column dictionary, and HTTP ranges of
    the full parquet, counting every byte of the parquet it sends."""

    def __init__(self, record, files, full_url, full_bytes):
        super().__init__(record, files)
        self.full_url = full_url
        self.full_bytes = full_bytes
        self.full_requests = 0

    def get(self, url, stream=False, headers=None, timeout=None):
        if url == self.full_url:
            self.full_requests += 1
            assert headers and "Range" in headers, "full download"
            a, b = headers["Range"].removeprefix("bytes=").split("-")
            return _RangeResponse(
                self.full_bytes[int(a) : int(b) + 1], 206
            )
        return super().get(url, stream=stream, headers=headers)


def _setup_full(tmp_path, n_genes=500, n_cols=120):
    import numpy as np

    rng = np.random.default_rng(0)
    sources = ["roadmap", "atac", "gtex", "mrt"]
    names = [f"{sources[i % 4]}_{i}" for i in range(n_cols)]
    full = pd.DataFrame(
        rng.normal(size=(n_genes, n_cols)),
        index=pd.Index(
            [f"ENSG{i:011d}" for i in range(n_genes)],
            name="ensembl_gene_id",
        ),
        columns=names,
    )
    path = tmp_path / "full.parquet"
    full.to_parquet(path)
    full_bytes = path.read_bytes()
    dictionary = (
        pd.DataFrame(
            {
                "column": names,
                "source": [n.split("_")[0] for n in names],
            }
        )
        .to_csv(index=False)
        .encode()
    )
    record = {
        "id": 998877,
        "files": [
            {
                "key": "cov_matrix_full.parquet",
                "size": len(full_bytes),
                "checksum": "md5:"
                + hashlib.md5(full_bytes).hexdigest(),
                "links": {"self": "http://zenodo/full"},
            },
            {
                "key": "cov_matrix_columns.csv",
                "size": len(dictionary),
                "checksum": "md5:"
                + hashlib.md5(dictionary).hexdigest(),
                "links": {"self": "http://zenodo/columns"},
            },
        ],
    }
    session = _RangeSession(
        record,
        {"http://zenodo/columns": dictionary},
        "http://zenodo/full",
        full_bytes,
    )
    config = {
        "api_url": "http://zenodo/api/records",
        "records": {"COAD": 998877},
        "layout_version": "v1",
    }
    return full, session, config


def test_columns_match_a_full_download(tmp_path):
    full, session, config = _setup_full(tmp_path)
    wanted = ["atac_1", "roadmap_8", "mrt_3"]
    got = fetch_covariate_matrix(
        "COAD",
        "full",
        columns=wanted,
        dest=tmp_path / "cache",
        config=config,
        session=session,
    )
    # dictionary order, gene index kept
    expected = full[["atac_1", "mrt_3", "roadmap_8"]]
    pd.testing.assert_frame_equal(got, expected)


def test_columns_transfer_a_small_fraction(tmp_path):
    _, session, config = _setup_full(tmp_path)
    got = fetch_covariate_matrix(
        "COAD",
        "full",
        columns=["gtex_2"],
        dest=tmp_path / "cache",
        config=config,
        session=session,
    )
    # the point is the transfer, so assert on the bytes moved
    assert got.attrs["bytes_read"] < 0.15 * len(session.full_bytes)


def test_blocks_select_by_recorded_source(tmp_path):
    full, session, config = _setup_full(tmp_path)
    got = fetch_covariate_matrix(
        "COAD",
        "full",
        blocks=["chromatin", "atac"],
        dest=tmp_path / "cache",
        config=config,
        session=session,
    )
    expected = [
        c for c in full.columns if c.startswith(("roadmap", "atac"))
    ]
    assert list(got.columns) == expected


def test_unknown_names_raise_before_any_transfer(tmp_path):
    _, session, config = _setup_full(tmp_path)
    with pytest.raises(ValueError, match="atac_1"):
        fetch_covariate_matrix(
            "COAD",
            "full",
            columns=["atac1"],
            dest=tmp_path / "cache",
            config=config,
            session=session,
        )
    with pytest.raises(ValueError, match="Unknown covariate blocks"):
        fetch_covariate_matrix(
            "COAD",
            "full",
            blocks=["histone"],
            dest=tmp_path / "cache",
            config=config,
            session=session,
        )
    assert session.full_requests == 0


def test_columns_need_the_full_matrix(tmp_path):
    _, session, config = _setup_full(tmp_path)
    with pytest.raises(ValueError, match="which='full'"):
        fetch_covariate_matrix(
            "COAD",
            "pca",
            columns=["atac_1"],
            config=config,
            session=session,
        )
