"""Fetch pre-computed covariate matrices from Zenodo.

Most users should not need to download ~30 GB of raw tracks and
rebuild: the matrices for supported projects are published on
Zenodo, one record per project. The default artifact is the
PCA-reduced matrix (a few MB); the full matrix (~0.5-1 GB per
project), the simple/tcga matrices, and the column dictionary are
opt-in.

Files published in each project's Zenodo record:

    build_manifest.json
    cov_matrix_pca.parquet
    cov_matrix_full.parquet
    cov_matrix_simple.csv
    cov_matrix_tcga.csv
    cov_matrix_columns.csv

Record ids per project live in ``data/zenodo.json`` (``records``,
mapping PROJECT -> Zenodo record id; a new version of a record gets
a new id, so this file points at the latest). Unlike a hosted index
file, this mapping is entirely local: whether a project is published
is known without a network call, and only fetching its file
list/checksums requires one (``GET {api_url}/{record_id}``, Zenodo's own REST API).
"""

import difflib
import io
import json
import logging
import os
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
import requests

from .covariates_locations import location_covariates_data
from .download import download_file

logger = logging.getLogger(__name__)

location_zenodo_config = location_covariates_data / "zenodo.json"

_ARTIFACTS = (
    "cov_matrix_pca",
    "cov_matrix_full",
    "cov_matrix_simple",
    "cov_matrix_tcga",
    "cov_matrix_columns",
    "build_manifest",
)
_FILENAMES = {
    "cov_matrix_pca": "cov_matrix_pca.parquet",
    "cov_matrix_full": "cov_matrix_full.parquet",
    "cov_matrix_simple": "cov_matrix_simple.csv",
    "cov_matrix_tcga": "cov_matrix_tcga.csv",
    "cov_matrix_columns": "cov_matrix_columns.csv",
    "build_manifest": "build_manifest.json",
}


class CovariateArtifactsUnavailable(FileNotFoundError):
    """Raised when pre-built covariates are not (yet) published."""


def _cache_dir() -> Path:
    base = os.environ.get(
        "XDG_CACHE_HOME", str(Path.home() / ".cache")
    )
    return Path(base) / "sigmutselcovs"


def load_zenodo_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path is not None else location_zenodo_config
    return json.loads(path.read_text())


def zenodo_record(
    project: str,
    *,
    config: dict | None = None,
    session: requests.Session | None = None,
    timeout: int = 30,
) -> dict:
    """Fetch a project's Zenodo record (file list, sizes, checksums)."""
    config = config or load_zenodo_config()
    project = project.upper()
    records = config.get("records", {})
    record_id = records.get(project)
    if not record_id:
        raise CovariateArtifactsUnavailable(
            f"No published covariate artifacts for {project}; "
            f"published: {sorted(records) or 'none'}.\n"
            "Build locally instead:\n"
            f"    sigmutselcovs download {project} --data-dir <dir>\n"
            f"    sigmutselcovs build {project} --data-dir <dir>"
        )
    api_url = config.get("api_url", "https://zenodo.org/api/records")
    session = session or requests.Session()
    response = session.get(f"{api_url}/{record_id}", timeout=timeout)
    if response.status_code == 404:
        # A record id in `records` isn't necessarily published yet --
        # publishing a Zenodo deposition is a separate, deliberate
        # step from creating/uploading to it (see publish.py), so a
        # configured id can genuinely 404 on the public records API
        # for a while. Same clean error as "not configured at all"
        # rather than a raw HTTPError, since the caller can't do
        # anything different either way.
        raise CovariateArtifactsUnavailable(
            f"Covariate artifacts for {project} are configured "
            f"(record id {record_id}) but not published yet (404 "
            f"from the Zenodo API).\n"
            "Build locally instead:\n"
            f"    sigmutselcovs download {project} --data-dir <dir>\n"
            f"    sigmutselcovs build {project} --data-dir <dir>"
        )
    response.raise_for_status()
    return response.json()


def zenodo_available(project: str, **kwargs) -> bool:
    """Whether a published, reachable record exists for a project."""
    try:
        zenodo_record(project, **kwargs)
    except CovariateArtifactsUnavailable:
        return False
    return True


def fetch_covariate_artifacts(
    project: str,
    *,
    artifacts: tuple[str, ...] = ("cov_matrix_pca",),
    dest: str | Path | None = None,
    force: bool = False,
    config: dict | None = None,
    session: requests.Session | None = None,
) -> dict[str, Path]:
    """Download published artifacts for a project.

    Parameters
    ----------
    artifacts : tuple[str, ...]
        Any of cov_matrix_pca, cov_matrix_full, cov_matrix_simple,
        cov_matrix_tcga, cov_matrix_columns, build_manifest.
    dest : str | Path | None
        Where to put the files; defaults to
        ``$XDG_CACHE_HOME/sigmutselcovs/<PROJECT>/<layout_version>/``.
        Pass a project data dir's ``covariate_matrices/`` to make a
        fetch and a local build interchangeable.

    Returns
    -------
    dict[str, Path]
        Artifact name -> local path.
    """
    unknown = set(artifacts) - set(_ARTIFACTS)
    if unknown:
        raise ValueError(
            f"Unknown artifacts: {sorted(unknown)}; "
            f"valid: {_ARTIFACTS}"
        )
    project = project.upper()
    config = config or load_zenodo_config()
    session = session or requests.Session()
    record = zenodo_record(project, config=config, session=session)
    files = {f["key"]: f for f in record.get("files", [])}

    layout = config.get("layout_version", "v1")
    dest = (
        Path(dest)
        if dest is not None
        else _cache_dir() / project / layout
    )
    dest.mkdir(parents=True, exist_ok=True)

    out: dict[str, Path] = {}
    for artifact in artifacts:
        filename = _FILENAMES[artifact]
        if filename not in files:
            raise CovariateArtifactsUnavailable(
                f"{project}: artifact {filename} not in the "
                f"published set {sorted(files)}"
            )
        meta = files[filename]
        checksum = meta.get("checksum")  # Zenodo format: "md5:<hex>"
        md5 = checksum.split(":", 1)[1] if checksum else None
        out[artifact] = download_file(
            meta["links"]["self"],
            dest / filename,
            expected_md5=md5,
            expected_size=meta.get("size"),
            force=force,
            session=session,
        )
    return out


class HTTPRangeFile(io.RawIOBase):
    """A read-only, seekable file over HTTP range requests.

    Enough of a file for ``pyarrow.parquet`` to read a few columns of
    a remote parquet without downloading the rest: each ``readinto``
    is one ``Range: bytes=a-b`` request. ``bytes_read`` counts the
    bytes actually transferred.

    Zenodo serves range requests (``206 Partial Content``) although
    it does not advertise ``Accept-Ranges`` on a HEAD request, so the
    file size is taken from the record's metadata rather than from
    HEAD.
    """

    def __init__(
        self,
        url: str,
        size: int,
        *,
        session: requests.Session | None = None,
        timeout: int = 60,
    ):
        super().__init__()
        self.url = url
        self.size = int(size)
        self.session = session or requests.Session()
        self.timeout = timeout
        self.bytes_read = 0
        self._pos = 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self._pos

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            pos = offset
        elif whence == io.SEEK_CUR:
            pos = self._pos + offset
        elif whence == io.SEEK_END:
            pos = self.size + offset
        else:
            raise ValueError(f"invalid whence: {whence}")
        if pos < 0:
            raise ValueError("negative seek position")
        self._pos = pos
        return pos

    def readinto(self, buffer) -> int:
        view = memoryview(buffer).cast("B")
        if self._pos >= self.size or len(view) == 0:
            return 0
        end = min(self._pos + len(view), self.size) - 1
        response = self.session.get(
            self.url,
            headers={"Range": f"bytes={self._pos}-{end}"},
            timeout=self.timeout,
        )
        response.raise_for_status()
        if response.status_code != 206:
            raise OSError(
                f"{self.url} ignored the range request (HTTP "
                f"{response.status_code}); refusing to download the "
                "whole file through a column read."
            )
        data = response.content
        n = len(data)
        view[:n] = data
        self._pos += n
        self.bytes_read += n
        return n


def _resolve_columns(
    dictionary: pd.DataFrame,
    columns: Sequence[str] | None,
    blocks: Sequence[str] | None,
) -> list[str]:
    """Named columns plus every column of the named blocks, in the
    dictionary's order, validated before anything large is fetched."""
    from .builder import SOURCE_BLOCKS

    known = list(dictionary["column"])
    wanted: set[str] = set()
    if columns is not None:
        missing = [c for c in columns if c not in set(known)]
        if missing:
            hints = {
                c: difflib.get_close_matches(c, known, n=3)
                for c in missing
            }
            raise ValueError(
                "Unknown covariate columns (closest matches): "
                + "; ".join(
                    f"{c} -> {hints[c] or 'none'}" for c in missing
                )
            )
        wanted.update(columns)
    if blocks is not None:
        block_of = dictionary["source"].map(SOURCE_BLOCKS)
        block_of = block_of.fillna("other")
        valid = sorted(set(block_of))
        unknown = sorted(set(blocks) - set(valid))
        if unknown:
            raise ValueError(
                f"Unknown covariate blocks {unknown}; this record "
                f"has: {valid}"
            )
        wanted.update(dictionary.loc[block_of.isin(blocks), "column"])
    return [c for c in known if c in wanted]


def _fetch_full_columns(
    project: str,
    columns: Sequence[str] | None,
    blocks: Sequence[str] | None,
    *,
    dest: str | Path | None,
    force: bool,
    config: dict | None = None,
    session: requests.Session | None = None,
) -> pd.DataFrame:
    import pyarrow.parquet as pq

    config = config or load_zenodo_config()
    session = session or requests.Session()
    dictionary_path = fetch_covariate_artifacts(
        project,
        artifacts=("cov_matrix_columns",),
        dest=dest,
        force=force,
        config=config,
        session=session,
    )["cov_matrix_columns"]
    dictionary = pd.read_csv(dictionary_path)
    selected = _resolve_columns(dictionary, columns, blocks)
    if not selected:
        raise ValueError("The selection matches no columns.")

    record = zenodo_record(project, config=config, session=session)
    files = {f["key"]: f for f in record.get("files", [])}
    meta = files.get(_FILENAMES["cov_matrix_full"])
    if meta is None:
        raise CovariateArtifactsUnavailable(
            f"{project.upper()}: cov_matrix_full.parquet not in the "
            f"published set {sorted(files)}"
        )
    remote = HTTPRangeFile(
        meta["links"]["self"], meta["size"], session=session
    )
    # Unbuffered on purpose: a read-ahead buffer would fetch bytes
    # beyond each column chunk, which is what this path avoids.
    table = pq.read_table(
        remote,
        columns=selected,
        use_pandas_metadata=True,
    )
    logger.info(
        "%s: read %d columns of cov_matrix_full.parquet, %.1f MB of "
        "%.1f MB",
        project.upper(),
        len(selected),
        remote.bytes_read / 1e6,
        remote.size / 1e6,
    )
    frame = table.to_pandas()
    frame.attrs["bytes_read"] = remote.bytes_read
    return frame


def fetch_covariate_matrix(
    project: str,
    which: str = "pca",
    *,
    columns: Sequence[str] | None = None,
    blocks: Sequence[str] | None = None,
    dest: str | Path | None = None,
    force: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Fetch one published covariate matrix as a DataFrame.

    Parameters
    ----------
    which : str
        'pca' (default, a few MB), 'full', 'simple', or 'tcga'.
    columns : sequence of str, optional
        Named columns of the full matrix to read, with ``which='full'``.
        Only those columns' bytes (plus the parquet footer) are
        transferred, by HTTP range requests, instead of the whole
        ~0.5-1 GB file; nothing is cached. Names are checked against
        the record's ``cov_matrix_columns.csv`` before any large
        transfer.
    blocks : sequence of str, optional
        Whole covariate blocks of the full matrix, e.g.
        ``["chromatin", "atac"]``, resolved through each column's
        recorded source (`builder.SOURCE_BLOCKS`; ``"other"`` for a
        column with none). Combines with ``columns``.
    """
    if columns is not None or blocks is not None:
        if which != "full":
            raise ValueError(
                "columns= and blocks= select from the full matrix; "
                "pass which='full'."
            )
        return _fetch_full_columns(
            project, columns, blocks, dest=dest, force=force, **kwargs
        )
    artifact = f"cov_matrix_{which}"
    paths = fetch_covariate_artifacts(
        project,
        artifacts=(artifact,),
        dest=dest,
        force=force,
        **kwargs,
    )
    path = paths[artifact]
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, index_col=0)
