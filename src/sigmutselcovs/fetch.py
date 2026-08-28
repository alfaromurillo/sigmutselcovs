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
mapping PROJECT -> Zenodo record id) -- currently empty, until the
artifacts are first uploaded. Unlike a hosted index file, this
mapping is entirely local: whether a project is published is known
without a network call, and only fetching its file list/checksums
requires one (``GET {api_url}/{record_id}``, Zenodo's own REST API).
"""

import json
import logging
import os
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


def fetch_covariate_matrix(
    project: str,
    which: str = "pca",
    *,
    dest: str | Path | None = None,
    force: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Fetch one published covariate matrix as a DataFrame.

    Parameters
    ----------
    which : str
        'pca' (default, a few MB), 'full', 'simple', or 'tcga'.
    """
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
