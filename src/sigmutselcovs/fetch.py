"""Fetch pre-computed covariate matrices from OSF.

Most users should not need to download ~30 GB of raw tracks and
rebuild: the matrices for supported projects are published on OSF.
The default artifact is the PCA-reduced matrix (a few MB); the full
matrix (~0.5-1 GB per project), the simple/tcga matrices, the column
dictionary, and the per-gene intermediates are opt-in.

Layout on OSF (osfstorage), addressed through ``covariates/
index.json``:

    covariates/index.json
    covariates/<PROJECT>/<layout_version>/build_manifest.json
    covariates/<PROJECT>/<layout_version>/cov_matrix_pca.parquet
    covariates/<PROJECT>/<layout_version>/cov_matrix_full.parquet
    covariates/<PROJECT>/<layout_version>/cov_matrix_simple.csv
    covariates/<PROJECT>/<layout_version>/cov_matrix_tcga.csv
    covariates/<PROJECT>/<layout_version>/cov_matrix_columns.csv
    covariates/<PROJECT>/<layout_version>/intermediates/...

``index.json`` maps each project to its published files (name, size,
md5, download URL).  The OSF project id / index URL live in
``data/osf.json`` — currently placeholders until the artifacts are
first uploaded.
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

location_osf_config = location_covariates_data / "osf.json"

_ARTIFACTS = ("cov_matrix_pca", "cov_matrix_full",
              "cov_matrix_simple", "cov_matrix_tcga",
              "cov_matrix_columns", "build_manifest")
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
    base = os.environ.get("XDG_CACHE_HOME",
                          str(Path.home() / ".cache"))
    return Path(base) / "sigmutselcovs"


def load_osf_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path is not None else location_osf_config
    return json.loads(path.read_text())


def osf_index(*,
              config: dict | None = None,
              session: requests.Session | None = None,
              timeout: int = 30) -> dict:
    """Fetch and return the OSF artifact index."""
    config = config or load_osf_config()
    index_url = config.get("index_url")
    if not index_url:
        raise CovariateArtifactsUnavailable(
            "Pre-built covariate artifacts are not published yet "
            "(index_url is not configured in sigmutselcovs/data/"
            "osf.json).\nBuild locally instead:\n"
            "    sigmutselcovs download <PROJECT> --data-dir <dir>\n"
            "    sigmutselcovs build <PROJECT> --data-dir <dir>")
    session = session or requests.Session()
    response = session.get(index_url, timeout=timeout)
    response.raise_for_status()
    return response.json()


def osf_available(project: str, **kwargs) -> bool:
    """Whether pre-built artifacts exist for a project."""
    try:
        index = osf_index(**kwargs)
    except CovariateArtifactsUnavailable:
        return False
    return project.upper() in index.get("projects", {})


def fetch_covariate_artifacts(
        project: str,
        *,
        artifacts: tuple[str, ...] = ("cov_matrix_pca",),
        dest: str | Path | None = None,
        force: bool = False,
        config: dict | None = None,
        session: requests.Session | None = None
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
        raise ValueError(f"Unknown artifacts: {sorted(unknown)}; "
                         f"valid: {_ARTIFACTS}")
    project = project.upper()
    config = config or load_osf_config()
    index = osf_index(config=config, session=session)
    projects = index.get("projects", {})
    if project not in projects:
        raise CovariateArtifactsUnavailable(
            f"No published covariate artifacts for {project}; "
            f"published: {sorted(projects) or 'none'}.\n"
            "Build locally instead:\n"
            f"    sigmutselcovs download {project} --data-dir <dir>\n"
            f"    sigmutselcovs build {project} --data-dir <dir>")
    entry = projects[project]
    files = {f["name"]: f for f in entry.get("files", [])}

    layout = entry.get("layout_version",
                       config.get("layout_version", "v1"))
    dest = (Path(dest) if dest is not None
            else _cache_dir() / project / layout)
    dest.mkdir(parents=True, exist_ok=True)

    session = session or requests.Session()
    out: dict[str, Path] = {}
    for artifact in artifacts:
        filename = _FILENAMES[artifact]
        if filename not in files:
            raise CovariateArtifactsUnavailable(
                f"{project}: artifact {filename} not in the "
                f"published set {sorted(files)}")
        meta = files[filename]
        out[artifact] = download_file(
            meta["url"], dest / filename,
            expected_md5=meta.get("md5"),
            expected_size=meta.get("size"),
            force=force, session=session)
    return out


def fetch_covariate_matrix(project: str,
                           which: str = "pca",
                           *,
                           dest: str | Path | None = None,
                           force: bool = False,
                           **kwargs) -> pd.DataFrame:
    """Fetch one published covariate matrix as a DataFrame.

    Parameters
    ----------
    which : str
        'pca' (default, a few MB), 'full', 'simple', or 'tcga'.
    """
    artifact = f"cov_matrix_{which}"
    paths = fetch_covariate_artifacts(
        project, artifacts=(artifact,), dest=dest, force=force,
        **kwargs)
    path = paths[artifact]
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, index_col=0)
