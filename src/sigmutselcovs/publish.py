"""Publish covariate-matrix artifacts to a Zenodo deposition.

Uploads the same file set :mod:`fetch` expects to download later
(``cov_matrix_pca.parquet``, ``cov_matrix_full.parquet``,
``cov_matrix_columns.csv``, ``build_manifest.json``,
``pca_manifest.json``) to an existing Zenodo *deposition* via its
bucket API, replacing any existing file of the same name -- safe to
call repeatedly against a still-draft deposition while getting the
content right.

Deliberately does **not** read credentials from any particular local
secret-storage convention (a `.authinfo.gpg`, a keyring, an env var
by a fixed name, ...) -- callers pass a token explicitly, which keeps
this module usable in tests and CI without assuming any one setup.

Publishing a deposition (:func:`publish_deposition`) is a *separate*,
deliberate call: once published, a Zenodo record cannot be
unpublished, only superseded by a new version. Uploading files does
not publish anything by itself.
"""

import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

DEFAULT_API_URL = "https://zenodo.org/api/deposit/depositions"

# Matches fetch.py's _FILENAMES -- the artifact set this package
# advertises as downloadable. Not every project publishes every one
# (e.g. cov_matrix_simple.csv/cov_matrix_tcga.csv are commonly left
# out as internal ablation experiments); a filename missing locally
# is skipped, not an error.
ARTIFACT_FILENAMES = (
    "cov_matrix_pca.parquet",
    "cov_matrix_full.parquet",
    "cov_matrix_columns.csv",
    "build_manifest.json",
    "pca_manifest.json",
)


def _headers(token: str) -> dict:
    return {"Authorization": f"Bearer {token}"}


def get_deposition(
    deposit_id: int | str,
    token: str,
    *,
    api_url: str = DEFAULT_API_URL,
    session: requests.Session | None = None,
    timeout: int = 30,
) -> dict:
    """Fetch a deposition's current state (metadata, files, bucket URL)."""
    session = session or requests.Session()
    response = session.get(
        f"{api_url}/{deposit_id}",
        headers=_headers(token),
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()


def upload_artifact_files(
    deposit_id: int | str,
    matrices_dir: str | Path,
    token: str,
    *,
    filenames: tuple[str, ...] = ARTIFACT_FILENAMES,
    api_url: str = DEFAULT_API_URL,
    session: requests.Session | None = None,
    timeout: int = 600,
) -> dict[str, dict]:
    """Upload/replace files in an existing (draft) deposition's bucket.

    Parameters
    ----------
    deposit_id : int | str
        The Zenodo deposition id (same numeric id the published
        record will use).
    matrices_dir : str | Path
        A project's ``covariate_matrices/`` directory (as written by
        ``build_covariate_matrix(..., cache_matrices=True)`` and
        :func:`pca_artifact.save_pca_artifact`).
    token : str
        Zenodo personal access token with ``deposit:write`` scope.
    filenames : tuple[str, ...]
        Which artifact filenames to look for and upload.

    Returns
    -------
    dict[str, dict]
        Filename -> {"size": int, "checksum": str} for each file
        actually uploaded. Filenames not found in ``matrices_dir``
        are logged and skipped, not raised.
    """
    session = session or requests.Session()
    deposition = get_deposition(
        deposit_id,
        token,
        api_url=api_url,
        session=session,
        timeout=timeout,
    )
    bucket_url = deposition["links"]["bucket"]

    matrices_dir = Path(matrices_dir)
    results: dict[str, dict] = {}
    for filename in filenames:
        path = matrices_dir / filename
        if not path.exists():
            logger.info(
                "Skipping %s: not found in %s", filename, matrices_dir
            )
            continue
        with open(path, "rb") as f:
            response = session.put(
                f"{bucket_url}/{filename}",
                data=f,
                headers=_headers(token),
                timeout=timeout,
            )
        response.raise_for_status()
        meta = response.json()
        results[filename] = {
            "size": meta.get("size"),
            "checksum": meta.get("checksum"),
        }
        logger.info(
            "Uploaded %s (%s bytes)", filename, meta.get("size")
        )
    return results


def publish_deposition(
    deposit_id: int | str,
    token: str,
    *,
    api_url: str = DEFAULT_API_URL,
    session: requests.Session | None = None,
    timeout: int = 30,
) -> dict:
    """Publish a deposition -- **irreversible**.

    Only call this once the uploaded files have been verified (e.g.
    via :func:`get_deposition`'s file listing/checksums). A published
    Zenodo record cannot be unpublished, only superseded by a new
    version.
    """
    session = session or requests.Session()
    response = session.post(
        f"{api_url}/{deposit_id}/actions/publish",
        headers=_headers(token),
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()
