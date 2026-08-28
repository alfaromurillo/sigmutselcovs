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


def create_deposition(
    token: str,
    *,
    metadata: dict | None = None,
    api_url: str = DEFAULT_API_URL,
    session: requests.Session | None = None,
    timeout: int = 30,
) -> dict:
    """Create a new, empty draft deposition (optionally with metadata
    already set). Returns the deposition dict (its ``id`` and
    ``links.bucket`` are what :func:`upload_artifact_files` needs).
    """
    session = session or requests.Session()
    response = session.post(
        api_url,
        json={"metadata": metadata} if metadata else {},
        headers=_headers(token),
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()


def update_deposition_metadata(
    deposit_id: int | str,
    metadata: dict,
    token: str,
    *,
    api_url: str = DEFAULT_API_URL,
    session: requests.Session | None = None,
    timeout: int = 30,
) -> dict:
    """Replace a draft deposition's metadata (title, description,
    etc.). Use to correct a previously-published draft's description
    after resizing its PCA artifact, or to set metadata on a freshly
    created deposition.
    """
    session = session or requests.Session()
    response = session.put(
        f"{api_url}/{deposit_id}",
        json={"metadata": metadata},
        headers=_headers(token),
        timeout=timeout,
    )
    response.raise_for_status()
    return response.json()


def default_metadata(
    project: str,
    tissue_name: str,
    *,
    n_components: int,
    cumulative_variance: float,
    related_identifiers: list[dict] | None = None,
) -> dict:
    """Build the standard metadata dict for a covariate-matrix
    deposition, following the shape of the first published records
    (COAD, BRCA) -- same creator/license/keywords template, with the
    PCA component count and cohort name filled in dynamically so the
    description text can't drift out of sync with the actual
    artifact the way the original COAD/BRCA records' "20 components"
    text did after the pipeline's own default moved on.

    Parameters
    ----------
    project : str
        TCGA study code, e.g. "COAD".
    tissue_name : str
        Full study name, e.g. "colon adenocarcinoma".
    n_components : int
        The PCA artifact's selected component count
        (``build_pca_artifact``'s ``n_components_selected``).
    cumulative_variance : float
        The cumulative explained-variance ratio at that component
        count (0-1), e.g. 0.99.
    related_identifiers : list[dict] | None
        Extra Zenodo ``related_identifiers`` entries (e.g. an
        ``isPartOf`` link to another record in the same series) on
        top of the standard GitHub-repo ``isSupplementTo`` entry.
    """
    project = project.upper()
    identifiers = [
        {
            "identifier": "https://github.com/alfaromurillo/sigmutselcovs",
            "relation": "isSupplementTo",
            "resource_type": "software",
            "scheme": "url",
        },
        *(related_identifiers or []),
    ]
    return {
        "title": f"sigmutselcovs covariate matrices: TCGA-{project}",
        "upload_type": "dataset",
        "description": (
            f"<p>Pre-built per-gene covariate matrices for "
            f"TCGA-{project} ({tissue_name}), built with "
            f'<a href="https://github.com/alfaromurillo/sigmutselcovs">'
            f"sigmutselcovs</a>: GTEx tissue-matched gene expression, "
            f"TCGA per-sample gene expression, replication timing "
            f"(mean replication time and CLR-transformed S-phase "
            f"fractions), Roadmap Epigenomics histone marks, and TCGA "
            f"ATAC-seq accessibility, indexed by Ensembl gene ID.</p>\n\n"
            f"<p>Files:</p>\n<ul>\n"
            f"<li><b>cov_matrix_full.parquet</b> &mdash; the full "
            f"matrix (all covariate columns)</li>\n"
            f"<li><b>cov_matrix_pca.parquet</b> &mdash; PCA-reduced "
            f"to {n_components} components, the smallest number "
            f"reaching {cumulative_variance:.0%} cumulative explained "
            f"variance (see pca_manifest.json)</li>\n"
            f"<li><b>cov_matrix_columns.csv</b> &mdash; a data "
            f"dictionary: one row per column of the full matrix, "
            f"with source, description, assembly, and units</li>\n"
            f"<li><b>build_manifest.json</b>, <b>pca_manifest.json</b> "
            f"&mdash; build provenance and PCA parameters</li>\n"
            f"</ul>\n\n"
            f"<p>Full source citations and methodology: "
            f'<a href="https://github.com/alfaromurillo/sigmutselcovs/'
            f'blob/master/SOURCES.md">SOURCES.md</a>.</p>'
        ),
        "access_right": "open",
        "creators": [
            {
                "name": "Alfaro-Murillo, Jorge A.",
                "affiliation": "University of Costa Rica",
                "orcid": "0000-0002-0481-0161",
            }
        ],
        "keywords": [
            "sigmutselcovs",
            "TCGA",
            f"TCGA-{project}",
            "covariate matrix",
            "mutation rate",
            "replication timing",
            "chromatin accessibility",
            "gene expression",
        ],
        "related_identifiers": identifiers,
        "license": "cc-by-4.0",
        "imprint_publisher": "Zenodo",
    }


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
