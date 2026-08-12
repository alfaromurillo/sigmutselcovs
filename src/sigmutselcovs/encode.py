"""ENCODE portal file resolution.

Resolves an ENCODE file accession (e.g. ``ENCFF001GSV``) to its
download URL and integrity metadata via the portal's JSON API.

BigWigs must always be downloaded to disk before opening —
``pyBigWig.open`` fails on both the portal ``@@download`` URLs and
the S3 URLs they redirect to.
"""

import logging

import requests

logger = logging.getLogger(__name__)

ENCODE_API = "https://www.encodeproject.org"


def resolve_encode_file(accession: str,
                        *,
                        session: requests.Session | None = None,
                        timeout: int = 60) -> dict:
    """Return metadata for an ENCODE file accession.

    Returns a dict with ``accession``, ``file_format``,
    ``output_type``, ``assembly``, ``md5sum``, ``file_size`` and
    ``url`` (direct S3 URL when the portal exposes it, otherwise the
    portal ``@@download`` URL, which redirects there).
    """
    session = session or requests.Session()
    response = session.get(
        f"{ENCODE_API}/files/{accession}/?format=json",
        timeout=timeout)
    response.raise_for_status()
    meta = response.json()
    url = (meta.get("cloud_metadata") or {}).get("url")
    if not url:
        url = ENCODE_API + meta["href"]
    return {
        "accession": accession,
        "file_format": meta.get("file_format"),
        "output_type": meta.get("output_type"),
        "assembly": meta.get("assembly"),
        "md5sum": meta.get("md5sum"),
        "file_size": meta.get("file_size"),
        "url": url,
    }
