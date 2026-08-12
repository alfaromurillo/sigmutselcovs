"""Download layer: fetch covariate source data for a project.

`download_covariates` orchestrates the per-source downloaders
(GDC gene expression, replication timing, Roadmap chromatin, TCGA
ATAC-seq) into a project data directory laid out per `ProjectPaths`.
All downloads are idempotent: present files are skipped, partial
files are resumed or replaced atomically.
"""

import hashlib
import logging
import os
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

_CHUNK = 1 << 20  # 1 MiB


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str,
                  dest: str | Path,
                  *,
                  expected_md5: str | None = None,
                  expected_size: int | None = None,
                  force: bool = False,
                  resume: bool = True,
                  session: requests.Session | None = None,
                  timeout: int = 600) -> Path:
    """Download a file idempotently and atomically.

    - Skips when ``dest`` exists (and matches ``expected_size`` when
      given) unless ``force``.
    - Streams to ``<dest>.part`` and renames into place, so an
      interrupted download never leaves a truncated ``dest`` (the
      old ``wget -O`` scripts did).
    - Resumes an existing ``.part`` with an HTTP Range request when
      ``resume`` (falls back to a full fetch if the server ignores
      it).
    - Verifies ``expected_md5`` when given.
    """
    dest = Path(dest)
    if dest.exists() and not force:
        if expected_size is None or dest.stat().st_size == expected_size:
            logger.info("Already present, skipping: %s", dest.name)
            return dest
        logger.warning("Size mismatch for %s (%d != %d); refetching",
                       dest.name, dest.stat().st_size, expected_size)
    dest.parent.mkdir(parents=True, exist_ok=True)

    session = session or requests.Session()
    part = dest.with_suffix(dest.suffix + ".part")
    headers = {}
    mode = "wb"
    if resume and part.exists() and part.stat().st_size > 0:
        headers["Range"] = f"bytes={part.stat().st_size}-"
        mode = "ab"

    with session.get(url, stream=True, headers=headers,
                     timeout=timeout) as response:
        if headers and response.status_code == 200:
            # Server ignored the Range request; start over.
            mode = "wb"
        elif headers and response.status_code == 206:
            logger.info("Resuming %s at %d bytes",
                        dest.name, part.stat().st_size)
        response.raise_for_status()
        with open(part, mode) as fh:
            for chunk in response.iter_content(chunk_size=_CHUNK):
                fh.write(chunk)

    if expected_size is not None and part.stat().st_size != expected_size:
        raise OSError(
            f"Downloaded size {part.stat().st_size} != expected "
            f"{expected_size} for {url}")
    if expected_md5 is not None:
        got = _md5(part)
        if got != expected_md5:
            part.unlink()
            raise OSError(f"md5 mismatch for {url}: {got} != "
                          f"{expected_md5}")
    os.replace(part, dest)
    logger.info("Downloaded %s (%d bytes)", dest.name,
                dest.stat().st_size)
    return dest
