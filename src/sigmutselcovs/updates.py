"""Check external covariate sources for updates.

Covariate data sources (GTEx releases, GDC data releases, ENCODE
files, Roadmap, GEO) change on their own schedules.  `check_updates`
compares each source's current remote state against the values
recorded in ``data/sources.json`` and reports what changed; run it
every ~6 months.

Never raises on network problems — unreachable sources are reported
as such.
"""

import json
import logging
from datetime import date
from pathlib import Path

import pandas as pd
import requests

from gdcfetch import get_data_size, search_files

from .covariates_locations import location_covariates_data

logger = logging.getLogger(__name__)

location_sources = location_covariates_data / "sources.json"


def _check_http_head(entry: dict, session, timeout: int) -> dict:
    url = entry.get("url") or entry["check"].get("url")
    response = session.head(
        url, timeout=timeout, allow_redirects=True
    )
    response.raise_for_status()
    return {
        key: response.headers.get(key)
        for key in entry["check"]["compare"]
    }


def _check_gdc_file_count(entry: dict, session, timeout: int) -> dict:
    return {
        project: len(
            search_files(project, session=session, timeout=timeout)
        )
        for project in entry["check"]["projects"]
    }


def _check_gdc_blob_size(entry: dict, session, timeout: int) -> dict:
    """Size-only check for GDC /data/<uuid> blobs.

    These UUIDs (e.g. the TCGA ATAC-seq tarballs) are not indexed in
    /files at all -- GET /files/<uuid> 404s for them even though the
    blob is real (verified 2026-08). A 1-byte ranged GET on
    /data/<uuid> is the only reliable way to get metadata, and it
    only gives size (the Content-MD5 header on a ranged response
    covers just the requested byte range, not the whole file, so
    md5 isn't checkable this way).
    """
    return {
        label: {
            "file_size": get_data_size(
                uuid, session=session, timeout=timeout
            )
        }
        for label, uuid in entry["check"]["uuids"].items()
    }


def _check_encode_files(entry: dict, session, timeout: int) -> dict:
    from .encode import resolve_encode_file

    current = {}
    for accession in entry["check"]["accessions"]:
        meta = resolve_encode_file(
            accession, session=session, timeout=timeout
        )
        current[accession] = {
            field: meta.get(field)
            for field in entry["check"]["compare"]
        }
    return current


_METHODS = {
    "http_head": _check_http_head,
    "gdc_file_count": _check_gdc_file_count,
    "gdc_blob_size": _check_gdc_blob_size,
    "encode_files": _check_encode_files,
}


def check_updates(
    *,
    sources: list[str] | None = None,
    sources_path: str | Path | None = None,
    update_file: bool = False,
    session: requests.Session | None = None,
    timeout: int = 30,
) -> pd.DataFrame:
    """Compare covariate sources against their recorded state.

    Parameters
    ----------
    sources : list[str] | None
        Restrict to these source names; None checks all.
    sources_path : str | Path | None
        Alternative sources.json; defaults to the packaged one.
    update_file : bool
        Rewrite the ``known`` blocks (and the ``checked`` date) with
        the current values, so the diff can be reviewed and
        committed.  Only useful on a writable checkout.
    timeout : int
        Per-request timeout in seconds.

    Returns
    -------
    pd.DataFrame
        Columns: source, status (ok|changed|unknown|unreachable),
        current, known, notes.  ``unknown`` means there was no
        recorded state yet (first run).
    """
    sources_path = (
        Path(sources_path)
        if sources_path is not None
        else location_sources
    )
    raw = json.loads(sources_path.read_text())
    session = session or requests.Session()

    rows = []
    for name, entry in raw["sources"].items():
        if sources is not None and name not in sources:
            continue
        method = entry["check"]["method"]
        try:
            current = _METHODS[method](entry, session, timeout)
        except KeyError:
            raise ValueError(
                f"Unknown check method {method!r} "
                f"for source {name}"
            )
        except (
            Exception
        ) as exc:  # noqa: BLE001 - network must not abort
            rows.append(
                {
                    "source": name,
                    "status": "unreachable",
                    "current": None,
                    "known": entry.get("known"),
                    "notes": str(exc),
                }
            )
            continue
        known = entry.get("known") or {}
        if not known:
            status = "unknown"
        elif current == known:
            status = "ok"
        else:
            status = "changed"
        rows.append(
            {
                "source": name,
                "status": status,
                "current": current,
                "known": known or None,
                "notes": entry.get("notes", ""),
            }
        )
        if update_file:
            entry["known"] = current

    if update_file:
        raw["checked"] = date.today().isoformat()
        sources_path.write_text(json.dumps(raw, indent=2) + "\n")
        logger.info("Updated known states in %s", sources_path)

    frame = pd.DataFrame(
        rows,
        columns=["source", "status", "current", "known", "notes"],
    )
    for _, row in frame.iterrows():
        level = (
            logging.WARNING
            if row["status"] in ("changed", "unreachable")
            else logging.INFO
        )
        logger.log(
            level,
            "check_updates %s: %s",
            row["source"],
            row["status"],
        )
    return frame


def print_update_report(frame: pd.DataFrame) -> None:
    """Print the update report in a compact form."""
    for _, row in frame.iterrows():
        print(f"[{row['status']:11s}] {row['source']}")
        if row["status"] == "changed":
            print(f"    known:   {row['known']}")
            print(f"    current: {row['current']}")
        if row["notes"]:
            print(f"    note: {row['notes']}")
