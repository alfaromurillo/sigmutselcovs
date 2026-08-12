"""GDC (Genomic Data Commons) files-API client.

Replaces the checked-in per-project ``gdc_manifest`` and
``gdc_sample_sheet`` files with live queries: `query_gdc_files`
fetches the file inventory for a project, `write_gdc_manifest` and
`write_gdc_sample_sheet` reproduce the exact formats the loaders and
``gdc-client`` expect, and `download_gdc_files` fetches the files
(via ``gdc-client`` when available, plain HTTPS otherwise).
"""

import hashlib
import logging
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

import pandas as pd
import requests

logger = logging.getLogger(__name__)

GDC_API = "https://api.gdc.cancer.gov"

_DEFAULT_FIELDS = (
    "file_id",
    "file_name",
    "md5sum",
    "file_size",
    "state",
    "data_category",
    "data_type",
    "cases.project.project_id",
    "cases.submitter_id",
    "cases.samples.submitter_id",
    "cases.samples.tissue_type",
    "cases.samples.tumor_descriptor",
    "cases.samples.specimen_type",
    "cases.samples.preservation_method",
)

_SAMPLE_SHEET_COLUMNS = [
    "File ID", "File Name", "Data Category", "Data Type", "Project ID",
    "Case ID", "Sample ID", "Tissue Type", "Tumor Descriptor",
    "Specimen Type", "Preservation Method"]


def build_files_filter(project_id: str,
                       *,
                       data_type: str,
                       workflow_type: str | None = None,
                       access: str = "open",
                       extra: list[dict] | None = None) -> dict:
    """Build the GDC files-endpoint filter for a project query."""
    clauses = [
        {"op": "in", "content": {
            "field": "cases.project.project_id",
            "value": [project_id]}},
        {"op": "in", "content": {
            "field": "data_type", "value": [data_type]}},
        {"op": "in", "content": {
            "field": "access", "value": [access]}},
    ]
    if workflow_type is not None:
        clauses.append({"op": "in", "content": {
            "field": "analysis.workflow_type",
            "value": [workflow_type]}})
    clauses.extend(extra or [])
    return {"op": "and", "content": clauses}


def query_gdc_files(project_id: str,
                    *,
                    data_type: str = "Gene Expression Quantification",
                    workflow_type: str | None = "STAR - Counts",
                    access: str = "open",
                    fields: tuple[str, ...] = _DEFAULT_FIELDS,
                    page_size: int = 1000,
                    session: requests.Session | None = None,
                    timeout: int = 60) -> list[dict]:
    """Return all file hits for a project, sorted by file_id.

    Paginates through the GDC files endpoint; the sort makes the
    result deterministic across runs (and diffable against the
    historical checked-in manifests).
    """
    session = session or requests.Session()
    filters = build_files_filter(project_id,
                                 data_type=data_type,
                                 workflow_type=workflow_type,
                                 access=access)
    hits: list[dict] = []
    start = 0
    while True:
        payload = {
            "filters": filters,
            "fields": ",".join(fields),
            "size": page_size,
            "from": start,
        }
        response = session.post(f"{GDC_API}/files", json=payload,
                                timeout=timeout)
        response.raise_for_status()
        data = response.json()["data"]
        hits.extend(data["hits"])
        pagination = data["pagination"]
        start += pagination["count"]
        if start >= pagination["total"] or pagination["count"] == 0:
            break
    hits.sort(key=lambda h: h["file_id"])
    logger.info("GDC query for %s: %d files", project_id, len(hits))
    return hits


def write_gdc_manifest(hits: list[dict], path: str | Path) -> Path:
    """Write a gdc-client compatible manifest (id filename md5 size state)."""
    path = Path(path)
    rows = [(h["file_id"], h["file_name"], h["md5sum"],
             str(h["file_size"]), h.get("state", "released"))
            for h in hits]
    lines = ["\t".join(("id", "filename", "md5", "size", "state"))]
    lines += ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    logger.info("Wrote GDC manifest (%d files) to %s", len(rows), path)
    return path


def _joined(values: list[str]) -> str:
    """Join multi-valued fields the way GDC sample sheets do."""
    return ", ".join(dict.fromkeys(v for v in values if v))


def _sample_sheet_row(hit: dict) -> list[str]:
    cases = hit.get("cases", [])
    samples = [s for c in cases for s in c.get("samples", [])]

    def sample_field(name: str) -> str:
        return _joined([s.get(name) or "Unknown" for s in samples])

    return [
        hit["file_id"],
        hit["file_name"],
        hit.get("data_category", ""),
        hit.get("data_type", ""),
        _joined([c.get("project", {}).get("project_id", "")
                 for c in cases]),
        _joined([c.get("submitter_id", "") for c in cases]),
        _joined([s.get("submitter_id", "") for s in samples]),
        sample_field("tissue_type"),
        sample_field("tumor_descriptor"),
        sample_field("specimen_type"),
        sample_field("preservation_method"),
    ]


def write_gdc_sample_sheet(hits: list[dict],
                           path: str | Path | None = None,
                           *,
                           directory: str | Path | None = None) -> Path:
    """Write a GDC sample sheet (the 11 columns the loaders parse).

    Either give an explicit ``path`` or a ``directory``, in which
    case the file is named ``gdc_sample_sheet.<YYYY-MM-DD>.tsv`` —
    date-suffixed names sort correctly for the loaders' newest-file
    glob.
    """
    if path is None:
        if directory is None:
            raise ValueError("Give either path or directory")
        path = (Path(directory)
                / f"gdc_sample_sheet.{date.today().isoformat()}.tsv")
    path = Path(path)
    frame = pd.DataFrame([_sample_sheet_row(h) for h in hits],
                         columns=_SAMPLE_SHEET_COLUMNS)
    frame.to_csv(path, sep="\t", index=False)
    logger.info("Wrote GDC sample sheet (%d rows) to %s",
                len(frame), path)
    return path


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_one(hit: dict, dest: Path, *, session: requests.Session,
                  verify_md5: bool, timeout: int) -> Path:
    """Fetch one GDC file into <dest>/<file_id>/<file_name>."""
    target_dir = dest / hit["file_id"]
    target = target_dir / hit["file_name"]
    if target.exists() and target.stat().st_size == hit["file_size"]:
        return target
    target_dir.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(target.suffix + ".part")
    with session.get(f"{GDC_API}/data/{hit['file_id']}",
                     stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with open(part, "wb") as fh:
            for chunk in response.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    if verify_md5 and _md5(part) != hit["md5sum"]:
        part.unlink()
        raise OSError(f"md5 mismatch for {hit['file_id']} "
                      f"({hit['file_name']})")
    part.replace(target)
    return target


def download_gdc_files(hits: list[dict],
                       dest: str | Path,
                       *,
                       use_gdc_client: bool | None = None,
                       manifest_path: str | Path | None = None,
                       workers: int = 6,
                       verify_md5: bool = True,
                       timeout: int = 600) -> list[Path]:
    """Download GDC files into ``<dest>/<file_id>/<file_name>``.

    That layout is required: `import_tcga_gene_expression` recovers
    the file id from the parent directory name.

    Parameters
    ----------
    use_gdc_client : bool | None
        None (default) auto-detects ``gdc-client`` on PATH and uses
        it when present (resumable, parallel); otherwise plain HTTPS
        streaming with md5 verification.
    manifest_path : str | Path | None
        Manifest to hand to ``gdc-client -m``; written to a temporary
        file inside ``dest`` when not given.
    """
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)

    pending = [h for h in hits
               if not ((dest / h["file_id"] / h["file_name"]).exists()
                       and (dest / h["file_id"] / h["file_name"])
                       .stat().st_size == h["file_size"])]
    logger.info("GDC download: %d of %d files to fetch",
                len(pending), len(hits))
    if not pending:
        return [dest / h["file_id"] / h["file_name"] for h in hits]

    if use_gdc_client is None:
        use_gdc_client = shutil.which("gdc-client") is not None

    if use_gdc_client:
        if manifest_path is None:
            manifest_path = dest / ".gdc_manifest.pending.tsv"
            write_gdc_manifest(pending, manifest_path)
        subprocess.run(
            ["gdc-client", "download",
             "-m", str(manifest_path), "-d", str(dest)],
            check=True)
    else:
        session = requests.Session()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(_download_one, h, dest, session=session,
                            verify_md5=verify_md5, timeout=timeout)
                for h in pending]
            for future in futures:
                future.result()

    return [dest / h["file_id"] / h["file_name"] for h in hits]
