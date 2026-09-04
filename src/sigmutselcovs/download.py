"""Download layer: fetch covariate source data for a project.

`download_covariates` orchestrates the per-source downloaders
(GDC gene expression, replication timing, Roadmap chromatin, TCGA
ATAC-seq) into a project data directory laid out per `ProjectPaths`.
All downloads are idempotent: present files are skipped, partial
files are resumed or replaced atomically.
"""

import gzip
import hashlib
import logging
import os
import shutil
import tarfile
from dataclasses import dataclass, field
from pathlib import Path

import requests
from gdcfetch import (
    download_files,
    search_files,
    write_manifest,
    write_sample_sheet,
)

from .encode import resolve_encode_file
from .paths import ProjectPaths, bigwig_files, project_paths
from .registry import (
    AtacSpec,
    EncodeChromatinSpec,
    GexpSpec,
    RepliseqSpec,
    RoadmapSpec,
    get_project,
)

GDC_DATA_URL = "https://api.gdc.cancer.gov/data"

logger = logging.getLogger(__name__)

_CHUNK = 1 << 20  # 1 MiB


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(
    url: str,
    dest: str | Path,
    *,
    expected_md5: str | None = None,
    expected_size: int | None = None,
    force: bool = False,
    resume: bool = True,
    session: requests.Session | None = None,
    timeout: int = 600,
) -> Path:
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
        if (
            expected_size is None
            or dest.stat().st_size == expected_size
        ):
            logger.info("Already present, skipping: %s", dest.name)
            return dest
        logger.warning(
            "Size mismatch for %s (%d != %d); refetching",
            dest.name,
            dest.stat().st_size,
            expected_size,
        )
    dest.parent.mkdir(parents=True, exist_ok=True)

    session = session or requests.Session()
    part = dest.with_suffix(dest.suffix + ".part")
    headers = {}
    mode = "wb"
    if resume and part.exists() and part.stat().st_size > 0:
        headers["Range"] = f"bytes={part.stat().st_size}-"
        mode = "ab"

    with session.get(
        url, stream=True, headers=headers, timeout=timeout
    ) as response:
        if headers and response.status_code == 200:
            # Server ignored the Range request; start over.
            mode = "wb"
        elif headers and response.status_code == 206:
            logger.info(
                "Resuming %s at %d bytes",
                dest.name,
                part.stat().st_size,
            )
        response.raise_for_status()
        with open(part, mode) as fh:
            fh.writelines(response.iter_content(chunk_size=_CHUNK))

    if (
        expected_size is not None
        and part.stat().st_size != expected_size
    ):
        raise OSError(
            f"Downloaded size {part.stat().st_size} != expected "
            f"{expected_size} for {url}"
        )
    if expected_md5 is not None:
        got = _md5(part)
        if got != expected_md5:
            part.unlink()
            raise OSError(
                f"md5 mismatch for {url}: {got} != " f"{expected_md5}"
            )
    os.replace(part, dest)
    logger.info(
        "Downloaded %s (%d bytes)", dest.name, dest.stat().st_size
    )
    return dest


def _gtex_source() -> dict:
    """GTEx GCT url/version from the packaged sources.json."""
    import json

    from .covariates_locations import location_covariates_data

    raw = json.loads(
        (location_covariates_data / "sources.json").read_text()
    )
    return raw["sources"]["gtex_gene_median_tpm"]


def _gtex_cache_dir() -> Path:
    base = os.environ.get(
        "XDG_CACHE_HOME", str(Path.home() / ".cache")
    )
    return Path(base) / "sigmutselcovs" / "gtex"


def resolve_gtex_gct(explicit: str | Path | None = None) -> Path:
    """Resolve the GTEx GCT path without downloading.

    Order: explicit path, then any ``*.gct`` in the packaged data
    directory, then any ``*.gct`` in the user cache
    (``$XDG_CACHE_HOME/sigmutselcovs/gtex/``).  When none exists,
    returns the cache path the GCT *would* have, so callers get a
    meaningful FileNotFoundError (or can hand it to
    `ensure_gtex_gct`).
    """
    from .covariates_locations import (
        location_cov_gene_expression_gtex,
    )

    if explicit is not None:
        return Path(explicit)
    if location_cov_gene_expression_gtex.exists():
        return location_cov_gene_expression_gtex
    packaged = sorted(
        location_cov_gene_expression_gtex.parent.glob("*.gct")
    )
    if packaged:
        return packaged[0]
    cached = sorted(_gtex_cache_dir().glob("*.gct"))
    if cached:
        return cached[0]
    return _gtex_cache_dir() / location_cov_gene_expression_gtex.name


def ensure_gtex_gct(
    dest: str | Path | None = None,
    *,
    force: bool = False,
    session: requests.Session | None = None,
) -> Path:
    """Return a usable GTEx GCT path, downloading it if needed.

    Downloads go to the user cache (or ``dest`` when given) — never
    into site-packages, unlike the old package setup.sh, which broke
    read-only installs.
    """
    resolved = resolve_gtex_gct(dest)
    if resolved.exists() and not force:
        return resolved
    source = _gtex_source()
    gz = resolved.with_suffix(resolved.suffix + ".gz")
    logger.info(
        "Fetching GTEx %s GCT to %s",
        source.get("version", "?"),
        resolved,
    )
    download_file(source["url"], gz, force=force, session=session)
    with gzip.open(gz, "rb") as src, open(resolved, "wb") as out:
        shutil.copyfileobj(src, out)
    gz.unlink()
    return resolved


_GENCODE_SOURCE_KEYS = {"hg19": "gencode_v19", "hg38": "gencode_v38"}


def _gencode_source(assembly: str) -> dict:
    """GENCODE GTF url/version from the packaged sources.json."""
    import json

    from .covariates_locations import location_covariates_data

    if assembly not in _GENCODE_SOURCE_KEYS:
        raise ValueError(
            f"Unknown assembly {assembly!r}; expected one of "
            f"{sorted(_GENCODE_SOURCE_KEYS)}"
        )
    raw = json.loads(
        (location_covariates_data / "sources.json").read_text()
    )
    return raw["sources"][_GENCODE_SOURCE_KEYS[assembly]]


def _gencode_data_home() -> Path:
    """Shared, tool-agnostic GENCODE location -- not namespaced under
    this package. GENCODE GTFs are a universal reference (sigmutsel
    uses them directly too, independent of covariates), not data
    specific to sigmutselcovs, so this deliberately isn't
    ``.../sigmutselcovs/gencode``. ``GENCODE_DATA_HOME`` overrides;
    otherwise ``$XDG_DATA_HOME/gencode`` (data, not cache -- these
    files are worth keeping, not disposable the way an HTTP cache
    is; matches genomepy's choice of XDG_DATA_HOME for the same
    reason). sigmutsel's own downloader checks the same location
    (see its setup.py), so whichever package downloads first, the
    other reuses it instead of fetching its own copy.
    """
    override = os.environ.get("GENCODE_DATA_HOME")
    if override:
        return Path(override)
    xdg_data_home = os.environ.get(
        "XDG_DATA_HOME", str(Path.home() / ".local" / "share")
    )
    return Path(xdg_data_home) / "gencode"


def resolve_gencode_gtf(
    assembly: str, explicit: str | Path | None = None
) -> Path:
    """Resolve a GENCODE GTF path without downloading.

    Order: explicit path, then this package's own data directory,
    then the shared external location (``GENCODE_DATA_HOME`` or
    ``$XDG_DATA_HOME/gencode``, see `_gencode_data_home`). When none
    exists, returns the path the GTF *would* have there, so callers
    get a meaningful FileNotFoundError (or can hand it to
    `ensure_gencode_gtf`).
    """
    from .covariates_locations import (
        location_gencode19_gtf,
        location_gencode38_gtf,
    )

    packaged = {
        "hg19": location_gencode19_gtf,
        "hg38": location_gencode38_gtf,
    }[assembly]
    if explicit is not None:
        return Path(explicit)
    if packaged.exists():
        return packaged
    return _gencode_data_home() / packaged.name


def ensure_gencode_gtf(
    assembly: str,
    dest: str | Path | None = None,
    *,
    force: bool = False,
    session: requests.Session | None = None,
) -> Path:
    """Return a usable GENCODE GTF path, downloading it if needed.

    ``assembly`` is 'hg19' (GENCODE v19) or 'hg38' (GENCODE v38) --
    the same two builds `build_covariate_matrix` already dispatches
    on per registry source. Downloads go to the shared external
    location (or ``dest`` when given), never into site-packages.
    Stored gzipped, same as GTF loading elsewhere in this package
    already handles gzip transparently.
    """
    resolved = resolve_gencode_gtf(assembly, dest)
    if resolved.exists() and not force:
        return resolved
    source = _gencode_source(assembly)
    logger.info(
        "Fetching GENCODE %s GTF to %s",
        source.get("version", "?"),
        resolved,
    )
    download_file(
        source["url"], resolved, force=force, session=session
    )
    return resolved


def download_gdc_gene_expression(
    spec: GexpSpec,
    paths: ProjectPaths,
    *,
    force: bool = False,
    session: requests.Session | None = None,
) -> int:
    """Query GDC and download a project's STAR count files.

    Writes a fresh manifest and sample sheet into the gene
    expression directory (replacing the need for checked-in ones)
    and downloads into ``star_gene_counts/<file_id>/<file_name>``.
    Returns the number of files in the inventory.
    """
    paths.gexp_tcga_dir.mkdir(parents=True, exist_ok=True)
    hits = search_files(
        spec.tcga_project_id,
        data_type=spec.data_type,
        workflow_type=spec.workflow_type,
        session=session,
    )
    if not hits:
        raise ValueError(
            "GDC returned no files for " f"{spec.tcga_project_id}"
        )
    from datetime import UTC, datetime

    today = datetime.now(tz=UTC).date().isoformat()
    manifest = paths.gexp_tcga_dir / f"gdc_manifest.{today}.txt"
    write_manifest(hits, manifest)
    write_sample_sheet(hits, directory=paths.gexp_tcga_dir)
    download_files(hits, paths.gexp_star_dir, session=session)
    return len(hits)


def download_repliseq(
    spec: RepliseqSpec,
    paths: ProjectPaths,
    *,
    force: bool = False,
    session: requests.Session | None = None,
) -> list[Path]:
    """Download the replication-timing source files for a project.

    - ``mat``: fetch the gzipped GEO table and gunzip it next to the
      cache files (skipped when the .mat is already present).
    - ``fraction_bigwigs`` / ``wavelet``: resolve each ENCODE track
      accession and download the bigWig to
      ``replication_timing/encode/<ACCESSION>.bigWig``, verifying
      the portal's md5.
    """
    session = session or requests.Session()
    if spec.type == "mat":
        target = paths.rt_dir / spec.filename
        if target.exists() and not force:
            logger.info(
                "Repli-seq MAT already present: %s", target.name
            )
            return [target]
        if spec.url is None:
            raise ValueError(
                "Registry gives no URL for the "
                f"{spec.cell_line} MAT file"
            )
        gz = target.with_suffix(target.suffix + ".gz")
        download_file(spec.url, gz, force=force, session=session)
        with gzip.open(gz, "rb") as src, open(target, "wb") as out:
            shutil.copyfileobj(src, out)
        gz.unlink()
        logger.info("Decompressed %s", target.name)
        return [target]

    if spec.type in ("fraction_bigwigs", "wavelet"):
        out: list[Path] = []
        for track in spec.tracks:
            dest = paths.rt_encode_dir / f"{track.accession}.bigWig"
            if dest.exists() and not force:
                logger.info("Track already present: %s", dest.name)
                out.append(dest)
                continue
            if track.url is not None:
                url, md5, size = track.url, None, None
            else:
                meta = resolve_encode_file(
                    track.accession, session=session
                )
                url = meta["url"]
                md5 = meta["md5sum"]
                size = meta["file_size"]
                if meta.get("assembly") not in (None, spec.assembly):
                    logger.warning(
                        "ENCODE %s assembly %s != registry %s",
                        track.accession,
                        meta["assembly"],
                        spec.assembly,
                    )
            out.append(
                download_file(
                    url,
                    dest,
                    expected_md5=md5,
                    expected_size=size,
                    force=force,
                    session=session,
                )
            )
        return out

    raise ValueError(f"Unknown repliseq type {spec.type!r}")


def download_tcga_atac(
    spec: AtacSpec,
    paths: ProjectPaths,
    *,
    force: bool = False,
    session: requests.Session | None = None,
) -> list[Path]:
    """Download and unpack the TCGA ATAC-seq bigWig tarball.

    The Corces et al. archives nest the bigWigs under a deep
    Stanford path; members are extracted flattened into the ATAC
    directory, keeping only ``*.bw``/``*.bigWig`` basenames and
    refusing absolute or parent-traversing names.  The tarball is
    removed on success.  Skipped entirely when bigWigs are already
    present.
    """
    existing = bigwig_files(paths.atac_dir)
    if existing and not force:
        logger.info(
            "ATAC bigWigs already present (%d files); skipping",
            len(existing),
        )
        return existing

    tarball = (
        paths.atac_dir / f"{spec.column_prefix.upper()}_bigWigs.tgz"
    )
    download_file(
        f"{GDC_DATA_URL}/{spec.gdc_uuid}",
        tarball,
        force=force,
        session=session,
    )

    extracted: list[Path] = []
    with tarfile.open(tarball, "r:gz") as tar:
        for member in tar:
            if not member.isfile():
                continue
            name = Path(member.name)
            if name.is_absolute() or ".." in name.parts:
                logger.warning(
                    "Refusing suspicious tar member: %s", member.name
                )
                continue
            if name.suffix.lower() not in (".bw", ".bigwig"):
                continue
            target = paths.atac_dir / name.name
            source = tar.extractfile(member)
            if source is None:
                continue
            with source, open(target, "wb") as out:
                shutil.copyfileobj(source, out)
            extracted.append(target)
    tarball.unlink()
    if not extracted:
        raise OSError(
            "No bigWig members found in the ATAC tarball "
            f"for {spec.gdc_uuid}"
        )
    logger.info("Extracted %d ATAC bigWigs", len(extracted))
    return sorted(extracted)


@dataclass
class DownloadReport:
    """Per-source outcome of a download_covariates run."""

    project: str
    sources: dict[str, dict] = field(default_factory=dict)

    def add(self, source: str, status: str, **info) -> None:
        self.sources[source] = {"status": status, **info}

    def summary(self) -> str:
        lines = [f"Download report for {self.project}:"]
        for source, info in self.sources.items():
            extra = ", ".join(
                f"{k}={v}" for k, v in info.items() if k != "status"
            )
            lines.append(
                f"  {source:10s} {info['status']}"
                + (f" ({extra})" if extra else "")
            )
        return "\n".join(lines)


_WHICH = ("gexp", "repliseq", "roadmap", "atac", "encode_chromatin")


def download_covariates(
    project: str,
    data_dir: str | Path,
    *,
    which: tuple[str, ...] = _WHICH,
    registry_path: str | Path | None = None,
    force: bool = False,
    dry_run: bool = False,
    session: requests.Session | None = None,
) -> DownloadReport:
    """Download all covariate source data for a registered project.

    Parameters
    ----------
    project : str
        TCGA study code from the registry (e.g. 'COAD', 'BRCA').
    data_dir : str | Path
        Project data directory (created as needed) laid out per
        `ProjectPaths`.
    which : tuple[str, ...]
        Sources to fetch, any of gexp, repliseq, roadmap, atac,
        encode_chromatin.
    dry_run : bool
        Only report what would be fetched.

    Notes
    -----
    MAF files are not covariates and stay with
    ``sigmutsel.download_tcga_data``.
    """
    unknown = set(which) - set(_WHICH)
    if unknown:
        raise ValueError(
            f"Unknown sources: {sorted(unknown)}; " f"valid: {_WHICH}"
        )
    spec = get_project(project, registry_path)
    paths = project_paths(data_dir)
    report = DownloadReport(project=spec.code)

    for source in _WHICH:
        if source not in which:
            continue
        source_spec = getattr(spec, source)
        if source_spec is None:
            report.add(source, "not-in-registry")
            continue
        if dry_run:
            report.add(source, "would-download")
            continue
        try:
            if source == "gexp":
                n = download_gdc_gene_expression(
                    source_spec, paths, force=force, session=session
                )
                report.add(source, "ok", n_files=n)
            elif source == "repliseq":
                files = download_repliseq(
                    source_spec, paths, force=force, session=session
                )
                report.add(source, "ok", n_files=len(files))
            elif source == "roadmap":
                files = download_roadmap_tracks(
                    source_spec, paths, force=force, session=session
                )
                expected = len(source_spec.eids) * len(
                    source_spec.marks
                )
                report.add(
                    source,
                    "ok",
                    n_files=len(files),
                    expected=expected,
                )
            elif source == "atac":
                files = download_tcga_atac(
                    source_spec, paths, force=force, session=session
                )
                report.add(source, "ok", n_files=len(files))
            elif source == "encode_chromatin":
                files = download_encode_chromatin_tracks(
                    source_spec, paths, force=force, session=session
                )
                report.add(
                    source,
                    "ok",
                    n_files=len(files),
                    expected=len(source_spec.tracks),
                )
        # One source's failure shouldn't abort the others; recorded
        # in the report instead of raised.
        except Exception as exc:  # noqa: BLE001
            logger.error("Download failed for %s: %s", source, exc)
            report.add(source, "failed", error=str(exc))

    logger.info("%s", report.summary())
    return report


def download_roadmap_tracks(
    spec: RoadmapSpec,
    paths: ProjectPaths,
    *,
    force: bool = False,
    session: requests.Session | None = None,
) -> list[Path]:
    """Download Roadmap fold-change signal bigWigs (eids x marks).

    Missing marks are common (e.g. E027 has no H3K27ac) — a 404 is
    logged and skipped unless the mark is listed in
    ``required_marks``.
    """
    session = session or requests.Session()
    out: list[Path] = []
    for eid in spec.eids:
        for mark in spec.marks:
            name = f"{eid}-{mark}.fc.signal.bigwig"
            dest = paths.roadmap_dir / name
            if dest.exists() and not force:
                logger.info("Track already present: %s", name)
                out.append(dest)
                continue
            url = spec.url_template.format(eid=eid, mark=mark)
            try:
                out.append(
                    download_file(
                        url, dest, force=force, session=session
                    )
                )
            except Exception as exc:
                if mark in spec.required_marks:
                    raise
                logger.warning("Skipping %s: %s", name, exc)
    return out


def download_encode_chromatin_tracks(
    spec: EncodeChromatinSpec,
    paths: ProjectPaths,
    *,
    force: bool = False,
    session: requests.Session | None = None,
) -> list[Path]:
    """Download ENCODE-native histone/DNase bigWigs (one per track).

    Same accession -> bigWig resolution as the ``fraction_bigwigs``/
    ``wavelet`` branch of `download_repliseq` (portal metadata via
    `encode.resolve_encode_file`, md5-verified download), just to
    ``chromatin/encode/<ACCESSION>.bigWig`` instead of
    ``replication_timing/encode/``.
    """
    session = session or requests.Session()
    out: list[Path] = []
    for track in spec.tracks:
        # "_fc_signal_" is a naming convention, not a literal claim
        # about DNase-seq's signal representation -- it lets the
        # resulting column reuse validate.py's existing "fc_signal"
        # chromatin-signal marker and _mark_mean's substring lookup
        # by track.label (e.g. "h3k23ac", "dnase") the same way
        # Roadmap's "{eid}-{mark}.fc.signal.bigwig" naming already
        # does, without a parallel source-specific validation path.
        dest = (
            paths.encode_chromatin_dir
            / f"{track.label}_fc_signal_{track.accession}.bigWig"
        )
        if dest.exists() and not force:
            logger.info("Track already present: %s", dest.name)
            out.append(dest)
            continue
        if track.url is not None:
            url, md5, size = track.url, None, None
        else:
            meta = resolve_encode_file(
                track.accession, session=session
            )
            url = meta["url"]
            md5 = meta["md5sum"]
            size = meta["file_size"]
            if meta.get("assembly") not in (None, spec.assembly):
                logger.warning(
                    "ENCODE %s assembly %s != registry %s",
                    track.accession,
                    meta["assembly"],
                    spec.assembly,
                )
        out.append(
            download_file(
                url,
                dest,
                expected_md5=md5,
                expected_size=size,
                force=force,
                session=session,
            )
        )
    return out
