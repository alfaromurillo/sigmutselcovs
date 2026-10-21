"""Parallel, resumable chromatin cache builds that keep no bigWigs.

A large pool (every Roadmap consolidated epigenome and mark, hundreds
of ENCODE tracks) does not fit on disk at once: each track is
downloaded, summarised over the source's windows and transforms, and
deleted, several tracks at a time. Each finished track is written to
its own partial file, so an interrupted build resumes where it
stopped, and the partials are merged into the cache file that
`build_covariate_matrix` reads for the same options
(`covariates_chromatin.chromatin_cache_path`). With no bigWigs on
disk the build then reads that cache as it is.

Columns are named as a build from bigWigs on disk names them (the
file stem, e.g. ``e075_h3k4me3_fc_signal_cds_span_logmean``), in the
same order (file names sorted), so the two routes give one matrix.
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import requests

from .covariates_chromatin import (
    _gene_table,
    check_transforms,
    chromatin_cache_path,
    summarize_bigwig_windows,
)
from .download import download_file
from .encode import resolve_encode_file
from .paths import ProjectPaths, project_paths
from .registry import ProjectSpec, get_project

logger = logging.getLogger(__name__)

STREAM_SOURCES = ("roadmap", "encode_chromatin")


@dataclass(frozen=True)
class TrackJob:
    """One track to fetch: where from, where to, and whether a
    failure is fatal (a required Roadmap mark)."""

    url: str
    dest: Path
    md5: str | None = None
    size: int | None = None
    required: bool = False


def track_jobs(
    spec: ProjectSpec,
    source: str,
    paths: ProjectPaths,
    *,
    session: requests.Session | None = None,
) -> list[TrackJob]:
    """The tracks of one chromatin source, named as the downloaders
    name them (``<eid>-<mark>.fc.signal.bigwig``,
    ``<label>_fc_signal_<accession>.bigWig``)."""
    if source == "roadmap":
        rm = spec.roadmap
        return [
            TrackJob(
                url=rm.url_template.format(eid=eid, mark=mark),
                dest=paths.roadmap_dir
                / f"{eid}-{mark}.fc.signal.bigwig",
                required=mark in rm.required_marks,
            )
            for eid in rm.eids
            for mark in rm.marks
        ]
    if source == "encode_chromatin":
        session = session or requests.Session()
        jobs = []
        for track in spec.encode_chromatin.tracks:
            dest = (
                paths.encode_chromatin_dir
                / f"{track.label}_fc_signal_{track.accession}.bigWig"
            )
            if track.url is not None:
                jobs.append(TrackJob(url=track.url, dest=dest))
                continue
            meta = resolve_encode_file(
                track.accession, session=session
            )
            jobs.append(
                TrackJob(
                    url=meta["url"],
                    dest=dest,
                    md5=meta["md5sum"],
                    size=meta["file_size"],
                )
            )
        return jobs
    raise ValueError(
        f"source must be one of {STREAM_SOURCES}, not {source!r}"
    )


def _summaries(
    path: Path, genes, windows, transforms
) -> pd.DataFrame:
    frames = [
        summarize_bigwig_windows(
            path,
            genes,
            windows=windows,
            transform=transform,
            label=path.stem,
        )
        for transform in transforms
    ]
    return pd.concat(frames, axis=1)


def stream_chromatin_cache(
    project: str,
    data_dir: str | Path,
    source: str,
    *,
    registry_path: str | Path | None = None,
    gencode_gtfs: Mapping[str, str | Path] | None = None,
    chromatin_windows: Sequence[str] | None = None,
    workers: int = 6,
    keep_files: bool = False,
    jobs: Sequence[TrackJob] | None = None,
) -> Path:
    """Build one chromatin source's cache by streaming its tracks.

    Parameters
    ----------
    project : str
        Registered project code.
    data_dir : str | Path
        Project data directory (`ProjectPaths` layout).
    source : {"roadmap", "encode_chromatin"}
        The chromatin source to build.
    registry_path, gencode_gtfs
        As for `build_covariate_matrix`.
    chromatin_windows : Sequence[str] | None
        Windows replacing the registry's, as
        `build_covariate_matrix` takes them (the cache then matches a
        build with the same override).
    workers : int, default 6
        Tracks downloaded and summarised at once (disk holds at most
        this many bigWigs).
    keep_files : bool, default False
        Keep each bigWig after summarising it.
    jobs : Sequence[TrackJob] | None
        The tracks, if not the registry's (`track_jobs`).

    Returns
    -------
    pathlib.Path
        The cache file written. A track that fails to download is
        skipped with a warning, unless it is a required Roadmap mark.
    """
    from .builder import (
        _default_gencode_gtfs,
        _with_chromatin_windows,
    )

    spec = get_project(project, registry_path)
    if chromatin_windows is not None:
        spec = _with_chromatin_windows(spec, chromatin_windows)
    chrom_spec = (
        getattr(spec, source, None)
        if source in STREAM_SOURCES
        else None
    )
    if chrom_spec is None:
        raise ValueError(f"{spec.code} has no {source!r} source")
    paths = project_paths(data_dir)
    base = (
        paths.roadmap_covs_csv
        if source == "roadmap"
        else paths.encode_chromatin_covs_csv
    )
    windows = tuple(chrom_spec.windows)
    transforms = check_transforms(chrom_spec.transform)
    cache = chromatin_cache_path(
        base, windows=windows, transform=chrom_spec.transform
    )
    partials = cache.with_name(cache.stem + ".partials")
    partials.mkdir(parents=True, exist_ok=True)
    cache.parent.mkdir(parents=True, exist_ok=True)

    gtfs = (
        dict(gencode_gtfs)
        if gencode_gtfs
        else _default_gencode_gtfs()
    )
    genes = _gene_table(
        gtfs[chrom_spec.assembly], ("protein_coding",), windows
    )
    jobs = (
        list(jobs)
        if jobs is not None
        else track_jobs(spec, source, paths)
    )
    todo = [
        j
        for j in jobs
        if not (partials / f"{j.dest.stem}.csv").exists()
    ]
    logger.info(
        "%s %s: %d of %d tracks to stream (%d done before)",
        spec.code,
        source,
        len(todo),
        len(jobs),
        len(jobs) - len(todo),
    )

    def one(job: TrackJob) -> str | None:
        job.dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            path = download_file(
                job.url,
                job.dest,
                expected_md5=job.md5,
                expected_size=job.size,
                session=requests.Session(),
            )
        except Exception as exc:
            if job.required:
                raise
            logger.warning("Skipping %s: %s", job.dest.name, exc)
            return None
        try:
            frame = _summaries(path, genes, windows, transforms)
            tmp = partials / f"{job.dest.stem}.csv.part"
            frame.to_csv(tmp)
            tmp.rename(partials / f"{job.dest.stem}.csv")
            return job.dest.name
        finally:
            if not keep_files:
                path.unlink(missing_ok=True)

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for k, _ in enumerate(ex.map(one, todo), start=1):
            if k % 25 == 0 or k == len(todo):
                logger.info(
                    "%s %s: %d/%d", spec.code, source, k, len(todo)
                )

    names = sorted(
        j.dest.name
        for j in jobs
        if (partials / f"{j.dest.stem}.csv").exists()
    )
    if not names:
        raise RuntimeError(
            f"No track of {source} could be summarised"
        )
    merged = pd.concat(
        [
            pd.read_csv(partials / f"{Path(n).stem}.csv", index_col=0)
            for n in names
        ],
        axis=1,
    )
    merged.to_csv(cache)
    record = cache.parent / f"{cache.name}.tracks.json"
    record.write_text(json.dumps(names, indent=1) + "\n")
    logger.info(
        "%s %s: %d tracks -> %s %s",
        spec.code,
        source,
        len(names),
        cache,
        merged.shape,
    )
    return cache
