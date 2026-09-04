"""Chromatin covariate utilities."""

import logging
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyBigWig

from .covariates_utilities import (
    fetch_bigwig_stat,
    load_gene_bodies_from_gtf,
    normalize_chromosome_name,
    sanitize_feature_label,
)

logger = logging.getLogger(__name__)


@dataclass(slots=True)
class TrackSpec:
    """Describe a chromatin feature track."""

    label: str
    path: Path
    statistic: str = "mean"


# Shared with the replication-timing bigWig ingestion.
_fetch_stat = fetch_bigwig_stat


def summarize_bigwig_over_genes(
    bw_path: str | Path,
    genes: pd.DataFrame,
    *,
    label: str | None = None,
    statistic: str = "mean",
) -> pd.Series:
    """Return bigWig signal summarised across each gene body.

    Parameters
    ----------
    bw_path : str | pathlib.Path
        Path to the bigWig track on disk.
    genes : pandas.DataFrame
        Gene table indexed by ``ensembl_gene_id`` with columns
        ``Chromosome``, ``start``, ``end`` (0-based, half-open).
    label : str, optional
        Column label to use in the resulting Series. Defaults to the
        file stem.
    statistic : str, default "mean"
        Summary statistic requested from ``pyBigWig.stats``.

    Returns
    -------
    pandas.Series
        Gene-indexed series containing the summarised signal.
    """
    path = Path(bw_path)
    if label is None:
        label = path.stem

    logger.info("Summarising %s over %d genes", path, genes.shape[0])

    with pyBigWig.open(str(path)) as bw:
        chroms = bw.chroms()
        data: dict[str, float] = {}
        missing_chroms: set[str] = set()

        for gene_id, row in genes.iterrows():
            chrom = str(row["Chromosome"])
            start = int(row["start"])
            end = int(row["end"])

            chrom_norm = normalize_chromosome_name(chrom, chroms)
            if chrom_norm is None:
                missing_chroms.add(chrom)
                data[gene_id] = float("nan")
                continue

            val = _fetch_stat(bw, chrom_norm, start, end, statistic)
            data[gene_id] = val

        if missing_chroms:
            missing = ", ".join(sorted(missing_chroms))
            logger.warning(
                "Chromosomes absent in %s: %s", path.name, missing
            )

    series = pd.Series(data, name=sanitize_feature_label(label))
    series.index.name = genes.index.name
    return series


def summarize_promoter_signal(
    bw_path: str | Path,
    genes: pd.DataFrame,
    *,
    upstream: int = 2000,
    downstream: int = 200,
    label: str | None = None,
    statistic: str = "mean",
) -> pd.Series:
    """Summarise bigWig signal across strand-aware promoters.

    Parameters
    ----------
    bw_path : str | pathlib.Path
        Path to the bigWig track on disk.
    genes : pandas.DataFrame
        Gene table with ``Chromosome``, ``start``, ``end``,
        ``strand`` columns.
    upstream : int, default 2000
        Number of bases to extend upstream of the TSS.
    downstream : int, default 200
        Number of bases to extend downstream of the TSS.
    label : str, optional
        Column label to use in the resulting Series. Defaults to the
        file stem plus ``_promoter``.
    statistic : str, default "mean"
        Summary statistic requested from ``pyBigWig.stats``.

    Returns
    -------
    pandas.Series
        Gene-indexed series containing promoter-level summaries.
    """
    path = Path(bw_path)
    if label is None:
        label = f"{path.stem}_promoter"

    logger.info("Summarising promoters for %s", path)

    with pyBigWig.open(str(path)) as bw:
        chroms = bw.chroms()
        values: dict[str, float] = {}
        missing_chroms: set[str] = set()

        for gene_id, row in genes.iterrows():
            chrom = str(row["Chromosome"])
            start = int(row["start"])
            end = int(row["end"])
            strand = str(row.get("strand", "+") or "+")

            if strand == "-":
                prom_start = max(end - downstream, 0)
                prom_end = end + upstream
            else:
                prom_start = max(start - upstream, 0)
                prom_end = start + downstream

            chrom_norm = normalize_chromosome_name(chrom, chroms)
            if chrom_norm is None:
                missing_chroms.add(chrom)
                values[gene_id] = float("nan")
                continue

            val = _fetch_stat(
                bw, chrom_norm, prom_start, prom_end, statistic
            )
            values[gene_id] = val

        if missing_chroms:
            missing = ", ".join(sorted(missing_chroms))
            logger.warning(
                "Chromosomes absent in %s: %s", path.name, missing
            )

    series = pd.Series(values, name=sanitize_feature_label(label))
    series.index.name = genes.index.name
    return series


def load_tracks(
    tracks: (
        Sequence[str | Path | TrackSpec] | Mapping[str, str | Path]
    ),
    *,
    default_statistic: str = "mean",
) -> list[TrackSpec]:
    """Normalise a collection of bigWig track definitions.

    Parameters
    ----------
    tracks : sequence or mapping
        Either an iterable of paths/``TrackSpec`` objects or a mapping
        from label to path.
    default_statistic : str, default "mean"
        Statistic to use for entries that do not specify one.

    Returns
    -------
    list[TrackSpec]
        List of fully populated ``TrackSpec`` instances.
    """
    specs: list[TrackSpec] = []

    if isinstance(tracks, Mapping):
        for label, path in tracks.items():
            specs.append(
                TrackSpec(
                    label=str(label),
                    path=Path(path),
                    statistic=default_statistic,
                )
            )
        return specs

    for item in tracks:
        if isinstance(item, TrackSpec):
            specs.append(item)
            continue
        specs.append(
            TrackSpec(
                label=Path(str(item)).stem,
                path=Path(item),
                statistic=default_statistic,
            )
        )

    return specs


def summarise_tracks_to_genes(
    tracks: (
        Sequence[str | Path | TrackSpec] | Mapping[str, str | Path]
    ),
    genes: pd.DataFrame,
    *,
    include_promoter: bool = True,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
) -> pd.DataFrame:
    """Stack gene-level summaries for each track into a DataFrame.

    Parameters
    ----------
    tracks : sequence or mapping
        Collection of tracks describable by ``TrackSpec``.
    genes : pandas.DataFrame
        Gene table indexed by ``ensembl_gene_id`` with the
        columns used by the summary helpers.
    include_promoter : bool, default True
        If True, add promoter summaries alongside body summaries.
    promoter_upstream : int, default 2000
        Promoter extension upstream of the TSS when relevant.
    promoter_downstream : int, default 200
        Promoter extension downstream of the TSS.

    Returns
    -------
    pandas.DataFrame
        Gene-indexed DataFrame containing one column per requested
        summary.
    """
    specs = load_tracks(tracks)
    out = pd.DataFrame(index=genes.index)

    for spec in specs:
        body_series = summarize_bigwig_over_genes(
            spec.path,
            genes,
            label=f"{spec.label}_body",
            statistic=spec.statistic,
        )

        out[body_series.name] = body_series

        if include_promoter:
            prom_series = summarize_promoter_signal(
                spec.path,
                genes,
                upstream=promoter_upstream,
                downstream=promoter_downstream,
                label=f"{spec.label}_promoter",
                statistic=spec.statistic,
            )

            out[prom_series.name] = prom_series

    return out


def summarize_tracks_to_genes_streaming(
    track_source: Iterable[tuple[str, str | Path]],
    genes: pd.DataFrame,
    *,
    include_promoter: bool = True,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    statistic: str = "mean",
) -> pd.DataFrame:
    """Like :func:`summarise_tracks_to_genes`, but consumes tracks
    lazily from an iterable of ``(label, path)`` pairs instead of a
    pre-materialized list of already-downloaded files.

    The point: pair this with a generator that downloads one track,
    ``yield``s it, and deletes it once this function's loop resumes
    (e.g. :func:`download.stream_roadmap_tracks`/
    ``stream_encode_chromatin_tracks``) -- then disk usage never
    exceeds ~1 bigwig regardless of how many tracks are pooled,
    unlike the batch download-everything-then-summarize-everything
    path :func:`load_or_generate_chromatin_covariates` normally
    takes. Built for large pools (e.g. ``GENERIC``) where downloading
    every track up front doesn't fit on disk.

    Parameters
    ----------
    track_source : iterable of (str, str | pathlib.Path)
        Yields ``(label, bigwig_path)`` one at a time. Each path must
        still exist when this function reads it (the *caller's*
        generator is responsible for deleting it only *after*
        control returns here, not before).
    genes, include_promoter, promoter_upstream, promoter_downstream,
    statistic
        Same as :func:`summarise_tracks_to_genes`.

    Returns
    -------
    pandas.DataFrame
        Same shape/columns as ``summarise_tracks_to_genes`` would
        produce from the same tracks pre-downloaded.
    """
    columns: dict[str, pd.Series] = {}
    for label, path in track_source:
        body_series = summarize_bigwig_over_genes(
            path,
            genes,
            label=f"{label}_body",
            statistic=statistic,
        )
        columns[body_series.name] = body_series

        if include_promoter:
            prom_series = summarize_promoter_signal(
                path,
                genes,
                upstream=promoter_upstream,
                downstream=promoter_downstream,
                label=f"{label}_promoter",
                statistic=statistic,
            )
            columns[prom_series.name] = prom_series

    return pd.DataFrame(columns, index=genes.index)


def _load_gene_bodies(
    gtf_path: str | Path, biotypes: Iterable[str] | None
) -> pd.DataFrame:
    logger.info("Preparing gene coordinates from %s", gtf_path)
    genes = load_gene_bodies_from_gtf(
        str(gtf_path),
        biotypes=list(biotypes) if biotypes is not None else None,
        add_chr_prefix_if_needed=True,
        autosomes_only=False,
    )
    if "start" not in genes.columns:
        genes = genes.rename(
            columns={"region_start": "start", "region_end": "end"}
        )
    return genes.loc[:, ["Chromosome", "start", "end", "strand"]]


def _collapse_by_assay(cov_df: pd.DataFrame) -> pd.DataFrame:
    """Average columns across tracks that share the same assay (e.g.
    all tracks matching ``H3K27ac``), separately for gene body and
    promoter summaries. Column names become, for example,
    ``h3k27ac_body`` and ``h3k27ac_promoter``."""
    if cov_df.empty:
        return cov_df
    import re

    def assay_key(col: str) -> str:
        base = col
        if base.endswith("_body"):
            suffix = "_body"
            base = base[: -len("_body")]
        elif base.endswith("_promoter"):
            suffix = "_promoter"
            base = base[: -len("_promoter")]
        else:
            suffix = ""

        tokens = base.split("_")
        pat = re.compile(r"^h[23]k\d+(?:ac|me\d)$")
        hit = None
        for t in tokens:
            if pat.match(t):
                hit = t
                break
        key = (hit or tokens[-1]) + suffix
        return key

    groups: dict[str, list[str]] = {}
    for c in cov_df.columns:
        k = assay_key(c)
        groups.setdefault(k, []).append(c)

    collapsed = {}
    for k, cols in groups.items():
        collapsed[k] = cov_df[cols].mean(axis=1)
    return pd.DataFrame(collapsed, index=cov_df.index)


def load_or_generate_chromatin_covariates(
    location_df: str | Path,
    tracks: (
        Sequence[str | Path | TrackSpec] | Mapping[str, str | Path]
    ),
    gtf_path: str | Path,
    *,
    biotypes: Iterable[str] | None = ("protein_coding",),
    include_promoter: bool = True,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    force_generation: bool = False,
    average_by_assay: bool = False,
) -> pd.DataFrame:
    """Build or cache gene-level chromatin covariates.

    Parameters
    ----------
    location_df : str | pathlib.Path
        Where to read/write the cached pickled DataFrame.
    tracks : sequence or mapping
        Collection of tracks describable by ``TrackSpec``.
    gtf_path : str | pathlib.Path
        Path to the reference GTF used to obtain gene loci.
    biotypes : iterable of str or None, optional
        Gene biotypes to retain; ``None`` keeps all gene entries.
    include_promoter : bool, default True
        Include promoter summaries alongside gene-body values.
    promoter_upstream : int, default 2000
        Promoter extension upstream of the TSS.
    promoter_downstream : int, default 200
        Promoter extension downstream of the TSS.
    force_generation : bool, default False
        Recompute even if the cache is present.
    average_by_assay : bool, default False
        If True, average columns across tissues that share the same
        assay (e.g., all tracks matching ``H3K27ac``), separately for
        gene body and promoter summaries. Column names become, for
        example, ``h3k27ac_body`` and ``h3k27ac_promoter``.

    Returns
    -------
    pandas.DataFrame
        Gene-indexed chromatin covariate table suitable for
        ``cov_matrix_full``.
    """
    location_df = Path(location_df)

    if location_df.exists() and not force_generation:
        logger.info(
            "Loading chromatin covariates from %s",
            location_df,
        )
        cached = pd.read_csv(location_df, index_col=0)
        if not cached.index.is_unique:
            cached = cached.groupby(level=0).mean()
        return cached

    genes = _load_gene_bodies(gtf_path, biotypes)

    cov_df = summarise_tracks_to_genes(
        tracks,
        genes,
        include_promoter=include_promoter,
        promoter_upstream=promoter_upstream,
        promoter_downstream=promoter_downstream,
    )
    if not cov_df.index.is_unique:
        cov_df = cov_df.groupby(level=0).mean()

    if average_by_assay:
        cov_df = _collapse_by_assay(cov_df)

    cov_df.to_csv(location_df)
    logger.info("Saved chromatin covariates to %s", location_df)
    return cov_df


def load_or_generate_chromatin_covariates_streaming(
    location_df: str | Path,
    track_source: Iterable[tuple[str, str | Path]],
    gtf_path: str | Path,
    *,
    biotypes: Iterable[str] | None = ("protein_coding",),
    include_promoter: bool = True,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    force_generation: bool = False,
    average_by_assay: bool = False,
) -> pd.DataFrame:
    """Streaming counterpart of :func:`load_or_generate_chromatin_covariates`.

    Same cache-check, gene-loading, and (optional) ``average_by_assay``
    collapsing behavior, but builds the raw per-track columns via
    :func:`summarize_tracks_to_genes_streaming` -- meant for a
    ``track_source`` generator that downloads one bigwig at a time and
    deletes it as this reads it (e.g. :func:`download.
    stream_roadmap_tracks`/``stream_encode_chromatin_tracks``), so
    building a covariate matrix from a large pool of tracks (e.g.
    ``GENERIC``) never needs more than ~1 bigwig on disk at once,
    regardless of how many are pooled.

    ``track_source`` is consumed exactly once (it's a generator, not
    a re-iterable collection) -- call this function again with a
    fresh generator if you need to rebuild.

    Parameters, return value: same as
    :func:`load_or_generate_chromatin_covariates`, except ``tracks``
    is replaced by ``track_source``.
    """
    location_df = Path(location_df)

    if location_df.exists() and not force_generation:
        logger.info(
            "Loading chromatin covariates from %s", location_df
        )
        cached = pd.read_csv(location_df, index_col=0)
        if not cached.index.is_unique:
            cached = cached.groupby(level=0).mean()
        return cached

    genes = _load_gene_bodies(gtf_path, biotypes)

    cov_df = summarize_tracks_to_genes_streaming(
        track_source,
        genes,
        include_promoter=include_promoter,
        promoter_upstream=promoter_upstream,
        promoter_downstream=promoter_downstream,
    )
    if not cov_df.index.is_unique:
        cov_df = cov_df.groupby(level=0).mean()

    if average_by_assay:
        cov_df = _collapse_by_assay(cov_df)

    cov_df.to_csv(location_df)
    logger.info("Saved chromatin covariates to %s", location_df)
    return cov_df
