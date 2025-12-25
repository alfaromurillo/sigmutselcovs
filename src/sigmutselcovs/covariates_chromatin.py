"""Chromatin covariate utilities."""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pyBigWig

import pandas as pd

from covariates_utilities import load_gene_bodies_from_gtf
from covariates_utilities import normalize_chromosome_name
from covariates_utilities import sanitize_feature_label


logger = logging.getLogger(__name__)


@dataclass(slots=True)
class TrackSpec:
    """Describe a chromatin feature track."""

    label: str
    path: Path
    statistic: str = "mean"


def _fetch_stat(
        bw: pyBigWig.pyBigWig,
        chrom: str,
        start: int,
        end: int,
        stat: str) -> float:
    """Return a summary statistic for a genomic interval.

    Parameters
    ----------
    bw : pyBigWig.pyBigWig
        Open bigWig handle.
    chrom : str
        Chromosome name recognised by the bigWig header.
    start : int
        Zero-based inclusive start coordinate.
    end : int
        Zero-based exclusive end coordinate.
    stat : str
        Summary statistic accepted by ``pyBigWig.stats``
        (e.g., ``mean``).

    Returns
    -------
    float
        Requested statistic, or ``nan`` when the query cannot be
        satisfied.
    """
    if end <= start:
        return float("nan")
    try:
        val = bw.stats(
            chrom,
            int(start),
            int(end),
            type=stat,
            exact=True)
    except RuntimeError:
        return float("nan")
    if not val:
        return float("nan")
    item = val[0]
    return float("nan") if item is None else float(item)


def summarize_bigwig_over_genes(
        bw_path: str | Path,
        genes: pd.DataFrame,
        *,
        label: str | None = None,
        statistic: str = "mean") -> pd.Series:
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
                "Chromosomes absent in %s: %s",
                path.name,
                missing)

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
        statistic: str = "mean") -> pd.Series:
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
                bw,
                chrom_norm,
                prom_start,
                prom_end,
                statistic)
            values[gene_id] = val

        if missing_chroms:
            missing = ", ".join(sorted(missing_chroms))
            logger.warning(
                "Chromosomes absent in %s: %s",
                path.name,
                missing)

    series = pd.Series(values, name=sanitize_feature_label(label))
    series.index.name = genes.index.name
    return series


def load_tracks(
        tracks: Sequence[str | Path | TrackSpec]
        | Mapping[str, str | Path],
        *,
        default_statistic: str = "mean") -> list[TrackSpec]:
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
                    statistic=default_statistic))
        return specs

    for item in tracks:
        if isinstance(item, TrackSpec):
            specs.append(item)
            continue
        specs.append(
            TrackSpec(
                label=Path(str(item)).stem,
                path=Path(item),
                statistic=default_statistic))

    return specs


def summarise_tracks_to_genes(
        tracks: Sequence[str | Path | TrackSpec]
        | Mapping[str, str | Path],
        genes: pd.DataFrame,
        *,
        include_promoter: bool = True,
        promoter_upstream: int = 2000,
        promoter_downstream: int = 200) -> pd.DataFrame:
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
            statistic=spec.statistic)

        out[body_series.name] = body_series

        if include_promoter:
            prom_series = summarize_promoter_signal(
                spec.path,
                genes,
                upstream=promoter_upstream,
                downstream=promoter_downstream,
                label=f"{spec.label}_promoter",
                statistic=spec.statistic)

            out[prom_series.name] = prom_series

    return out


def load_or_generate_chromatin_covariates(
        location_df: str | Path,
        tracks: Sequence[str | Path | TrackSpec]
        | Mapping[str, str | Path],
        gtf_path: str | Path,
        *,
        biotypes: Iterable[str] | None = ("protein_coding",),
        include_promoter: bool = True,
        promoter_upstream: int = 2000,
        promoter_downstream: int = 200,
        force_generation: bool = False,
        average_by_assay: bool = False) -> pd.DataFrame:
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
            location_df,)
        cached = pd.read_csv(location_df, index_col=0)
        if not cached.index.is_unique:
            cached = cached.groupby(level=0).mean()
        return cached

    logger.info("Preparing gene coordinates from %s", gtf_path)
    genes = load_gene_bodies_from_gtf(
        str(gtf_path),
        biotypes=list(biotypes) if biotypes is not None else None,
        add_chr_prefix_if_needed=True,
        autosomes_only=False)

    if "start" not in genes.columns:
        genes = genes.rename(columns={"region_start": "start",
                                      "region_end": "end"})
    genes = genes.loc[:, ["Chromosome", "start", "end", "strand"]]

    cov_df = summarise_tracks_to_genes(
        tracks,
        genes,
        include_promoter=include_promoter,
        promoter_upstream=promoter_upstream,
        promoter_downstream=promoter_downstream)
    if not cov_df.index.is_unique:
        cov_df = cov_df.groupby(level=0).mean()

    if average_by_assay and not cov_df.empty:
        # Collapse columns by assay (e.g., h3k27ac) while preserving
        # body vs promoter suffixes.
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
        cov_df = pd.DataFrame(collapsed, index=cov_df.index)

    cov_df.to_csv(location_df)
    logger.info("Saved chromatin covariates to %s", location_df)
    return cov_df
