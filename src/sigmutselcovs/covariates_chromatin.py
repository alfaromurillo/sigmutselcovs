"""Chromatin covariate utilities."""

import logging
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

from .covariates_utilities import (
    fetch_bigwig_stat,
    load_gene_bodies_from_gtf,
    normalize_chromosome_name,
    sanitize_feature_label,
)
from .registry import (
    CHROMATIN_TRANSFORMS,
    CHROMATIN_WINDOWS,
    DEFAULT_CHROMATIN_WINDOWS,
    PEAK_TYPES,
    TRANSFORM_SUFFIXES,
    chromatin_transforms,
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


# Windows whose start/end come from the CDS span columns that
# `load_gene_bodies_from_gtf(..., with_cds_span=True)` adds.
_CDS_WINDOWS = ("cds_span", "cds_span2k")
_CDS_FLANK = 2000
_LOG_SUFFIX = TRANSFORM_SUFFIXES["log1p"]
# Bases of one bigWig read in `_log1p_interval_means`: bounds memory
# (a float64 log and a cumulative sum per base) at ~200 MB however
# long the chromosome, while a gene of several Mb still fits in one.
_LOG1P_CHUNK_BASES = 1 << 23


def check_windows(windows: Iterable[str]) -> tuple[str, ...]:
    """Validate a sequence of window names.

    Parameters
    ----------
    windows : iterable of str
        Window names, each one of ``CHROMATIN_WINDOWS``.

    Returns
    -------
    tuple[str, ...]
        The windows, in the order given.

    Raises
    ------
    TypeError
        If ``windows`` is a bare string.
    ValueError
        If ``windows`` is empty, repeats a window or names an unknown
        one.
    """
    if isinstance(windows, str):
        raise TypeError(
            f"windows must be a sequence of names, not {windows!r}"
        )
    windows = tuple(windows)
    if not windows:
        raise ValueError("windows must name at least one window")
    unknown = [w for w in windows if w not in CHROMATIN_WINDOWS]
    if unknown:
        raise ValueError(
            f"Unknown chromatin windows {unknown}; valid: "
            f"{CHROMATIN_WINDOWS}"
        )
    if len(set(windows)) != len(windows):
        raise ValueError(f"Repeated chromatin windows in {windows}")
    return windows


def check_transform(transform: str | None) -> str | None:
    """Validate a chromatin signal transform.

    ``None`` (the linear mean), ``"log1p"`` (the mean of
    ``log1p(signal)`` per base) or ``"log1p_mean"`` (``log1p`` of the
    linear mean, the mean clipped at 0 first).
    """
    if (
        transform is not None
        and transform not in CHROMATIN_TRANSFORMS
    ):
        raise ValueError(
            f"Unknown chromatin transform {transform!r}; valid: "
            f"None, {CHROMATIN_TRANSFORMS}"
        )
    return transform


def check_transforms(transform) -> tuple[str | None, ...]:
    """Validate a scalar transform or a sequence of them.

    Returns
    -------
    tuple
        The transforms, as `registry.chromatin_transforms` gives them.

    Raises
    ------
    ValueError
        If a transform is unknown, the sequence is empty or it
        repeats one.
    """
    transforms = chromatin_transforms(transform)
    if not transforms:
        raise ValueError("transform must name at least one transform")
    for t in transforms:
        check_transform(t)
    if len(set(transforms)) != len(transforms):
        raise ValueError(
            f"Repeated chromatin transforms in {transforms}"
        )
    return transforms


def gene_windows(
    genes: pd.DataFrame,
    window: str,
    *,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
) -> pd.DataFrame:
    """Per-gene coordinates of one summary window.

    Parameters
    ----------
    genes : pandas.DataFrame
        Gene table with ``Chromosome``, ``start``, ``end`` (0-based,
        half-open) and ``strand``; the CDS windows also need
        ``cds_start``/``cds_end`` (see `load_gene_bodies_from_gtf`'s
        ``with_cds_span``).
    window : str
        One of ``CHROMATIN_WINDOWS``:

        - ``"body"``: the gene, ``[start, end)``.
        - ``"promoter"``: strand-aware ``[TSS - promoter_upstream,
          TSS + promoter_downstream)``, as
          `summarize_promoter_signal` computes it.
        - ``"cds_span"``: the first to the last coding base,
          ``[cds_start, cds_end)``; NaN for a gene with no CDS.
        - ``"cds_span2k"``: the CDS span extended by 2 kb on both
          sides.
    promoter_upstream, promoter_downstream : int
        Promoter extent around the TSS, in bases.

    Returns
    -------
    pandas.DataFrame
        Indexed like ``genes``, with ``Chromosome`` and float
        ``start``/``end`` (NaN where the window is undefined). Starts
        are clipped at 0; ends are not clipped to the chromosome
        length, which the gene table does not know.
    """
    check_windows([window])
    start = genes["start"].to_numpy(dtype=float)
    end = genes["end"].to_numpy(dtype=float)
    if window == "body":
        w_start, w_end = start, end
    elif window == "promoter":
        if "strand" in genes.columns:
            strand = genes["strand"].fillna("+").astype(str)
            minus = (strand == "-").to_numpy()
        else:
            minus = np.zeros(len(genes), dtype=bool)
        w_start = np.where(
            minus,
            np.maximum(end - promoter_downstream, 0),
            np.maximum(start - promoter_upstream, 0),
        )
        w_end = np.where(
            minus,
            end + promoter_upstream,
            start + promoter_downstream,
        )
    else:
        missing = {"cds_start", "cds_end"} - set(genes.columns)
        if missing:
            raise ValueError(
                f"Window {window!r} needs the gene table's cds_start/"
                "cds_end columns: load it with "
                "load_gene_bodies_from_gtf(..., with_cds_span=True)"
            )
        flank = _CDS_FLANK if window == "cds_span2k" else 0
        cds_start = genes["cds_start"].to_numpy(
            dtype=float, na_value=np.nan
        )
        cds_end = genes["cds_end"].to_numpy(
            dtype=float, na_value=np.nan
        )
        w_start = np.maximum(cds_start - flank, 0)
        w_end = cds_end + flank
    return pd.DataFrame(
        {
            "Chromosome": genes["Chromosome"].astype(str).to_numpy(),
            "start": w_start,
            "end": w_end,
        },
        index=genes.index,
    )


def window_column(
    label: str, window: str, transform: str | None = None
) -> str:
    """Column name of a bigWig summary: ``<label>_<window>`` plus the
    transform's suffix (``_log`` for ``"log1p"``, ``_logmean`` for
    ``"log1p_mean"``; see `TRANSFORM_SUFFIXES`)."""
    name = sanitize_feature_label(f"{label}_{window}")
    return name + TRANSFORM_SUFFIXES[transform]


def split_transform_suffix(column: str) -> tuple[str, str | None]:
    """A summary column's name without its transform suffix, and the
    transform (the inverse of `window_column`'s suffix)."""
    for transform in ("log1p_mean", "log1p"):
        suffix = TRANSFORM_SUFFIXES[transform]
        if column.endswith(suffix) and column[
            : -len(suffix)
        ].endswith(_WINDOW_SUFFIXES):
            return column[: -len(suffix)], transform
    return column, None


def _gene_series(
    genes: pd.DataFrame, values: Sequence[float], name: str
) -> pd.Series:
    """Gene-indexed series, keyed like `summarize_bigwig_over_genes`
    (a gene id listed twice keeps its last row's value)."""
    series = pd.Series(dict(zip(genes.index, values)), name=name)
    series.index.name = genes.index.name
    return series


def _resolve_chroms(
    chrom_names: Iterable[str], available: Mapping, path: Path
) -> list[str | None]:
    """Each name as the track spells it (None if absent), warning once
    about the absent ones."""
    resolved = []
    missing: set[str] = set()
    for chrom in chrom_names:
        norm = normalize_chromosome_name(chrom, available)
        if norm is None:
            missing.add(chrom)
        resolved.append(norm)
    if missing:
        logger.warning(
            "Chromosomes absent in %s: %s",
            path.name,
            ", ".join(sorted(missing)),
        )
    return resolved


def _linear_interval_means(
    bw, chroms: Sequence[str | None], starts, ends, statistic: str
) -> np.ndarray:
    """``pyBigWig.stats(..., exact=True)`` per interval, ends clipped
    to the chromosome length."""
    lengths = bw.chroms()
    out = np.full(len(chroms), np.nan)
    for i, (chrom, s, e) in enumerate(zip(chroms, starts, ends)):
        if chrom is None or not (np.isfinite(s) and np.isfinite(e)):
            continue
        e = min(e, lengths[chrom])
        out[i] = _fetch_stat(bw, chrom, int(s), int(e), statistic)
    return out


def _log1p_interval_means(
    bw, chroms: Sequence[str | None], starts, ends
) -> np.ndarray:
    """Mean of ``log1p(signal)`` per base over each interval.

    Averages over the bases the track covers (the bases
    ``bw.values`` does not return as NaN, which are the ones
    ``pyBigWig.stats(..., exact=True)`` counts), so a constant signal
    ``c`` gives ``log1p(c)`` whatever the gaps. Intervals are clipped
    to ``[0, chromosome length)``; one that covers no base is NaN.

    Each chromosome is read in chunks of at most
    ``_LOG1P_CHUNK_BASES`` holding whole intervals (an interval
    longer than that is a chunk of its own); within a chunk, interval
    sums come from cumulative sums of the per-base log and of the
    covered-base indicator.
    """
    lengths = bw.chroms()
    chroms = np.asarray(chroms, dtype=object)
    starts = np.asarray(starts, dtype=float)
    ends = np.asarray(ends, dtype=float)
    out = np.full(len(chroms), np.nan)
    for chrom in {c for c in chroms if c is not None}:
        length = lengths[chrom]
        idx = np.flatnonzero(
            (chroms == chrom)
            & np.isfinite(starts)
            & np.isfinite(ends)
        )
        s = np.clip(starts[idx], 0, length).astype(np.int64)
        e = np.clip(ends[idx], 0, length).astype(np.int64)
        keep = e > s
        idx, s, e = idx[keep], s[keep], e[keep]
        order = np.argsort(s, kind="stable")
        idx, s, e = idx[order], s[order], e[order]

        n = len(idx)
        first = 0
        while first < n:
            chunk_start, chunk_end = s[first], e[first]
            last = first + 1
            while (
                last < n
                and max(chunk_end, e[last]) - chunk_start
                <= _LOG1P_CHUNK_BASES
            ):
                chunk_end = max(chunk_end, e[last])
                last += 1
            values = bw.values(
                chrom, int(chunk_start), int(chunk_end), numpy=True
            )
            covered = ~np.isnan(values)
            logs = np.zeros(len(values), dtype=np.float64)
            logs[covered] = np.log1p(
                values[covered].astype(np.float64)
            )
            log_sum = np.concatenate(([0.0], np.cumsum(logs)))
            n_covered = np.concatenate(
                ([0], np.cumsum(covered, dtype=np.int64))
            )
            lo = s[first:last] - chunk_start
            hi = e[first:last] - chunk_start
            count = n_covered[hi] - n_covered[lo]
            with np.errstate(invalid="ignore", divide="ignore"):
                means = (log_sum[hi] - log_sum[lo]) / count
            out[idx[first:last]] = np.where(count > 0, means, np.nan)
            first = last
    return out


def summarize_bigwig_windows(
    bw_path: str | Path,
    genes: pd.DataFrame,
    *,
    windows: Sequence[str] = DEFAULT_CHROMATIN_WINDOWS,
    transform: str | None = None,
    label: str | None = None,
    statistic: str = "mean",
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
) -> pd.DataFrame:
    """Summarise one bigWig over several per-gene windows.

    Parameters
    ----------
    bw_path : str | pathlib.Path
        Path to the bigWig track on disk.
    genes : pandas.DataFrame
        Gene table as `gene_windows` takes it.
    windows : sequence of str, default ("body", "promoter")
        Windows to summarise over, each one of ``CHROMATIN_WINDOWS``.
    transform : {None, "log1p", "log1p_mean"}, default None
        ``None`` gives ``statistic`` of the signal. ``"log1p"`` gives
        the mean over the window of ``log1p(signal)`` per base (not
        the log of the mean), over the bases the track covers; it
        assumes a signal above -1 (fold change, p-values, coverage).
        ``"log1p_mean"`` gives ``log1p`` of the window's linear mean,
        the mean clipped at 0 first (a fold change is not negative).
        Both need ``statistic="mean"``.
    label : str, optional
        Column label stem; defaults to the file stem.
    statistic : str, default "mean"
        Summary statistic requested from ``pyBigWig.stats`` when
        ``transform`` is None.
    promoter_upstream, promoter_downstream : int
        Promoter extent around the TSS, in bases.

    Returns
    -------
    pandas.DataFrame
        Indexed like ``genes``, one column per window, named by
        `window_column` (``<label>_<window>``, plus the transform's
        suffix), in the order of ``windows``.

    Notes
    -----
    Without a transform, the ``body`` and ``promoter`` columns are
    exactly `summarize_bigwig_over_genes` and
    `summarize_promoter_signal`, which give NaN for a window running
    past the chromosome end. The CDS windows and the log transforms
    clip such a window to the chromosome instead.
    """
    windows = check_windows(windows)
    check_transform(transform)
    if transform is not None and statistic != "mean":
        raise ValueError(
            f"transform {transform!r} needs statistic='mean', "
            f"not {statistic!r}"
        )
    path = Path(bw_path)
    if label is None:
        label = path.stem

    columns: dict[str, pd.Series] = {}
    computed = []
    for window in windows:
        name = window_column(label, window, transform)
        if transform is None and window == "body":
            columns[name] = summarize_bigwig_over_genes(
                path, genes, label=name, statistic=statistic
            )
        elif transform is None and window == "promoter":
            columns[name] = summarize_promoter_signal(
                path,
                genes,
                upstream=promoter_upstream,
                downstream=promoter_downstream,
                label=name,
                statistic=statistic,
            )
        else:
            columns[name] = None
            computed.append((window, name))

    if computed:
        logger.info(
            "Summarising %s over %s for %d genes",
            path,
            ", ".join(w for w, _ in computed),
            genes.shape[0],
        )
        coords = [
            gene_windows(
                genes,
                window,
                promoter_upstream=promoter_upstream,
                promoter_downstream=promoter_downstream,
            )
            for window, _ in computed
        ]
        with pyBigWig.open(str(path)) as bw:
            chroms = _resolve_chroms(
                genes["Chromosome"].astype(str), bw.chroms(), path
            )
            all_chroms = chroms * len(coords)
            starts = np.concatenate([c["start"] for c in coords])
            ends = np.concatenate([c["end"] for c in coords])
            if transform == "log1p":
                values = _log1p_interval_means(
                    bw, all_chroms, starts, ends
                )
            else:
                values = _linear_interval_means(
                    bw, all_chroms, starts, ends, statistic
                )
                if transform == "log1p_mean":
                    values = np.log1p(np.clip(values, 0, None))
        for k, (_, name) in enumerate(computed):
            block = values[k * len(genes) : (k + 1) * len(genes)]
            columns[name] = _gene_series(genes, block, name)

    out = pd.DataFrame(index=genes.index)
    for name, series in columns.items():
        out[name] = series
    return out


def peak_file_parts(path: str | Path) -> tuple[str, str]:
    """Label and peak type of a peak file, from its name.

    ``E075-H3K4me3.narrowPeak.gz`` gives ``("E075-H3K4me3",
    "narrowPeak")``; the type suffix is matched case-insensitively
    and an optional ``.gz`` is ignored.

    Raises
    ------
    ValueError
        If the name does not end in one of ``PEAK_TYPES``.
    """
    name = Path(path).name
    if name.lower().endswith(".gz"):
        name = name[: -len(".gz")]
    for peak_type in PEAK_TYPES:
        suffix = "." + peak_type
        if name.lower().endswith(suffix.lower()):
            return name[: -len(suffix)], peak_type
    raise ValueError(
        f"{Path(path).name}: not a peak file (expected a name ending "
        f"in one of {[f'.{t}[.gz]' for t in PEAK_TYPES]})"
    )


def read_peak_intervals(
    peak_path: str | Path,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Read a BED-like peak file into merged intervals per chromosome.

    Parameters
    ----------
    peak_path : str | pathlib.Path
        narrowPeak, broadPeak or any BED file, gzipped or not
        (detected from the content, not the name). Only the first
        three columns are read; ``#`` comments and ``track``/
        ``browser`` lines are skipped.

    Returns
    -------
    dict[str, tuple[numpy.ndarray, numpy.ndarray]]
        Chromosome -> ``(starts, ends)``: the union of the peaks as
        sorted, disjoint 0-based half-open intervals (overlapping or
        abutting peaks merged).
    """
    path = Path(peak_path)
    with open(path, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype=str,
        comment="#",
        compression="gzip" if gzipped else None,
    )
    table = table[
        ~table["chrom"].str.startswith(("track", "browser"))
    ]
    starts_all = table["start"].astype(np.int64).to_numpy()
    ends_all = table["end"].astype(np.int64).to_numpy()
    chrom_all = table["chrom"].to_numpy()

    merged: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom in pd.unique(chrom_all):
        on = chrom_all == chrom
        s, e = starts_all[on], ends_all[on]
        order = np.lexsort((e, s))
        s, e = s[order], e[order]
        reach = np.maximum.accumulate(e)
        new = np.concatenate(([True], s[1:] > reach[:-1]))
        heads = np.flatnonzero(new)
        merged[str(chrom)] = (s[heads], np.maximum.reduceat(e, heads))
    return merged


def _covered_before(
    starts: np.ndarray, ends: np.ndarray, cum: np.ndarray, x
) -> np.ndarray:
    """Bases of the disjoint intervals that lie in ``[0, x)``."""
    k = np.searchsorted(starts, x, side="left")
    excess = np.where(
        k > 0, np.maximum(ends[np.maximum(k - 1, 0)] - x, 0), 0
    )
    return cum[k] - excess


def summarize_peak_coverage(
    peak_path: str | Path,
    genes: pd.DataFrame,
    *,
    windows: Sequence[str] = DEFAULT_CHROMATIN_WINDOWS,
    label: str | None = None,
    peak_type: str | None = None,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    chrom_sizes: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """Fraction of each gene window covered by a track's peaks.

    Parameters
    ----------
    peak_path : str | pathlib.Path
        Peak calls (narrowPeak/broadPeak BED, gzipped or not).
    genes : pandas.DataFrame
        Gene table as `gene_windows` takes it.
    windows : sequence of str, default ("body", "promoter")
        Windows to summarise over, each one of ``CHROMATIN_WINDOWS``.
    label : str, optional
        Column label stem; defaults to the file name without its
        peak-type suffix (see `peak_file_parts`).
    peak_type : str, optional
        Peak type for the column names; defaults to the one in the
        file name.
    promoter_upstream, promoter_downstream : int
        Promoter extent around the TSS, in bases.
    chrom_sizes : mapping of str to int, optional
        Chromosome lengths (e.g. `load_chrom_sizes("hg19")`). When
        given, windows are clipped to ``[0, length)`` of their
        chromosome, as the CDS windows and the log transform of
        `summarize_bigwig_windows` are; a chromosome missing from it
        is left unclipped. ``None`` clips at 0 only.

    Returns
    -------
    pandas.DataFrame
        Indexed like ``genes``, one column per window named
        ``<label>_<peaktype>_<window>`` (sanitised, so e.g.
        ``e075_h3k4me3_narrowpeak_body``): the number of the window's
        bases inside the union of the peaks, divided by the window's
        length. NaN where the window is undefined (no CDS) or empty.

    Notes
    -----
    A BED file does not list its chromosomes, so a chromosome with no
    peak is taken to have none: its genes get 0 (a warning names such
    chromosomes). Windows are clipped at the chromosome end only when
    ``chrom_sizes`` is given.
    """
    windows = check_windows(windows)
    path = Path(peak_path)
    if label is None or peak_type is None:
        file_label, file_type = peak_file_parts(path)
        label = file_label if label is None else label
        peak_type = file_type if peak_type is None else peak_type

    logger.info(
        "Summarising peaks %s over %d genes", path, len(genes)
    )
    intervals = read_peak_intervals(path)
    chroms = genes["Chromosome"].astype(str)
    resolved = [
        normalize_chromosome_name(c, intervals) for c in chroms
    ]
    absent = {c for c, r in zip(chroms, resolved) if r is None}
    if absent:
        logger.warning(
            "No peaks in %s on %s; coverage 0 there",
            path.name,
            ", ".join(sorted(absent)),
        )
    resolved = np.asarray(resolved, dtype=object)
    limit = np.full(len(genes), np.inf)
    if chrom_sizes is not None:
        unsized = set()
        for i, chrom in enumerate(chroms):
            norm = normalize_chromosome_name(chrom, chrom_sizes)
            if norm is None:
                unsized.add(chrom)
            else:
                limit[i] = chrom_sizes[norm]
        if unsized:
            logger.warning(
                "No chromosome size for %s; peak windows there are "
                "not clipped at the chromosome end",
                ", ".join(sorted(unsized)),
            )

    out = pd.DataFrame(index=genes.index)
    for window in windows:
        coords = gene_windows(
            genes,
            window,
            promoter_upstream=promoter_upstream,
            promoter_downstream=promoter_downstream,
        )
        w_start = coords["start"].to_numpy()
        w_end = coords["end"].to_numpy()
        if chrom_sizes is not None:
            w_start = np.minimum(w_start, limit)
            w_end = np.minimum(w_end, limit)
        covered = np.zeros(len(genes))
        for chrom, (p_start, p_end) in intervals.items():
            on = np.flatnonzero(
                (resolved == chrom)
                & np.isfinite(w_start)
                & np.isfinite(w_end)
            )
            if not len(on):
                continue
            cum = np.concatenate(([0], np.cumsum(p_end - p_start)))
            covered[on] = _covered_before(
                p_start, p_end, cum, w_end[on]
            ) - _covered_before(p_start, p_end, cum, w_start[on])
        width = w_end - w_start
        with np.errstate(invalid="ignore", divide="ignore"):
            fraction = np.where(width > 0, covered / width, np.nan)
        name = sanitize_feature_label(f"{label}_{peak_type}_{window}")
        out[name] = _gene_series(genes, fraction, name)
    return out


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
    windows: Sequence[str] | None = None,
    transform: str | None = None,
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
    windows : sequence of str, optional
        Windows to summarise each track over (see
        `summarize_bigwig_windows`); when given, ``include_promoter``
        is ignored. ``None`` keeps the body (and promoter) summaries.
    transform : {None, "log1p", "log1p_mean"} or a sequence, optional
        Transform of the signal summaries (see
        `summarize_bigwig_windows`). A sequence such as
        ``[None, "log1p"]`` gives each track's columns once per
        transform, in that order.

    Returns
    -------
    pandas.DataFrame
        Gene-indexed DataFrame containing one column per requested
        summary.
    """
    specs = load_tracks(tracks)
    out = pd.DataFrame(index=genes.index)

    if windows is not None or transform is not None:
        if windows is None:
            windows = (
                ("body", "promoter")
                if include_promoter
                else ("body",)
            )
        transforms = check_transforms(transform)
        for spec in specs:
            for each in transforms:
                frame = summarize_bigwig_windows(
                    spec.path,
                    genes,
                    windows=windows,
                    transform=each,
                    label=spec.label,
                    statistic=spec.statistic,
                    promoter_upstream=promoter_upstream,
                    promoter_downstream=promoter_downstream,
                )
                for name in frame.columns:
                    out[name] = frame[name].to_numpy()
        return out

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
    windows: Sequence[str] | None = None,
    transform: str | None = None,
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
    statistic, windows, transform
        Same as :func:`summarise_tracks_to_genes`.

    Returns
    -------
    pandas.DataFrame
        Same shape/columns as ``summarise_tracks_to_genes`` would
        produce from the same tracks pre-downloaded.
    """
    columns: dict[str, pd.Series] = {}
    if windows is not None or transform is not None:
        if windows is None:
            windows = (
                ("body", "promoter")
                if include_promoter
                else ("body",)
            )
        transforms = check_transforms(transform)
        for label, path in track_source:
            for each in transforms:
                frame = summarize_bigwig_windows(
                    path,
                    genes,
                    windows=windows,
                    transform=each,
                    label=label,
                    statistic=statistic,
                    promoter_upstream=promoter_upstream,
                    promoter_downstream=promoter_downstream,
                )
                for name in frame.columns:
                    columns[name] = frame[name].to_numpy()
        return pd.DataFrame(columns, index=genes.index)

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
    gtf_path: str | Path,
    biotypes: Iterable[str] | None,
    with_cds_span: bool = False,
) -> pd.DataFrame:
    logger.info("Preparing gene coordinates from %s", gtf_path)
    genes = load_gene_bodies_from_gtf(
        str(gtf_path),
        biotypes=list(biotypes) if biotypes is not None else None,
        add_chr_prefix_if_needed=True,
        autosomes_only=False,
        with_cds_span=with_cds_span,
    )
    if "start" not in genes.columns:
        genes = genes.rename(
            columns={"region_start": "start", "region_end": "end"}
        )
    columns = ["Chromosome", "start", "end", "strand"]
    if with_cds_span:
        columns += ["cds_start", "cds_end"]
    return genes.loc[:, columns]


def _gene_table(
    gtf_path: str | Path,
    biotypes: Iterable[str] | None,
    windows: Iterable[str] | None,
) -> pd.DataFrame:
    """The gene table, with the CDS span columns only when a window
    needs them."""
    if windows is not None and any(
        w in _CDS_WINDOWS for w in windows
    ):
        return _load_gene_bodies(
            gtf_path, biotypes, with_cds_span=True
        )
    return _load_gene_bodies(gtf_path, biotypes)


def chromatin_cache_path(
    base_csv: str | Path,
    *,
    windows: Sequence[str] = DEFAULT_CHROMATIN_WINDOWS,
    transform: str | None = None,
    peak_types: Sequence[str] | None = None,
) -> Path:
    """Cache file for one set of chromatin summary options.

    The signal summaries with the default options (body and promoter
    means, no transform) keep ``base_csv`` itself, so existing caches
    stay valid. Any other option set gets its own file next to it,
    named after the options, so a build never reads a cache computed
    under different ones; peak coverage always has its own file.

    Parameters
    ----------
    base_csv : str | pathlib.Path
        The source's default cache, e.g.
        ``chromatin/roadmap/chromatin_covs_roadmap.csv``.
    windows : sequence of str, default ("body", "promoter")
        Summary windows, in column order.
    transform : {None, "log1p", "log1p_mean"} or a sequence, optional
        Signal transform(s). A one-element sequence names the same
        file as its scalar; a longer one is spelled out in order,
        ``None`` as ``linear`` (e.g. ``..._linear-log1p.csv``).
    peak_types : sequence of str, optional
        For the peak-coverage cache: the peak types summarised.

    Returns
    -------
    pathlib.Path
        E.g. ``chromatin_covs_roadmap_cds_span-cds_span2k_log1p.csv``
        or ``chromatin_covs_roadmap_peaks_narrowpeak_body-promoter.csv``.
    """
    base = Path(base_csv)
    windows = check_windows(windows)
    transforms = check_transforms(transform)
    window_tag = "-".join(windows)
    if peak_types is not None:
        types = "-".join(t.lower() for t in peak_types)
        tag = f"peaks_{types}_{window_tag}"
    elif windows == DEFAULT_CHROMATIN_WINDOWS and transforms == (
        None,
    ):
        return base
    elif len(transforms) == 1:
        tag = window_tag + (
            f"_{transforms[0]}" if transforms[0] else ""
        )
    else:
        tag = (
            window_tag
            + "_"
            + "-".join(t or "linear" for t in transforms)
        )
    return base.with_name(f"{base.stem}_{tag}{base.suffix}")


# Longest first, so "_cds_span2k" is not read as "_cds_span".
_WINDOW_SUFFIXES = tuple(
    sorted(
        (f"_{w}" for w in CHROMATIN_WINDOWS), key=len, reverse=True
    )
)
_PEAK_SUFFIXES = tuple(f"_{t.lower()}" for t in PEAK_TYPES)


def _collapse_by_assay(cov_df: pd.DataFrame) -> pd.DataFrame:
    """Average columns across tracks that share the same assay (e.g.
    all tracks matching ``H3K27ac``), separately for each summary:
    window, transform and, for peak coverage, peak type. Column names
    become, for example, ``h3k27ac_body``, ``h3k27ac_promoter``,
    ``h3k27ac_cds_span_log`` and ``h3k27ac_narrowpeak_body``."""
    if cov_df.empty:
        return cov_df
    import re

    def assay_key(col: str) -> str:
        base = col
        tail = ""
        base, transform = split_transform_suffix(base)
        tail = TRANSFORM_SUFFIXES[transform]
        suffix = ""
        for window_suffix in _WINDOW_SUFFIXES:
            if base.endswith(window_suffix):
                suffix = window_suffix
                base = base[: -len(window_suffix)]
                break
        for peak_suffix in _PEAK_SUFFIXES:
            if base.endswith(peak_suffix):
                suffix = peak_suffix + suffix
                base = base[: -len(peak_suffix)]
                break
        suffix += tail

        tokens = base.split("_")
        # h2ak9ac (H2AK9ac) doesn't fit h[23]k... -- the histone
        # subunit token is "2a", not a bare digit, unlike h3k*/h4k*
        # marks. Found live when GENERIC picked up H2AK9ac tracks
        # for the first time (2026-09-07): without this alternation
        # it silently fell through to keying by accession instead of
        # collapsing, the same failure mode "dnase" had before.
        pat = re.compile(r"^h(?:2a|2b|[34])k\d+(?:ac|me\d)$|^dnase$")
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
    windows: Sequence[str] | None = None,
    transform: str | None = None,
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
    windows : sequence of str, optional
        Windows to summarise each track over (see
        `summarize_bigwig_windows`); overrides ``include_promoter``.
        The CDS windows read the GTF's ``CDS`` features as well.
    transform : {None, "log1p", "log1p_mean"} or a sequence, optional
        Transform(s) of the signal summaries, as in
        `summarise_tracks_to_genes`.

    Returns
    -------
    pandas.DataFrame
        Gene-indexed chromatin covariate table suitable for
        ``cov_matrix_full``.

    Notes
    -----
    The cache at ``location_df`` is read whatever options it was built
    with; give each option set its own file (see
    `chromatin_cache_path`).
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

    genes = _gene_table(gtf_path, biotypes, windows)

    cov_df = summarise_tracks_to_genes(
        tracks,
        genes,
        include_promoter=include_promoter,
        promoter_upstream=promoter_upstream,
        promoter_downstream=promoter_downstream,
        windows=windows,
        transform=transform,
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
    windows: Sequence[str] | None = None,
    transform: str | None = None,
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

    genes = _gene_table(gtf_path, biotypes, windows)

    cov_df = summarize_tracks_to_genes_streaming(
        track_source,
        genes,
        include_promoter=include_promoter,
        promoter_upstream=promoter_upstream,
        promoter_downstream=promoter_downstream,
        windows=windows,
        transform=transform,
    )
    if not cov_df.index.is_unique:
        cov_df = cov_df.groupby(level=0).mean()

    if average_by_assay:
        cov_df = _collapse_by_assay(cov_df)

    cov_df.to_csv(location_df)
    logger.info("Saved chromatin covariates to %s", location_df)
    return cov_df


def summarise_peaks_to_genes(
    peak_files: Sequence[str | Path],
    genes: pd.DataFrame,
    *,
    windows: Sequence[str] = DEFAULT_CHROMATIN_WINDOWS,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    chrom_sizes: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """Stack `summarize_peak_coverage` columns for each peak file.

    Parameters
    ----------
    peak_files : sequence of str | pathlib.Path
        Peak files, named ``<label>.<peakType>[.gz]``.
    genes : pandas.DataFrame
        Gene table as `gene_windows` takes it.
    windows : sequence of str, default ("body", "promoter")
        Windows to summarise over.
    promoter_upstream, promoter_downstream : int
        Promoter extent around the TSS, in bases.
    chrom_sizes : mapping of str to int, optional
        Chromosome lengths to clip windows to (see
        `summarize_peak_coverage`).

    Returns
    -------
    pandas.DataFrame
        Gene-indexed, ``<label>_<peaktype>_<window>`` columns, file by
        file in the order given.
    """
    out = pd.DataFrame(index=genes.index)
    for path in peak_files:
        frame = summarize_peak_coverage(
            path,
            genes,
            windows=windows,
            promoter_upstream=promoter_upstream,
            promoter_downstream=promoter_downstream,
            chrom_sizes=chrom_sizes,
        )
        for name in frame.columns:
            out[name] = frame[name].to_numpy()
    return out


def load_or_generate_peak_covariates(
    location_df: str | Path,
    peak_files: Sequence[str | Path],
    gtf_path: str | Path,
    *,
    windows: Sequence[str] = DEFAULT_CHROMATIN_WINDOWS,
    biotypes: Iterable[str] | None = ("protein_coding",),
    promoter_upstream: int = 2000,
    promoter_downstream: int = 200,
    force_generation: bool = False,
    average_by_assay: bool = False,
    chrom_sizes: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """Build or cache gene-level peak-coverage covariates.

    The peak counterpart of `load_or_generate_chromatin_covariates`:
    same cache-first pattern, gene table and ``average_by_assay``
    collapsing, with `summarise_peaks_to_genes` in place of the
    bigWig summaries.

    Parameters
    ----------
    location_df : str | pathlib.Path
        Where to read/write the cached CSV.
    peak_files : sequence of str | pathlib.Path
        Peak files (narrowPeak/broadPeak, gzipped or not).
    gtf_path : str | pathlib.Path
        Reference GTF used to obtain gene loci.
    windows : sequence of str, default ("body", "promoter")
        Windows to summarise over.
    biotypes : iterable of str or None, optional
        Gene biotypes to retain; ``None`` keeps all gene entries.
    promoter_upstream, promoter_downstream : int
        Promoter extent around the TSS, in bases.
    force_generation : bool, default False
        Recompute even if the cache is present.
    average_by_assay : bool, default False
        Average columns across tracks of the same mark, separately
        for each peak type and window (e.g. ``h3k27ac_narrowpeak_body``).
    chrom_sizes : mapping of str to int, optional
        Chromosome lengths to clip windows to (see
        `summarize_peak_coverage`); the builder passes the packaged
        table of the source's assembly.

    Returns
    -------
    pandas.DataFrame
        Gene-indexed peak-coverage fractions.
    """
    location_df = Path(location_df)
    windows = check_windows(windows)

    if location_df.exists() and not force_generation:
        logger.info("Loading peak covariates from %s", location_df)
        cached = pd.read_csv(location_df, index_col=0)
        if not cached.index.is_unique:
            cached = cached.groupby(level=0).mean()
        return cached

    genes = _gene_table(gtf_path, biotypes, windows)
    cov_df = summarise_peaks_to_genes(
        peak_files,
        genes,
        windows=windows,
        promoter_upstream=promoter_upstream,
        promoter_downstream=promoter_downstream,
        chrom_sizes=chrom_sizes,
    )
    if not cov_df.index.is_unique:
        cov_df = cov_df.groupby(level=0).mean()

    if average_by_assay:
        cov_df = _collapse_by_assay(cov_df)

    cov_df.to_csv(location_df)
    logger.info("Saved peak covariates to %s", location_df)
    return cov_df
