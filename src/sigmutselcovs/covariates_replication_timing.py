"""Replication timing covariate processing.

This module provides functions to load and process replication timing
data, including mean replication time per gene and various
replication-related features like initiation zones and termination
sites. These features are used as covariates in mutation rate
modeling.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from collections.abc import Sequence

from .covariates_utilities import read_bed_file

from .covariates_utilities import annotate_indicator_in_region
from .covariates_utilities import annotate_with_binned_features
from .covariates_utilities import fetch_bigwig_stat

from .covariates_utilities import load_gene_bodies_from_gtf
from .covariates_utilities import normalize_chromosome_name

import functools
import logging
import warnings


logger = logging.getLogger(__name__)


def _legacy_kwarg(old: str, new: str):
    """Map a renamed keyword argument with a DeprecationWarning."""

    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if old in kwargs:
                if new in kwargs:
                    raise TypeError(
                        f"{func.__name__}() got both {old!r} and {new!r}"
                    )
                warnings.warn(
                    f"{old!r} is deprecated; use {new!r}",
                    DeprecationWarning,
                    stacklevel=2,
                )
                kwargs[new] = kwargs.pop(old)
            return func(*args, **kwargs)

        return wrapper

    return deco


def annotate_rt_izs(
    variants_df: pd.DataFrame, loc_cov_rt_iz: str | Path
) -> pd.DataFrame:
    """Annotate variants with 0 or 1 flag for RT IZs.

    Gives a 1 if the variant occurs in any replication-timing (RT)
    initiation zone (IZ) present in the supplied BED file.

    A new column 'cov_iz_<IZ>' is created for each distinct IZ (named
    <IZ>). The column is

    - 1 if the variant falls inside any interval belonging to
      that IZ
    - 0 if it is outside all intervals of that IZ

    - NaN if the chromosome is absent from the BED (e.g. chrX/chrY) or
      if the variant's coordinate is ``NaN``.

    Parameters
    ----------
    variants_df : pandas.DataFrame
        Variant catalogue. Must contain

        - 'Chromosome'
        - 'Start_Position'
        - optionally an index named 'variant' (preserved).

    loc_cov_rt_iz : str or pathlib.Path
        Path to the BED-like file that maps genomic intervals to IZ names.
        The file is read once and expected to have a column 'IZ' with
        the zone identifiers.

    Returns
    -------
    pandas.DataFrame
        `variants_df` with one extra column per initiation zone, each
        named ``cov_iz_<IZ>``.

    """
    izs = read_bed_file(loc_cov_rt_iz, "IZ")
    out_df = variants_df.copy()
    for iz in izs["IZ"].unique():
        logger.info(f"Covariate for RT initiation zone: {iz}")
        out_df = annotate_indicator_in_region(
            out_df, izs[izs["IZ"] == iz], f"iz_{iz}"
        )
    return out_df


def annotate_rt_left_right(
    variants_df: pd.DataFrame,
    loc_cov_rt_left: str | Path,
    loc_cov_rt_right: str | Path,
) -> pd.DataFrame:
    """Add indicator columns for leftward and rightward TTRs.

    Columns ``cov_left`` and ``cov_right`` are set to 1 when a variant
    falls inside the provided *left* or *right* regions respectively
    and 0 otherwise. Chromosomes absent from the BED files (e.g. chrX)
    or variants with missing coordinates receive ``NaN``.
    """
    out_df = variants_df.copy()

    logger.info("Covariate for RT leftward transition (TTR‐left)")
    out_df = annotate_indicator_in_region(
        out_df, read_bed_file(loc_cov_rt_left), "left"
    )

    logger.info("Covariate for RT rightward transition (TTR‐right)")
    out_df = annotate_indicator_in_region(
        out_df, read_bed_file(loc_cov_rt_right), "right"
    )

    return out_df


def annotate_rt_terms(
    variants_df: pd.DataFrame, loc_cov_rt_terms: str | Path
) -> pd.DataFrame:
    """Add ``cov_terms`` indicating proximity to RT termination sites."""
    logger.info("Covariate for RT termination sites (TTR‐terms)")
    return annotate_indicator_in_region(
        variants_df.copy(),
        read_bed_file(loc_cov_rt_terms, has_index_col=True),
        "terms",
    )


def annotate_rt_twidth(
    variants_df: pd.DataFrame,
    loc_cov_rt_twidth: str | Path,
    *,
    bin_size: int = 50000,
) -> pd.DataFrame:
    """Bin *Twidth* values and annotate variants with ``cov_twidth``.

    The BED‐like file at *loc_cov_rt_twidth* must contain a numeric column
    (conventionally called ``cov_twidth``) holding the *Twidth* value for
    each interval.
    """
    logger.info("Covariate for RT Twidth (binned)")
    return annotate_with_binned_features(
        variants_df.copy(),
        read_bed_file(loc_cov_rt_twidth, "cov_twidth"),
        "cov_twidth",
        bin_size=bin_size,
    )


def load_or_generate_rt_cov_df(
    location_df: str | Path,
    variants_df: pd.DataFrame,
    *,
    loc_cov_rt_izs: str | Path | None = None,
    loc_cov_rt_left: str | Path | None = None,
    loc_cov_rt_right: str | Path | None = None,
    loc_cov_rt_terms: str | Path | None = None,
    loc_cov_rt_twidth: str | Path | None = None,
    extra_cols_to_keep: list | None = None,
    force_generation: bool = False,
) -> pd.DataFrame:
    """Load or create a DataFrame containing only covariate columns.

    Parameters
    ----------
    location_df : str or pathlib.Path
        Where the pickled covariate DataFrame is or will be stored.
    variants_df : pandas.DataFrame
        Must index variants by a column or index named ``'variant'`` and
        include ``'Chromosome'`` and ``'Start_Position'``.
    location_* : str or pathlib.Path, optional
        Paths for each RT feature.  Omit or set to ``None`` to skip the
        corresponding annotation.
    extra_cols_to_keep : list, optional
        Keep these columns too besides the cov columns.
    force_generation : bool, default False
        Rebuild the covariate DataFrame even if *location_df* exists.
    """
    if os.path.exists(location_df) and not force_generation:
        logger.info(f"Loading covariates from {location_df}")
        cov_df = pd.read_pickle(location_df)
        return cov_df

    logger.info("Generating covariate DataFrame...")
    covariates = variants_df.copy()

    if loc_cov_rt_izs is not None:
        covariates = annotate_rt_izs(covariates, loc_cov_rt_izs)

    if loc_cov_rt_left is not None and loc_cov_rt_right is not None:
        covariates = annotate_rt_left_right(
            covariates, loc_cov_rt_left, loc_cov_rt_right
        )

    if loc_cov_rt_terms is not None:
        covariates = annotate_rt_terms(covariates, loc_cov_rt_terms)

    if loc_cov_rt_twidth is not None:
        covariates = annotate_rt_twidth(covariates, loc_cov_rt_twidth)

    if extra_cols_to_keep is None:
        extra_cols_to_keep = []

    def keep_cov_cols_only(
        df: pd.DataFrame, extra_cols: list
    ) -> pd.DataFrame:
        """Return only the variant index and *cov_* columns."""
        cov_cols = [c for c in df.columns if c.startswith("cov_")]
        return df.loc[:, extra_cols + cov_cols].copy()

    cov_df = keep_cov_cols_only(covariates, extra_cols_to_keep)
    cov_df.to_pickle(location_df)
    logger.info(f"Saved covariates to {location_df}")
    logger.info("... done.")
    return cov_df


def print_info_about_izs(variants_df):
    import numpy as np

    iz_cols = [
        c for c in variants_df.columns if c.startswith("cov_iz_")
    ]

    # ------------------------------------------------------------------
    # 1. For every variant say whether it is inside *any* IZ
    # ------------------------------------------------------------------
    in_any_iz = variants_df[iz_cols].any(
        axis=1
    )  # True / False per row

    print(
        "Variants in at least one IZ: "
        f"{len(variants_df[in_any_iz])} out of {len(variants_df)} "
        f"({round(100*len(variants_df[in_any_iz])/len(variants_df))}%)"
    )

    # ------------------------------------------------------------------
    # 2. For every gene ask: “do *any* of its variants hit an IZ?”
    # ------------------------------------------------------------------
    has_iz_per_gene = in_any_iz.groupby(variants_df["gene"]).any()

    print(
        "\nGenes in at least one IZ: "
        f"{sum(has_iz_per_gene)} out of {variants_df['gene'].nunique()} "
        f"({round(100*sum(has_iz_per_gene)/variants_df['gene'].nunique())}%)"
    )

    # ------------------------------------------------------------------
    # 3. For every variant pick the IZ column that is 1 (which is the
    # max), rows with all-zero flags → NaN so they don’t count
    # ------------------------------------------------------------------
    iz_of_variant = variants_df[iz_cols].idxmax(axis=1)
    iz_of_variant[variants_df[iz_cols].sum(axis=1) == 0] = np.nan

    # ------------------------------------------------------------------
    # 4.  Genes that have *both* IZ and non-IZ variants
    # ------------------------------------------------------------------
    has_non_iz_per_gene = (
        (~in_any_iz).groupby(variants_df["gene"]).any()
    )
    mixed_genes = has_iz_per_gene & has_non_iz_per_gene  # both True

    print(
        "\nGenes with variants inside *and* outside IZs: "
        f"{mixed_genes.sum()} out of {variants_df['gene'].nunique()} "
        f"({round(100*mixed_genes.sum()/variants_df['gene'].nunique())}%)"
    )
    # if mixed_genes.any():
    #     print(", ".join(mixed_genes[mixed_genes].index.tolist()))

    # ------------------------------------------------------------------
    # 4. Now, per gene, count distinct IZs (NaN ignored)
    # ------------------------------------------------------------------
    n_iz_per_gene = iz_of_variant.groupby(
        variants_df["gene"]
    ).nunique()

    # genes with variants in ≥ 2 different IZs
    genes_with_multi_iz = n_iz_per_gene[n_iz_per_gene >= 2]
    print("\nReturning genes in multiple IZs:")
    return genes_with_multi_iz


_AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))


def _bin_bigwigs(
    bigwig_paths: Sequence[str | Path],
    *,
    bin_size: int = 50_000,
    chromosomes: Sequence[str] | None = None,
    statistic: str = "mean",
) -> pd.DataFrame:
    """Bin one or more bigWig tracks into fixed genomic windows.

    Bins are 0-based half-open windows of `bin_size` starting at 0 on
    each chromosome (the convention of the Repli-seq MAT tables and
    of `annotate_with_binned_features`).  The value of track *i* goes
    into column ``track_i`` (1-based, in input order).

    Chromosome naming ('chr1' vs '1') is reconciled per track via
    `normalize_chromosome_name`; output names are the requested ones
    (default 'chr1'..'chr22').  When tracks disagree on a chromosome
    length, the minimum is used and a warning is logged.
    """
    import pyBigWig

    handles = [pyBigWig.open(str(p)) for p in bigwig_paths]
    try:
        wanted = (
            list(chromosomes)
            if chromosomes is not None
            else list(_AUTOSOMES)
        )
        chrom_col: list[str] = []
        start_col: list[int] = []
        end_col: list[int] = []
        values: list[list[float]] = [[] for _ in handles]

        for chrom in wanted:
            names = [
                normalize_chromosome_name(chrom, bw.chroms())
                for bw in handles
            ]
            lengths = [
                bw.chroms()[name]
                for bw, name in zip(handles, names)
                if name is not None
            ]
            if not lengths:
                logger.warning(
                    "Chromosome %s not found in any track; "
                    "skipping",
                    chrom,
                )
                continue
            if len(set(lengths)) > 1:
                logger.warning(
                    "Chromosome %s length differs across tracks "
                    "(%s); using the minimum",
                    chrom,
                    sorted(set(lengths)),
                )
            length = min(lengths)
            for start in range(0, length, bin_size):
                end = min(start + bin_size, length)
                chrom_col.append(chrom)
                start_col.append(start)
                end_col.append(end)
                for i, (bw, name) in enumerate(zip(handles, names)):
                    values[i].append(
                        fetch_bigwig_stat(
                            bw, name, start, end, statistic
                        )
                        if name is not None
                        else float("nan")
                    )
    finally:
        for bw in handles:
            bw.close()

    out = pd.DataFrame(
        {
            "Chromosome": chrom_col,
            "region_start": start_col,
            "region_end": end_col,
        }
    )
    for i, vals in enumerate(values, start=1):
        out[f"track_{i}"] = vals
    return out


def load_repliseq_fractions_bins_from_bigwigs(
    bigwig_paths: Sequence[str | Path],
    *,
    bin_size: int = 50_000,
    chromosomes: Sequence[str] | None = None,
    statistic: str = "mean",
) -> pd.DataFrame:
    """Bin N Repli-seq fraction bigWigs into the fractions table.

    The tracks must be given in **early-to-late** order (e.g. the UW
    6-fraction Repli-seq: G1b, S1, S2, S3, S4, G2).  Returns the same
    schema as :func:`load_repliseq_fractions_bins`:

    - 'Chromosome', 'region_start', 'region_end'
    - 'rt_s1'..'rt_sN' (float in 0..1; NaN if a bin has no signal)

    Each bin's fraction vector is normalized to sum to 1, so tracks
    that are already percentage-normalized (like the UW ones) and raw
    signal tracks are both handled.

    Parameters
    ----------
    bigwig_paths : sequence of paths
        One bigWig per S-phase fraction, earliest first.
    bin_size : int
        Genomic window size (default 50 kb, matching the MAT tables).
    chromosomes : sequence of str | None
        Chromosomes to bin; default autosomes chr1..chr22 (matching
        `load_gene_bodies_from_gtf(autosomes_only=True)` downstream).
    statistic : str
        Per-bin summary statistic (``pyBigWig.stats`` type).
    """
    if len(bigwig_paths) < 2:
        raise ValueError(
            "Need at least two fraction bigWigs; got "
            f"{len(bigwig_paths)} (for a single wavelet "
            "track see generate_rt_wavelet_per_gene)"
        )
    bins = _bin_bigwigs(
        bigwig_paths,
        bin_size=bin_size,
        chromosomes=chromosomes,
        statistic=statistic,
    )
    frac_cols = [
        f"track_{i}" for i in range(1, len(bigwig_paths) + 1)
    ]
    return _normalize_fraction_bins(bins, frac_cols)


def _normalize_fraction_bins(
    bins: pd.DataFrame, frac_cols: list[str]
) -> pd.DataFrame:
    """Normalize per-bin fraction signal to rt_s1..rt_sN on [0, 1].

    Each bin's fraction vector is scaled so it sums to 100 (per Zhao
    et al., Genome Biology 2020) and then divided by 100.  Bins with
    zero or missing total signal are returned as NaN.

    Parameters
    ----------
    bins : pd.DataFrame
        Must contain 'Chromosome', 'region_start', 'region_end' and
        the columns listed in `frac_cols` (raw fraction signal,
        earliest S-phase fraction first).
    frac_cols : list[str]
        The fraction columns, in early-to-late order.

    Returns
    -------
    pd.DataFrame
        'Chromosome', 'region_start', 'region_end',
        'rt_s1'..'rt_sN' (float in 0..1; NaN if bin has no signal).
    """
    out_cols = [f"rt_s{x}" for x in range(1, len(frac_cols) + 1)]

    row_sum = bins[frac_cols].sum(axis=1, min_count=1)
    scale = (100.0 / row_sum).where(row_sum > 0)
    scaled = bins[frac_cols].mul(scale, axis=0)

    out = bins[["Chromosome", "region_start", "region_end"]].copy()
    for src, dst in zip(frac_cols, out_cols):
        out[dst] = scaled[src].where(row_sum > 0).astype(float) / 100
    return out


def _resolve_source_type(source, source_type: str) -> str:
    """Resolve 'auto' to a concrete repliseq source type."""
    if source_type != "auto":
        if source_type not in ("mat", "fraction_bigwigs", "wavelet"):
            raise ValueError(
                f"Unknown repliseq source_type {source_type!r}"
            )
        return source_type
    if isinstance(source, (list, tuple)):
        return "fraction_bigwigs"
    return "mat"


def load_repliseq_fractions_bins(
    source, *, source_type: str = "auto", bin_size: int = 50_000
) -> pd.DataFrame:
    """Load multi-fraction Repli-seq and return per-bin normalized fractions.

    Two source types are supported (``source_type='auto'`` infers
    from the argument):

    - ``mat``: a single *wide, transposed* Repli-seq table (Zhao et
      al., 2020); ``source`` is its path.
    - ``fraction_bigwigs``: N per-fraction bigWig tracks in
      early-to-late order; ``source`` is a list/tuple of paths (see
      :func:`load_repliseq_fractions_bins_from_bigwigs`).

    Each bin is column-normalized so its S-phase fractions sum to 100
    and expressed on a [0, 1] scale (dividing by 100).

    Parameters
    ----------
    source : str | Path | Sequence[str | Path]
        Path to the transposed Repli-seq file, or the fraction bigWig
        paths. For the ``mat`` type, expected columns after parsing:

        - 'Chromosome', 'region_start', 'region_end'
        - fraction columns (any count ≥ 1).
    source_type : str
        'auto' (default), 'mat', or 'fraction_bigwigs'.
    bin_size : int
        Genomic window size for the bigWig path (default 50 kb).

    Returns
    -------
    pandas.DataFrame
        Columns:
        - 'Chromosome' (str)
        - 'region_start' (int, 0-based)
        - 'region_end'   (int, exclusive)
        - 'rt_s1'..'rt_sN' (float in 0..1; NaN if bin has no signal)

    Notes
    -----
    - Each bin's fraction vector F is scaled so sum(F) = 100 before
      dividing by 100 (per Zhao et al., Genome Biology 2020).
    - Bins with zero or missing total signal are returned as NaN.
    - 'rt_s1' is the earliest S-phase fraction; 'rt_sN' is the latest.

    """
    resolved = _resolve_source_type(source, source_type)
    if resolved == "fraction_bigwigs":
        return load_repliseq_fractions_bins_from_bigwigs(
            source, bin_size=bin_size
        )
    if resolved == "wavelet":
        raise ValueError(
            "Fractions are undefined for a wavelet track; use "
            "generate_rt_wavelet_per_gene instead"
        )

    repli_seq = read_bed_file(
        source,
        feature_name=None,
        has_index_col=False,
        has_header=False,
        file_is_transposed=True,
    )

    n_phases = len(repli_seq.columns) - 3
    frac_cols = [
        f"fraction_signal_s{x}" for x in range(1, n_phases + 1)
    ]
    repli_seq.columns = list(repli_seq.columns[:3]) + frac_cols

    num_cols = ["region_start", "region_end"] + frac_cols
    repli_seq[num_cols] = repli_seq[num_cols].apply(
        pd.to_numeric, errors="coerce"
    )

    return _normalize_fraction_bins(repli_seq, frac_cols)


def load_repliseq_mrt_bins(
    source,
    *,
    source_type: str = "auto",
    bin_size: int = 50_000,
    mrt_fraction_cols: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Load multi-fraction Repli-seq and compute per-bin MRT.

    Calls :func:`load_repliseq_fractions_bins`, then collapses the
    normalized fractions into a single mean replication time (MRT) via
    a weighted average over S-phase midpoints.  Output MRT is on
    [0, 1] where 0 is earliest and 1 is latest replication.

    Parameters
    ----------
    source : str | Path | Sequence[str | Path]
        Transposed Repli-seq file, or fraction bigWig paths (see
        :func:`load_repliseq_fractions_bins`).
    source_type : str
        'auto' (default), 'mat', or 'fraction_bigwigs'.
    bin_size : int
        Genomic window size for the bigWig path (default 50 kb).
    mrt_fraction_cols : Sequence[str] | None
        Restrict the weighted average to these ``rt_s*`` columns
        (renormalizing over the subset).  None (default) uses all
        fractions — for the UW 6-fraction data that reproduces the
        canonical weighted average of Hansen et al. (2010).

    Returns
    -------
    pandas.DataFrame
        Columns:
        - 'Chromosome' (str)
        - 'region_start' (int, 0-based)
        - 'region_end'   (int, exclusive)
        - 'mrt'       (float in 0..1; NaN if a bin has no signal)

    Notes
    -----
    - MRT is computed as sum(F * t) / sum(F), where t are S-phase
      midpoints ( (i+0.5)/N , i=0..N-1 ).
    - Bins with zero or missing total signal are returned as NaN.
    - For the individual per-fraction columns see
      :func:`load_repliseq_fractions_bins`.

    """
    fracs = load_repliseq_fractions_bins(
        source, source_type=source_type, bin_size=bin_size
    )
    rt_cols = [c for c in fracs.columns if c.startswith("rt_s")]
    if mrt_fraction_cols is not None:
        missing = [c for c in mrt_fraction_cols if c not in rt_cols]
        if missing:
            raise ValueError(
                f"mrt_fraction_cols not in the fractions: {missing}; "
                f"available: {rt_cols}"
            )
        rt_cols = list(mrt_fraction_cols)
    n_phases = len(rt_cols)
    t = (np.arange(n_phases, dtype=float) + 0.5) / n_phases
    F = fracs[rt_cols].to_numpy(dtype=float)
    has_signal = np.isfinite(F).any(axis=1)
    if mrt_fraction_cols is None:
        # full set: fractions already sum to 1 per bin
        mrt = np.nansum(F * t, axis=1)
    else:
        # subset: renormalize over the selected fractions
        denom = np.nansum(F, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mrt = np.nansum(F * t, axis=1) / denom
        has_signal &= denom > 0

    out = fracs[["Chromosome", "region_start", "region_end"]].copy()
    out["mrt"] = pd.Series(mrt, index=fracs.index).where(has_signal)
    return out


@_legacy_kwarg("repli_seq_hct", "repliseq_source")
def generate_rt_fractions_per_gene(
    repliseq_source,
    gencode_annotation: str | Path,
    *,
    source_type: str = "auto",
    bin_size: int = 50_000,
) -> pd.DataFrame:
    """Compute gene-level Repli-seq fractions from multi-fraction data.

    Returns one column per S-phase fraction (``rt_s1``…``rt_sN``),
    each being a length-weighted average of the normalized bin values
    over the gene body.  Use this to add individual replication-timing
    fractions as separate covariates in a model.

    For a scalar summary see :func:`generate_mrt_per_gene`.

    Parameters
    ----------
    repliseq_source : str | Path
        Path to the transposed Repli-seq file for the cell line
        (e.g., HCT116), compatible with
        :func:`load_repliseq_fractions_bins`.
    gencode_annotation : str | Path
        Path to a GENCODE/Ensembl GTF (GRCh38 to match hg38 bins).

    Returns
    -------
    pd.DataFrame
        Gene-level fraction values on 0..1 scale, indexed by
        ``ensembl_gene_id``.  Columns are ``rt_s1``…``rt_sN``.
        Genes with insufficient bin coverage yield NaN rows.

    """
    fracs = load_repliseq_fractions_bins(
        repliseq_source, source_type=source_type, bin_size=bin_size
    )
    gene_bodies = load_gene_bodies_from_gtf(gencode_annotation)
    rt_cols = [c for c in fracs.columns if c.startswith("rt_s")]
    annotated = annotate_with_binned_features(
        gene_bodies, fracs, feature_cols=rt_cols
    )
    return annotated[rt_cols]


@_legacy_kwarg("repli_seq_hct", "repliseq_source")
def load_or_generate_rt_fractions(
    location_csv: str | Path,
    repliseq_source,
    gencode_annotation: str | Path,
    *,
    source_type: str = "auto",
    bin_size: int = 50_000,
    force_generation: bool = False,
    float_format: str = "%.6g",
) -> pd.DataFrame:
    """Load or generate gene-level Repli-seq fractions.

    If the CSV exists at ``location_csv`` and ``force_generation`` is
    False, load it; otherwise call
    :func:`generate_rt_fractions_per_gene` and cache the result.

    Parameters
    ----------
    location_csv : str | Path
        Path to the CSV to read/write.  Append ``.gz`` for transparent
        compression.
    repliseq_source : str | Path
        Path to the transposed multi-fraction Repli-seq file.
    gencode_annotation : str | Path
        Path to a GENCODE/Ensembl GTF (hg38/GRCh38).
    force_generation : bool
        Regenerate even if the CSV already exists.
    float_format : str
        Format for writing floats, default ``'%.6g'``.

    Returns
    -------
    pd.DataFrame
        Index: ``ensembl_gene_id``; columns: ``rt_s1``…``rt_sN``
        (normalized S-phase fractions on 0..1 scale).

    """
    location_csv = Path(location_csv)

    if location_csv.exists() and not force_generation:
        logger.info(
            "Loading RT fractions per gene from %s", location_csv
        )
        df = pd.read_csv(location_csv, index_col=0)
        df.index.name = "ensembl_gene_id"
        logger.info("... done loading RT fractions per gene.")
        return df.astype(float)

    logger.info(
        "Generating RT fractions per gene from %s and %s",
        repliseq_source,
        gencode_annotation,
    )
    df = generate_rt_fractions_per_gene(
        repliseq_source,
        gencode_annotation,
        source_type=source_type,
        bin_size=bin_size,
    )
    df.to_csv(location_csv, float_format=float_format)
    logger.info("Saved RT fractions per gene to %s", location_csv)
    logger.info("... done generating RT fractions per gene.")
    return df


def generate_rt_wavelet_per_gene(
    bigwig_path: str | Path,
    gencode_annotation: str | Path,
    *,
    bin_size: int = 50_000,
    chromosomes: Sequence[str] | None = None,
    statistic: str = "mean",
) -> pd.Series:
    """Compute a per-gene RT value from a single smoothed signal track.

    For cell lines whose only processed replication-timing output is
    one wavelet-smoothed early/late signal bigWig (the Gilbert-lab
    ENCODE series: LNCaP, Caki2, A549, NCI-H460, G401, SK-N-MC,
    T47D), there are no fractions to normalize — the track itself is
    the covariate.  The signal is binned and length-weight averaged
    over gene bodies, like MRT.

    Note the orientation is the track's own (for log2(early/late)
    signal, larger = **earlier** replication — opposite of ``mrt``).
    As a model covariate the sign is absorbed by the coefficient.

    Returns
    -------
    pd.Series
        Indexed by 'ensembl_gene_id'; name 'rt_wavelet'.
    """
    bins = _bin_bigwigs(
        [bigwig_path],
        bin_size=bin_size,
        chromosomes=chromosomes,
        statistic=statistic,
    )
    bins = bins.rename(columns={"track_1": "rt_wavelet"})
    gene_bodies = load_gene_bodies_from_gtf(gencode_annotation)
    annotated = annotate_with_binned_features(
        gene_bodies, bins, feature_cols="rt_wavelet"
    )
    return annotated["rt_wavelet"].astype(float).rename("rt_wavelet")


def load_or_generate_rt_wavelet(
    location_csv: str | Path,
    bigwig_path: str | Path,
    gencode_annotation: str | Path,
    *,
    bin_size: int = 50_000,
    force_generation: bool = False,
    float_format: str = "%.6g",
) -> pd.Series:
    """Load or generate the per-gene wavelet RT covariate.

    Same caching contract as :func:`load_or_generate_mrt`.
    """
    location_csv = Path(location_csv)

    if location_csv.exists() and not force_generation:
        logger.info(
            "Loading wavelet RT per gene from %s", location_csv
        )
        tbl = pd.read_csv(location_csv, index_col=0)
        ser = tbl.iloc[:, 0].astype(float).rename("rt_wavelet")
        logger.info("... done loading wavelet RT per gene.")
        return ser

    logger.info(
        "Generating wavelet RT per gene from %s and %s",
        bigwig_path,
        gencode_annotation,
    )
    ser = generate_rt_wavelet_per_gene(
        bigwig_path, gencode_annotation, bin_size=bin_size
    )
    ser.to_frame().to_csv(location_csv, float_format=float_format)
    logger.info("Saved wavelet RT per gene to %s", location_csv)
    logger.info("... done generating wavelet RT per gene.")
    return ser


@_legacy_kwarg("repli_seq_hct", "repliseq_source")
def generate_mrt_per_gene(
    repliseq_source,
    gencode_annotation,
    *,
    source_type: str = "auto",
    bin_size: int = 50_000,
    mrt_fraction_cols: Sequence[str] | None = None,
):
    """Compute gene-level MRT from multi-fraction Repli-seq.

    This function aggregates per-bin mean replication timing (MRT)
    into a single value for each gene by taking a length-weighted
    average of bin values over the gene body.

    Parameters
    ----------
    repliseq_source : str | Path
        Path to the transposed Repli-seq file for the cell line (e.g., HCT116),
        compatible with `load_repliseq_mrt_bins`.
    gencode_annotation : str | Path
        Path to a GENCODE/Ensembl GTF (GRCh38 to match hg38 bins).

    Returns
    -------
    pd.Series
        Gene-level MRT on 0..1 scale, indexed by 'ensembl_gene_id'.
        Larger values indicate later replication.
        Name: 'rt_mrt'.

    Notes
    -----
    - Coordinate convention is 0-based half-open for both bins and
      gene bodies.
    - Assembly should match between the Repli-seq bins and the GTF
      (e.g., hg38).
    - Orientation: if you prefer "earliness", use `1 - result`.
    - Genes with insufficient or missing bin coverage yield NaN.

    """
    cov_mrt = load_repliseq_mrt_bins(
        repliseq_source,
        source_type=source_type,
        bin_size=bin_size,
        mrt_fraction_cols=mrt_fraction_cols,
    )

    gene_bodies = load_gene_bodies_from_gtf(gencode_annotation)

    mrt_per_gene = annotate_with_binned_features(
        gene_bodies, cov_mrt
    )["mrt"]

    return mrt_per_gene


@_legacy_kwarg("repli_seq_hct", "repliseq_source")
def load_or_generate_mrt(
    location_csv: str | Path,
    repliseq_source,
    gencode_annotation: str | Path,
    *,
    source_type: str = "auto",
    bin_size: int = 50_000,
    mrt_fraction_cols: Sequence[str] | None = None,
    force_generation: bool = False,
    float_format: str = "%.6g",
) -> pd.Series:
    """Load or generate gene-level MRT (mean replication time, 0..1).

    If the CSV exists at `location_csv` and `force_generation` is
    False, load it.

    Parameters
    ----------
    location_csv : str | Path
        Path to the CSV to read/write (e.g., 'mrt_per_gene.csv').
        If the filename ends with '.gz', pandas will transparently compress.
    repliseq_source : str | Path
        Path to the transposed multi-fraction Repli-seq file for the cell line.
    gencode_annotation : str | Path
        Path to a GENCODE/Ensembl GTF (hg38/GRCh38 to match the
        Repli-seq bins).
    force_generation : bool
        If True, regenerate even if the CSV already exists.
    float_format : str
        Format for writing floats to CSV, default '%.6g'.

    Returns
    -------
    pd.Series
        Index: `ensembl_gene_id`; values: MRT in 0..1 (larger = later
        replication).
        Name: 'mrt'.

    Notes
    -----
    - Coordinate system is 0-based half-open for both bins and gene bodies.
    - `load_gene_bodies_from_gtf` defaults will add 'chr' and keep autosomes.
    - If you prefer “earliness”, transform later: `1 - returned_series`.

    """
    location_csv = Path(location_csv)

    if location_csv.exists() and not force_generation:
        logger.info("Loading MRT per gene from %s", location_csv)
        tbl = pd.read_csv(location_csv, index_col=0)
        if tbl.shape[1] == 1:
            ser = tbl.iloc[:, 0]
        else:
            col = "mrt" if "mrt" in tbl.columns else tbl.columns[0]
            ser = tbl[col]
        ser = ser.astype(float).rename("mrt")
        logger.info("... done loading MRT per gene.")
        return ser

    logger.info(
        f"Generating MRT per gene from {repliseq_source} "
        f"and {gencode_annotation}"
    )

    cov_mrt = load_repliseq_mrt_bins(
        repliseq_source,
        source_type=source_type,
        bin_size=bin_size,
        mrt_fraction_cols=mrt_fraction_cols,
    )
    gene_bodies = load_gene_bodies_from_gtf(gencode_annotation)
    annotated = annotate_with_binned_features(
        gene_bodies, cov_mrt, feature_cols="mrt"
    )
    ser = annotated["mrt"].astype(float).rename("mrt")
    ser.to_frame().to_csv(location_csv, float_format=float_format)
    logger.info("Saved MRT per gene to %s", location_csv)
    logger.info("... done generating MRT per gene.")
    return ser
