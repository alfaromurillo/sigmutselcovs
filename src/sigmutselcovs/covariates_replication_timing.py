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

from .covariates_utilities import read_bed_file

from .covariates_utilities import annotate_indicator_in_region
from .covariates_utilities import annotate_with_binned_features

from .covariates_utilities import load_gene_bodies_from_gtf

import logging


logger = logging.getLogger(__name__)


def annotate_rt_izs(variants_df: pd.DataFrame,
                    loc_cov_rt_iz: str | Path) -> pd.DataFrame:
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
    for iz in izs['IZ'].unique():
        logger.info(f"Covariate for RT initiation zone: {iz}")
        out_df = annotate_indicator_in_region(
            out_df,
            izs[izs['IZ'] == iz],
            f"iz_{iz}")
    return out_df


def annotate_rt_left_right(
        variants_df: pd.DataFrame,
        loc_cov_rt_left: str | Path,
        loc_cov_rt_right: str | Path) -> pd.DataFrame:
    """Add indicator columns for leftward and rightward TTRs.

    Columns ``cov_left`` and ``cov_right`` are set to 1 when a variant
    falls inside the provided *left* or *right* regions respectively
    and 0 otherwise. Chromosomes absent from the BED files (e.g. chrX)
    or variants with missing coordinates receive ``NaN``.
    """
    out_df = variants_df.copy()

    logger.info("Covariate for RT leftward transition (TTR‐left)")
    out_df = annotate_indicator_in_region(
        out_df,
        read_bed_file(loc_cov_rt_left),
        "left")

    logger.info("Covariate for RT rightward transition (TTR‐right)")
    out_df = annotate_indicator_in_region(
        out_df,
        read_bed_file(loc_cov_rt_right),
        "right")

    return out_df


def annotate_rt_terms(
        variants_df: pd.DataFrame,
        loc_cov_rt_terms: str | Path) -> pd.DataFrame:
    """Add ``cov_terms`` indicating proximity to RT termination sites."""
    logger.info("Covariate for RT termination sites (TTR‐terms)")
    return annotate_indicator_in_region(
        variants_df.copy(),
        read_bed_file(loc_cov_rt_terms, has_index_col=True),
        "terms")


def annotate_rt_twidth(
        variants_df: pd.DataFrame,
        loc_cov_rt_twidth: str | Path,
        *,
        bin_size: int = 50000) -> pd.DataFrame:
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
        bin_size=50000)


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
        force_generation: bool = False) -> pd.DataFrame:
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
            covariates,
            loc_cov_rt_left,
            loc_cov_rt_right)

    if loc_cov_rt_terms is not None:
        covariates = annotate_rt_terms(covariates, loc_cov_rt_terms)

    if loc_cov_rt_twidth is not None:
        covariates = annotate_rt_twidth(covariates, loc_cov_rt_twidth)

    if extra_cols_to_keep is None:
        extra_cols_to_keep = []

    def keep_cov_cols_only(df: pd.DataFrame,
                           extra_cols: list) -> pd.DataFrame:
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
    iz_cols = [c for c in variants_df.columns if c.startswith("cov_iz_")]

    # ------------------------------------------------------------------
    # 1. For every variant say whether it is inside *any* IZ
    # ------------------------------------------------------------------
    in_any_iz = variants_df[iz_cols].any(axis=1)     # True / False per row

    print("Variants in at least one IZ: "
          f"{len(variants_df[in_any_iz])} out of {len(variants_df)} "
          f"({round(100*len(variants_df[in_any_iz])/len(variants_df))}%)")

    # ------------------------------------------------------------------
    # 2. For every gene ask: “do *any* of its variants hit an IZ?”
    # ------------------------------------------------------------------
    has_iz_per_gene = in_any_iz.groupby(variants_df["gene"]).any()

    print(
        "\nGenes in at least one IZ: "
        f"{sum(has_iz_per_gene)} out of {variants_df['gene'].nunique()} "
        f"({round(100*sum(has_iz_per_gene)/variants_df['gene'].nunique())}%)")

    # ------------------------------------------------------------------
    # 3. For every variant pick the IZ column that is 1 (which is the
    # max), rows with all-zero flags → NaN so they don’t count
    # ------------------------------------------------------------------
    iz_of_variant = (
        variants_df[iz_cols].idxmax(axis=1))
    iz_of_variant[variants_df[iz_cols].sum(axis=1) == 0] = np.nan

    # ------------------------------------------------------------------
    # 4.  Genes that have *both* IZ and non-IZ variants
    # ------------------------------------------------------------------
    has_non_iz_per_gene = (~in_any_iz).groupby(variants_df["gene"]).any()
    mixed_genes = has_iz_per_gene & has_non_iz_per_gene      # both True

    print(
        "\nGenes with variants inside *and* outside IZs: "
        f"{mixed_genes.sum()} out of {variants_df['gene'].nunique()} "
        f"({round(100*mixed_genes.sum()/variants_df['gene'].nunique())}%)")
    # if mixed_genes.any():
    #     print(", ".join(mixed_genes[mixed_genes].index.tolist()))

    # ------------------------------------------------------------------
    # 4. Now, per gene, count distinct IZs (NaN ignored)
    # ------------------------------------------------------------------
    n_iz_per_gene = iz_of_variant.groupby(variants_df["gene"]).nunique()

    # genes with variants in ≥ 2 different IZs
    genes_with_multi_iz = n_iz_per_gene[n_iz_per_gene >= 2]
    print("\nReturning genes in multiple IZs:")
    return genes_with_multi_iz


def load_repliseq_fractions_bins(path: str) -> pd.DataFrame:
    """Load multi-fraction Repli-seq and return per-bin normalized fractions.

    Reads a *wide, transposed* Repli-seq table and column-normalizes
    each bin so its S-phase fractions sum to 100 (Zhao et al., 2020).
    Each fraction column is then expressed on a [0, 1] scale (dividing
    by 100).

    Parameters
    ----------
    path : str
        Path to the transposed Repli-seq file. Expected columns after
        parsing:

        - 'Chromosome', 'region_start', 'region_end'
        - fraction columns (any count ≥ 1).

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
    repli_seq = read_bed_file(path,
                              feature_name=None,
                              has_index_col=False,
                              has_header=False,
                              file_is_transposed=True)

    n_phases = len(repli_seq.columns) - 3
    frac_cols = [f"fraction_signal_s{x}" for x in range(1, n_phases + 1)]
    out_cols = [f"rt_s{x}" for x in range(1, n_phases + 1)]
    repli_seq.columns = list(repli_seq.columns[:3]) + frac_cols

    num_cols = ["region_start", "region_end"] + frac_cols
    repli_seq[num_cols] = repli_seq[num_cols].apply(pd.to_numeric,
                                                    errors="coerce")

    row_sum = repli_seq[frac_cols].sum(axis=1, min_count=1)
    scale = (100.0 / row_sum).where(row_sum > 0)
    scaled = repli_seq[frac_cols].mul(scale, axis=0)

    out = repli_seq[["Chromosome", "region_start", "region_end"]].copy()
    for src, dst in zip(frac_cols, out_cols):
        out[dst] = scaled[src].where(row_sum > 0).astype(float) / 100
    return out


def load_repliseq_mrt_bins(path: str) -> pd.DataFrame:
    """Load multi-fraction Repli-seq and compute per-bin MRT.

    Calls :func:`load_repliseq_fractions_bins`, then collapses the
    normalized fractions into a single mean replication time (MRT) via
    a weighted average over S-phase midpoints.  Output MRT is on
    [0, 1] where 0 is earliest and 1 is latest replication.

    Parameters
    ----------
    path : str
        Path to the transposed Repli-seq file (same format as
        :func:`load_repliseq_fractions_bins`).

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
    fracs = load_repliseq_fractions_bins(path)
    rt_cols = [c for c in fracs.columns if c.startswith("rt_s")]
    n_phases = len(rt_cols)
    t = (np.arange(n_phases, dtype=float) + 0.5) / n_phases
    F = fracs[rt_cols].to_numpy(dtype=float)
    has_signal = np.isfinite(F).any(axis=1)
    mrt = np.nansum(F * t, axis=1)

    out = fracs[["Chromosome", "region_start", "region_end"]].copy()
    out["mrt"] = pd.Series(mrt, index=fracs.index).where(has_signal)
    return out


def generate_rt_fractions_per_gene(
        repli_seq_hct: str | Path,
        gencode_annotation: str | Path) -> pd.DataFrame:
    """Compute gene-level Repli-seq fractions from multi-fraction data.

    Returns one column per S-phase fraction (``rt_s1``…``rt_sN``),
    each being a length-weighted average of the normalized bin values
    over the gene body.  Use this to add individual replication-timing
    fractions as separate covariates in a model.

    For a scalar summary see :func:`generate_mrt_per_gene`.

    Parameters
    ----------
    repli_seq_hct : str | Path
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
    fracs = load_repliseq_fractions_bins(repli_seq_hct)
    gene_bodies = load_gene_bodies_from_gtf(gencode_annotation)
    rt_cols = [c for c in fracs.columns if c.startswith("rt_s")]
    annotated = annotate_with_binned_features(
        gene_bodies, fracs, feature_cols=rt_cols)
    return annotated[rt_cols]


def load_or_generate_rt_fractions(
        location_csv: str | Path,
        repli_seq_hct: str | Path,
        gencode_annotation: str | Path,
        *,
        force_generation: bool = False,
        float_format: str = "%.6g"
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
    repli_seq_hct : str | Path
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
        logger.info("Loading RT fractions per gene from %s", location_csv)
        df = pd.read_csv(location_csv, index_col=0)
        df.index.name = "ensembl_gene_id"
        logger.info("... done loading RT fractions per gene.")
        return df.astype(float)

    logger.info("Generating RT fractions per gene from %s and %s",
                repli_seq_hct, gencode_annotation)
    df = generate_rt_fractions_per_gene(repli_seq_hct, gencode_annotation)
    df.to_csv(location_csv, float_format=float_format)
    logger.info("Saved RT fractions per gene to %s", location_csv)
    logger.info("... done generating RT fractions per gene.")
    return df


def generate_mrt_per_gene(repli_seq_hct,
                          gencode_annotation):
    """Compute gene-level MRT from multi-fraction Repli-seq.

    This function aggregates per-bin mean replication timing (MRT)
    into a single value for each gene by taking a length-weighted
    average of bin values over the gene body.

    Parameters
    ----------
    repli_seq_hct : str | Path
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
    cov_mrt = load_repliseq_mrt_bins(repli_seq_hct)

    gene_bodies = load_gene_bodies_from_gtf(gencode_annotation)

    mrt_per_gene = annotate_with_binned_features(gene_bodies, cov_mrt)['mrt']

    return mrt_per_gene


def load_or_generate_mrt(
        location_csv: str | Path,
        repli_seq_hct: str | Path,
        gencode_annotation: str | Path,
        *,
        force_generation: bool = False,
        float_format: str = "%.6g"
        ) -> pd.Series:
    """Load or generate gene-level MRT (mean replication time, 0..1).

    If the CSV exists at `location_csv` and `force_generation` is
    False, load it.

    Parameters
    ----------
    location_csv : str | Path
        Path to the CSV to read/write (e.g., 'mrt_per_gene.csv').
        If the filename ends with '.gz', pandas will transparently compress.
    repli_seq_hct : str | Path
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

    logger.info(f"Generating MRT per gene from {repli_seq_hct} "
                f"and {gencode_annotation}")

    cov_mrt = load_repliseq_mrt_bins(repli_seq_hct)
    gene_bodies = load_gene_bodies_from_gtf(gencode_annotation)
    annotated = annotate_with_binned_features(
        gene_bodies,
        cov_mrt,
        feature_cols="mrt")
    ser = annotated["mrt"].astype(float).rename("mrt")
    ser.to_frame().to_csv(location_csv, float_format=float_format)
    logger.info("Saved MRT per gene to %s", location_csv)
    logger.info("... done generating MRT per gene.")
    return ser
