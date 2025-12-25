"""Utility functions for processing genomic covariate data.

This module provides functions to extract genomic features from
BigWig files and GTF annotations, including signal extraction for
gene bodies and promoter regions. These utilities are used to build
covariate matrices for mutation rate modeling.
"""

import numpy as np
import pandas as pd

from collections.abc import Sequence
from pathlib import Path

import logging

logger = logging.getLogger(__name__)


def sanitize_feature_label(label: str) -> str:
    """Return a lowercase, underscore-safe label for feature columns.

    Parameters
    ----------
    label : str
        Original name of the feature or track.

    Returns
    -------
    str
        Normalised label containing only lowercase alphanumerics and
        underscores.
    """
    cleaned = label.lower().replace(" ", "_").replace("-", "_")
    cleaned = cleaned.replace(".", "_")
    return "".join(ch for ch in cleaned if ch.isalnum() or ch == "_")


def normalize_chromosome_name(
        chrom: str,
        available: dict[str, int] | dict[str, tuple[int, int]]
        ) -> str | None:
    """Return a chromosome name present in *available* if possible.

    Parameters
    ----------
    chrom : str
        Chromosome identifier to normalise (e.g., ``'chr1'`` or ``'1'``).
    available : dict[str, int] | dict[str, tuple[int, int]]
        Mapping of chromosome names present in a resource such as a bigWig
        header.

    Returns
    -------
    str | None
        Matching chromosome name from *available*, or ``None`` when no match
        is found.
    """
    if chrom in available:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom.removeprefix("chr")
        if alt in available:
            return alt
    else:
        alt = f"chr{chrom}"
        if alt in available:
            return alt
    return None


def run_pca_on_covariates(
        cov_df: pd.DataFrame,
        columns: list[str] | None = None,
        n_components: int | None = None,
        *,
        standardize: bool = True,
        dropna: str = "any",
        **pca_kwargs) -> pd.DataFrame:
    """Compute PCA over gene-level covariates and return scores.

    Parameters
    ----------
    cov_df : pandas.DataFrame
        Gene-indexed covariates (index = ensembl_gene_id).
    columns : list[str] | None
        Subset of columns to include. If None, use all numeric cols.
    standardize : bool, default True
        If True, z-score features before PCA.
    dropna : {'any','all','none'}, default 'any'
        How to handle NaNs across selected columns:
        - 'any': drop rows with any NaN
        - 'all': drop rows with all NaN
        - 'none': fill remaining NaNs with column means
    n_components : int | None, default None
        Number of principal components. If None, PCA decides based on
        provided parameters (e.g., if n_samples is larger than the
        number of covariates, then the number of covariates).
    **pca_kwargs
        Extra keyword arguments forwarded to
        sklearn.decomposition.PCA (e.g., whiten=True,
        svd_solver='full').

    Returns
    -------
    scores : pandas.DataFrame
        Gene-indexed PC scores with columns 'PC1', 'PC2', ...

        The returned DataFrame contains PCA metadata in
        ``scores.attrs``:
        - ``explained_variance_ratio`` : ndarray of shape (k,)
          Fraction of total variance explained by each principal
          component.
        - ``components`` : ndarray of shape (k, n_features)
          Principal axes (loadings) in feature space.

    """
    from sklearn.decomposition import PCA

    if columns is None:
        # keep only numeric columns
        cols = [c for c in cov_df.columns
                if pd.api.types.is_numeric_dtype(cov_df[c])]
    else:
        cols = list(columns)

    X = cov_df[cols].copy()

    if dropna == "any":
        X = X.dropna(how="any")
    elif dropna == "all":
        X = X.dropna(how="all")
        X = X.fillna(X.mean(numeric_only=True))
    elif dropna == "none":
        X = X.fillna(X.mean(numeric_only=True))
    else:
        raise ValueError("dropna must be one of {'any','all','none'}")

    means = X.mean(axis=0)
    scales = X.std(axis=0, ddof=0)
    if standardize:
        denom = scales.replace(0.0, 1.0)
        X_proc = (X - means) / denom
    else:
        X_proc = X

    if n_components is not None:
        pca = PCA(n_components=n_components, **pca_kwargs)
    else:
        pca = PCA(**pca_kwargs)
    Z = pca.fit_transform(X_proc.values)
    k = Z.shape[1]
    score_cols = [f"PC{i+1}" for i in range(k)]
    scores = pd.DataFrame(Z, index=X.index, columns=score_cols)

    return scores


def read_bed_file(bed_file,
                  feature_name=None,
                  has_index_col=False,
                  has_header=False,
                  file_is_transposed=False):
    """Read a BED or BED-like tab-delimited file into a DataFrame.

    This function assumes the file contains at least three columns:
    chromosome, start, and end. Additional columns (e.g., scores,
    labels, metadata) are supported and can be automatically named or
    manually specified.

    Args:
        bed_file (str): Path to the BED or BED-like file.

        feature_name (str | list[str] | None): How to name additional
            columns beyond the first three (Chromosome, region_start,
            region_end). Options:
            - If None (default), generic names are assigned:
              "feature_1", "feature_2", etc.
            - If a string is provided, it's used as the name for the
              fourth column. If there are other feature columns they
              will be named feature_2, etc.
            - If a list is provided, it's used directly for column
              names starting from column 4. If there are other feature
              columns beyond the names provided, they will be named
              feature_n, feature_n+1, etc.

        has_index_col (bool): Whether the first column should be used
            as index. Default is False.

        has_header (bool): Whether the file includes a header row.
            Default is False.

        file_is_transposed (bool): Whether the file has instead of
            columns for features, rows.

    Returns:
        pd.DataFrame: BED-like DataFrame with columns:
            - 'Chromosome'
            - 'region_start'
            - 'region_end'
            - Additional columns (if present and renamed)

    Raises:
        ValueError: If file has fewer than 3 columns, or if
            feature_name is a list of incorrect length.

    """
    regions = pd.read_csv(bed_file,
                          sep="\t",
                          header=0 if has_header else None,
                          index_col=0 if has_index_col else None)
    if file_is_transposed:
        regions = regions.T

    if regions.shape[1] < 3:
        logger.error(f"{bed_file} has only {regions.shape[1]} columns; "
                     "at least 3 are required for BED format.")
        raise ValueError("Input file must have at least 3 columns: "
                         "Chromosome, region_start, and region_end.")

    col_names = list(regions.columns)
    col_names[0:3] = ["Chromosome", "region_start", "region_end"]

    col_names[3:len(col_names)+1] = [f"feature_{i + 1}"
                                     for i in range(len(col_names) - 3)]
    if isinstance(feature_name, str):
        col_names[3] = feature_name
    elif isinstance(feature_name, list):
        col_names[3:len(feature_name)+3] = feature_name
    elif feature_name is not None:
        error_info = "feature_name must be either a string, a list, or None."
        logger.error(error_info)
        raise TypeError(error_info)

    regions.columns = col_names

    return regions


def load_gene_bodies_from_gtf(
    gtf_path: str | Path,
    biotypes: list[str] | None = None,
    add_chr_prefix_if_needed: bool = True,
    autosomes_only: bool = True,
    drop_version_suffix: bool = True) -> pd.DataFrame:
    """Load gene bodies from a GTF into a BED-like DataFrame.

    Parses only `feature == "gene"` entries from a GTF (e.g., GENCODE),
    returning 0-based, half-open gene intervals. The resulting DataFrame
    is indexed by `ensembl_gene_id`.

    Parameters
    ----------
    gtf_path : str or pathlib.Path
        Path to a GTF file (plain text or `.gz`).
    biotypes : list[str] or None, optional
        If provided, keep only genes whose `gene_type`/`gene_biotype`
        is in this list (e.g., ["protein_coding"]).
    add_chr_prefix_if_needed : bool, default True
        If chromosomes in the GTF are like "1", "2", ... add "chr"
        so they become "chr1", "chr2", ... (to match typical RT tracks).
    autosomes_only : bool, default True
        Keep only autosomes "chr1".."chr22".
    drop_version_suffix : bool, default True
        If True, strip the version from Ensembl IDs (e.g., ".17" in
        "ENSG00000123456.17" → "ENSG00000123456").

    Returns
    -------
    pandas.DataFrame
        Columns: ["Chromosome", "start", "end", "strand"], indexed by
        `ensembl_gene_id`. Coordinates are 0-based, half-open [start, end).

    Notes
    -----
    - GTF coordinates are 1-based inclusive; we convert to 0-based,
      half-open (BED-like) by doing `start-1, end`.
    - If `drop_version_suffix` is True, versioned Ensembl IDs are
      de-versioned to facilitate joins with resources that omit versions.

    Raises
    ------
    KeyError
        If the GTF lacks required attributes like `gene_id`.
    """
    import re

    # Convert to Path if string
    gtf_path = Path(gtf_path)

    # Open plain or gz
    if gtf_path.suffix == ".gz":
        import gzip
        op = gzip.open
    else:
        op = open

    rows = []
    with op(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            chrom, src, feat, start, end, score, strand, frame, attrs = (
                line.rstrip("\n").split("\t"))
            if feat != "gene":
                continue

            # Parse attributes: gene_id and gene_type/biotype
            m_id = re.search(r'gene_id "([^"]+)"', attrs)
            if not m_id:
                # Malformed line; skip
                continue
            gene_id = m_id.group(1)
            if drop_version_suffix:
                gene_id = gene_id.split(".")[0]

            m_type = (re.search(r'gene_type "([^"]+)"', attrs)
                      or re.search(r'gene_biotype "([^"]+)"', attrs))
            gtype = m_type.group(1) if m_type else None
            if biotypes is not None and gtype not in set(biotypes):
                continue

            # 1-based inclusive -> 0-based half-open
            s0 = int(start) - 1
            e1 = int(end)
            rows.append((gene_id, chrom, s0, e1, strand))

    df = pd.DataFrame(
        rows, columns=["ensembl_gene_id",
                       "Chromosome",
                       "start",
                       "end",
                       "strand"])

    if df.empty:
        # Return an empty, correctly-shaped DataFrame with index set
        return df.set_index("ensembl_gene_id")

    # Ensure 'chr' style to match typical RT bins if requested
    if add_chr_prefix_if_needed:
        # Add 'chr' only if the current values do not already start with it
        if not df["Chromosome"].astype(str).str.startswith("chr").all():
            df["Chromosome"] = "chr" + df["Chromosome"].astype(str)

    if autosomes_only:
        autos = {f"chr{i}" for i in range(1, 23)}
        df = df[df["Chromosome"].isin(autos)]

    # Set index to Ensembl gene ID
    df = df.set_index("ensembl_gene_id").sort_index()

    return df


def annotate_with_binned_features(df: pd.DataFrame,
                                  binned_df: pd.DataFrame,
                                  feature_cols=None,
                                  bin_size: int | None = None) -> pd.DataFrame:
    """Annotate a DataFrame with fixed-bin genomic features.

    Modes
    -----
    1) Variant mode (point):
       `df` index name is 'variant' and columns
       {'Chromosome','Start_Position'} exist.
       → Each row gets the feature(s) from the bin containing
       Start_Position.

    2) Gene mode (interval, weighted average):
       `df` index name is 'ensembl_gene_id' or 'gene', and columns
       {'Chromosome','start','end'} exist.
       → For each gene interval [start, end), compute a
         length-weighted average of feature(s) across overlapping
         bins.

    Parameters
    ----------
    df : pd.DataFrame
        Variant-like or gene-like input (see modes above).
    binned_df : pd.DataFrame
        BED-like fixed-bin table with columns:
        - 'Chromosome', 'region_start' (0-based), 'region_end'
          (exclusive)
        - plus one or more feature columns.
    feature_cols : str | list[str] | None, default None
        Which feature columns from `binned_df` to attach.
        If None, uses the all columns after the required three:
        ['Chromosome', 'region_start', 'region_end'].
    bin_size : int | None, default None
        Expected bin width. If None, inferred as the most common width
        among all **non-terminal** bins across chromosomes. The final
        bin on each chromosome may be shorter.

    Returns
    -------
    pd.DataFrame
        Copy of `df` with the requested feature columns added.
        Index is preserved (gene IDs stay as index in gene mode).

    Raises
    ------
    ValueError
        If required columns are missing, no feature columns are found,
        or bin sizes are inconsistent.

    """
    # ----- Validate binned_df -----
    required_bins = {"Chromosome", "region_start", "region_end"}
    if not required_bins.issubset(binned_df.columns):
        raise ValueError(f"binned_df must include {required_bins}")

    binned_df = binned_df.copy()
    binned_df = binned_df.sort_values(["Chromosome", "region_start"])

    # Pick default feature column(s) if None
    if feature_cols is None:
        candidate_cols = [c for c in binned_df.columns
                          if c not in ("Chromosome",
                                       "region_start",
                                       "region_end")]
        if not candidate_cols:
            raise ValueError("No feature columns found in binned_df.")
        feature_cols = candidate_cols
    elif isinstance(feature_cols, str):
        feature_cols = [feature_cols]

    # Ensure features are numeric (non-numeric → NaN)
    binned_df[feature_cols] = binned_df[feature_cols].apply(
        pd.to_numeric, errors="coerce")

    # ----- Infer/validate bin_size -----
    # Use all widths except the last bin per chromosome to infer size
    widths_all = []
    for _, grp in binned_df.groupby("Chromosome", sort=False):
        w = (grp["region_end"].values[:-1] - grp["region_start"].values[:-1])
        if w.size:
            widths_all.extend(w.tolist())

    if not widths_all:
        raise ValueError("Could not infer bin sizes (not enough bins).")

    inferred = pd.Series(widths_all).mode().iloc[0]
    if bin_size is None:
        bin_size = int(inferred)

    # Check consistency (non-terminal bins must equal bin_size)
    for chrom, grp in binned_df.groupby("Chromosome", sort=False):
        w = (grp["region_end"].values[:-1] - grp["region_start"].values[:-1])
        if w.size and not np.all(w == bin_size):
            raise ValueError(
                f"Inconsistent non-terminal bin sizes on {chrom}; "
                f"expected {bin_size}, saw {np.unique(w)}")

    # ----- Detect mode -----
    has_variant_cols = ({"Chromosome", "Start_Position"}.issubset(df.columns)
                        and (df.index.name == 'variant'))
    has_gene_cols = ({"Chromosome", "start", "end"}.issubset(df.columns)
                     and (df.index.name in ("ensembl_gene_id", "gene")))
    if not has_variant_cols and not has_gene_cols:
        raise ValueError(
            "Input `df` must be either:\n"
            " - Variant mode: columns {'Chromosome','Start_Position'}; or\n"
            " - Gene mode: index 'ensembl_gene_id' or 'gene' and "
            "columns {'Chromosome','start','end'}.")

    out = df.copy()

    # ---------- Variant mode ----------
    if has_variant_cols:
        tmp = out.copy()
        tmp["__order"] = np.arange(len(tmp))  # to preserve row order
        tmp = tmp.reset_index()

        pos = (pd.to_numeric(tmp["Start_Position"], errors="coerce")
               .fillna(-1)
               .astype(np.int64))
        region_start = (np.maximum(pos, 0) // bin_size) * bin_size
        tmp["region_start"] = region_start

        merged = tmp.merge(
            binned_df[["Chromosome", "region_start"] + feature_cols],
            on=["Chromosome", "region_start"],
            how="left")

        merged = (merged.sort_values("__order")
                  .drop(columns=["__order", "region_start"]))
        merged = merged.set_index("variant")
        return merged

    # ---------- Gene mode ----------
    idx_name = out.index.name  # preserve
    out = out.sort_values(["Chromosome", "start"]).copy()
    out["start"] = (pd.to_numeric(out["start"], errors="coerce")
                    .astype("Int64"))
    out["end"] = (pd.to_numeric(out["end"], errors="coerce")
                  .astype("Int64"))

    # Initialize output columns
    for c in feature_cols:
        out[c] = np.nan

    # Pre-sort bins and prepare arrays
    for chrom, bins in binned_df.groupby("Chromosome", sort=False):
        genes_chr = out[out["Chromosome"] == chrom]
        if genes_chr.empty:
            continue

        bs = bins["region_start"].to_numpy(np.int64)
        be = bins["region_end"].to_numpy(np.int64)
        feats = np.column_stack([
            pd.to_numeric(bins[c], errors="coerce").to_numpy(float)
            for c in feature_cols])  # shape: n_bins × n_feats

        for gid, row in genes_chr.iterrows():
            gs = row["start"]
            ge = row["end"]
            if pd.isna(gs) or pd.isna(ge) or int(ge) <= int(gs):
                continue
            gs = int(gs)
            ge = int(ge)

            # Overlapping bin indices [i0, i1)
            i0 = np.searchsorted(be, gs, side="right")
            i1 = np.searchsorted(bs, ge, side="left")
            if i0 >= i1:
                continue

            s = np.maximum(bs[i0:i1], gs)
            e = np.minimum(be[i0:i1], ge)
            w = (e - s).astype(np.int64)  # overlap lengths
            v = feats[i0:i1, :]           # values per overlapping bin
            ok = (w > 0)[:, None] & np.isfinite(v)

            w2 = np.where(ok, w[:, None], 0)
            num = (v * w2).sum(axis=0)
            den = w2.sum(axis=0)
            val = np.divide(num,
                            den,
                            out=np.full(len(feature_cols), np.nan),
                            where=den > 0)

            out.loc[gid, feature_cols] = val

    out = out.sort_index()
    out.index.name = idx_name
    return out


def annotate_min_dist_to_regions(df, regions_df, label=None):
    """Annotate a DataFrame with minimum distance to genomic regions.

    The DataFrame can by either a mutation or a variant DataFrame.

    For each mutation or variant in `df`, this function computes the
    shortest Euclidean distance to a set of genomic regions specified
    in a BED-like file.  If a mutation lies within any of the regions,
    a distance of 0 is assigned.  Otherwise, the distance to the
    nearest region boundary (start or end) is used.

    The returned DataFrame has a new column named 'cov_dist_<label>'
    indicating this distance. If `df` is indexed by 'variant', the
    function preserves that.

    Args:
        df (pd.DataFrame): DataFrame containing mutations or
            variants. Must have:
            - 'Chromosome': str, e.g., 'chr1'
            - 'Start_Position': int, genomic coordinate (1-based)
            - Index optionally named 'variant'

        regions_df : pd.DataFrame
            BED-like DataFrame as returned by `read_bed_file`.
            Must contain:
            - 'Chromosome'
            - 'region_start' (0-based)
            - 'region_end' (exclusive)

        label (str or None): Label to distinguish this region type in
            the output column. If None, the column is named
            'cov_dist_'. If 'term_site', the column becomes
            'cov_dist_term_site', etc.

    Returns:
        pd.DataFrame: A new DataFrame with an additional column:
            - 'cov_dist_<label>': float, minimum distance to any
              region, or 0 if overlapping, or NaN if chromosome not
              present.

            If input was indexed by 'variant', output preserves this.

    """
    new_df = df.copy()
    is_variant_df = (new_df.index.name == "variant")

    if is_variant_df:
        new_df = new_df.reset_index()

    distances = []

    # Group regions by chromosome for fast access
    term_by_chr = {
        chrom: group[["region_start", "region_end"]].values
        for chrom, group in regions_df.groupby("Chromosome")
    }

    for idx, row in df.iterrows():
        chrom = row["Chromosome"]
        # logic below still applies if pos is NaN, it will return NaN
        pos = row["Start_Position"]

        if chrom not in term_by_chr:
            # No data for this chromosome probably chrX or chrY
            distances.append(np.nan)
            continue

        intervals = term_by_chr[chrom]
        starts = intervals[:, 0]
        ends = intervals[:, 1]

        # Check if within any termination interval
        within_any = np.any((starts <= pos) & (pos <= ends))
        if within_any:
            distances.append(0)
        else:
            min_start_dist = np.min(np.abs(pos - starts))
            min_end_dist = np.min(np.abs(pos - ends))
            distances.append(min(min_start_dist, min_end_dist))

    # Add final column
    final_col_name = "cov_dist"
    if label is not None:
        final_col_name = final_col_name + f"_{label}"

    new_df[final_col_name] = distances
    return (new_df.set_index('variant')
            if is_variant_df
            else new_df)


def annotate_indicator_in_region(
        df: "pd.DataFrame",
        regions_df: "pd.DataFrame",
        label: str | None = None) -> "pd.DataFrame":
    """Add a 0 or 1 if mutation belongs or not to a region.

    Parameters
    ----------
    df : pandas.DataFrame
        Input table of mutations or aggregated variants.
        Must contain the columns

        - 'Chromosome'
        - 'Start_Position'

        If the index is named 'variant' it will be preserved.

    regions_df : pandas.DataFrame
        BED-like table returned by :func:`read_bed_file` with the
        columns

        - 'Chromosome'
        - 'region_start' (0-based, inclusive)
        - 'region_end'   (0-based, exclusive)

    label : str or None, optional
        Suffix for the output column name.
        - None: column name 'cov'
        - 'term': column name is f'cov_{term}'

    Returns
    -------
    pandas.DataFrame
        Copy of df with an added integer column (0 / 1 / NaN):

        - 1: coordinate falls inside **any** region on the same chr
        - 0: coordinate is outside all regions, including the case
          when the chromosome itself does not have any regions (except
          if it is chrX or ChrY)
        - NaN: chromosome absent from regions_df (happens for example
          for chrX and chrY in replication timing data) or
          'Start_Position' is NaN

        The 'variant' index is restored if it existed in *df*.

    """
    out = df.copy()
    has_variant_index = (out.index.name == "variant")
    if has_variant_index:
        out = out.reset_index()

    # Pre-group region intervals by chromosomed
    by_chr = {
        chrom: group[["region_start", "region_end"]].values
        for chrom, group in regions_df.groupby("Chromosome")}

    flags = []

    for _, row in out.iterrows():
        chrom = row["Chromosome"]
        pos = row["Start_Position"]

        # NaN if position missing
        if pd.isna(pos):
            flags.append(np.nan)
            continue

        # chromosome entirely absent from region file
        if chrom not in by_chr:
            flags.append(
                np.nan if chrom in {"chrX", "chrY", "X", "Y"}
                else 0)
            continue

        intervals = by_chr[chrom]
        starts, ends = intervals[:, 0], intervals[:, 1]

        inside = ((starts <= pos) & (pos < ends)).any()
        flags.append(1 if inside else 0)

    col_name = "cov" if label is None else f"cov_{label}"
    out[col_name] = flags

    return out.set_index("variant") if has_variant_index else out


def annotate_with_gene_features(
        df: pd.DataFrame,
        gene_feature_df: pd.DataFrame,
        feature_cols,
        gene_col: str = "ensembl_gene_id",
        strip_version: bool = True,
        labels: str | Sequence[str] | None = None) -> pd.DataFrame:
    """Add gene-level features to a mutation/variant DataFrame.

    Each row in *df* is matched to *gene_feature_df* on *gene_col*
    (optionally after removing the trailing “.v” version). The
    requested *feature_cols* are copied across and renamed according
    to *labels*:

    - ``labels=None`` (default): output columns are ``cov_{feature}``
      for each feature in *feature_cols*.
    - ``labels=str``: if one feature, output is ``cov_{labels}``; if
      multiple features, outputs are ``cov_{labels}_{feature}``.
    - ``labels=sequence[str]``: must match the length of *feature_cols*;
      outputs are ``cov_{label}`` for each provided label.

    The function also collapses duplicate versioned Ensembl IDs in
    *gene_feature_df* by keeping the highest numeric version.

    Parameters
    ----------
    df
        Table with at least a column named *gene_col*. If indexed by
        ``variant``, that index is preserved.
    gene_feature_df
        Gene-level feature matrix. Its index or one of its columns
        must hold *gene_col*. Version suffix ``.v`` is allowed. It can
        be a Series, if only one `feature_cols` is provided.
    feature_cols
        One or several column names in *gene_feature_df* to append.
    gene_col
        Join key name. Default ``ensembl_gene_id``.
    strip_version
        If True (default), remove the ``.v`` version from Ensembl IDs
        in both tables before matching. When collapsing duplicates in
        *gene_feature_df*, the highest version is kept.
    labels
        See description above for naming behavior.

    Returns
    -------
    pandas.DataFrame
        Copy of *df* with the requested features appended and renamed.

    Notes
    -----
    The merge is validated as many-to-one (many variants per gene, at
    most one feature row per gene). If this is violated, pandas raises.

    """
    # ── Normalize feature_cols to list ───────────────────────────────
    if isinstance(feature_cols, str):
        feature_cols = [feature_cols]

    if isinstance(gene_feature_df, pd.Series):
        gene_feature_df = pd.DataFrame(gene_feature_df)

    # ── Check requested columns exist ────────────────────────────────
    missing = [c for c in feature_cols
               if c not in gene_feature_df.columns]
    if missing:
        raise ValueError(f"Missing columns in gene_feature_df: {missing}")

    # ── Build output column names from labels spec ───────────────────
    if labels is None:
        out_names = [f"cov_{c}" for c in feature_cols]
    elif isinstance(labels, str):
        if len(feature_cols) == 1:
            out_names = [f"cov_{labels}"]
        else:
            out_names = [f"cov_{labels}_{c}" for c in feature_cols]
    else:
        # sequence of labels
        if len(labels) != len(feature_cols):
            raise ValueError("labels must have same length as feature_cols.")
        out_names = [f"cov_{lab}" for lab in labels]
    rename_map = dict(zip(feature_cols, out_names))

    # ── Prepare gene_feature_df for a clean many-to-one merge ────────
    gfeat = gene_feature_df.copy()

    # Ensure join key is an ordinary column
    if gene_col in gfeat.columns:
        pass
    elif gfeat.index.name:  # named index
        gfeat = (gfeat
                 .reset_index()
                 .rename(columns={gfeat.index.name: gene_col}))
    else:  # unnamed index
        gfeat = (gfeat
                 .reset_index()
                 .rename(columns={"index": gene_col}))

    if strip_version:
        key = gfeat[gene_col].astype(str)
        sp = key.str.split(".", n=1, expand=True)
        gfeat["_stable"] = sp[0]
        gfeat["_version"] = pd.to_numeric(sp[1], errors="coerce")

        # Keep highest version per stable ID
        gfeat = gfeat.sort_values("_version", na_position="first")
        gfeat = gfeat.drop_duplicates("_stable", keep="last")

        # Use stable ID as merge key
        gfeat = (gfeat
                 .drop(columns=[gene_col, "_version"])
                 .rename(columns={"_stable": gene_col}))

    gfeat = gfeat[[gene_col] + feature_cols]

    # ── Prepare df for merge while keeping index semantics ───────────
    out = df.copy()
    has_variant_index = (out.index.name == "variant")
    if has_variant_index:
        out = out.reset_index()

    if gene_col not in out.columns:
        raise KeyError(f"'df' must contain join column '{gene_col}'.")

    if strip_version:
        out[gene_col] = out[gene_col].astype(str).str.split(".", n=1).str[0]

    # ── Merge (many variants → one gene row), then rename features ───
    before = len(out)
    out = out.merge(gfeat, on=gene_col, how="left", validate="many_to_one")
    out = out.rename(columns=rename_map)

    # Warn about unmatched genes (use first appended column if present)
    first_new = out_names[0] if out_names else None
    if first_new and first_new in out.columns:
        unmatched = out[first_new].isna().sum()
        if unmatched:
            logger.warning(f"{unmatched} of {before} rows had no gene "
                           "feature match.")

    # ── Restore variant index if present ─────────────────────────────
    return out.set_index("variant") if has_variant_index else out
