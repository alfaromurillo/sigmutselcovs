"""Quality checks and preprocessing fixes for covariate matrices.

Typical workflow
----------------
1. Build a covariate matrix (genes × features) from whichever sources
   are relevant (GTEx, TCGA ATAC-seq, Roadmap, replication timing …).
2. Call :func:`check_all` to get a summary of potential issues.
3. Call :func:`fix_all` to apply recommended corrections and get a
   clean matrix ready for PCA or Riemannian dimensionality reduction.

Individual check/fix functions are available for finer control.
Missingness and collinearity are checked but not auto-fixed, as the
right action (imputation, dropping, etc.) is context-dependent.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Variance
# ---------------------------------------------------------------------------

def check_variance(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        cv_threshold: float = 0.01,
        dominant_pct_threshold: float = 0.95,
) -> pd.DataFrame:
    """Flag constant and near-constant columns.

    A column is flagged as **constant** when it has exactly one unique
    non-NaN value (variance zero).  It is flagged as **near-constant**
    when either:

    * the coefficient of variation ``std / |mean|`` falls below
      *cv_threshold* (evaluated only when ``|mean| > 0``), or
    * one value accounts for more than *dominant_pct_threshold* of all
      non-NaN entries.

    Constant columns cause division-by-zero during PCA standardisation.
    Near-constant columns add noise without contributing signal.

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns : list[str] or None
        Columns to inspect.  ``None`` checks all numeric columns.
    cv_threshold : float, default 0.01
        Coefficient-of-variation threshold for near-constant detection.
    dominant_pct_threshold : float, default 0.95
        Fraction of non-NaN entries that must share the most common
        value to flag a column as near-constant by dominance.

    Returns
    -------
    pd.DataFrame
        Index = column name.  Fields:

        ``variance``
            Sample variance (``NaN`` if fewer than two values).
        ``cv``
            ``std / |mean|``; ``NaN`` when mean is zero.
        ``pct_most_common``
            Fraction of non-NaN entries equal to the most common value.
        ``is_constant``
            True if the column has at most one unique non-NaN value.
        ``is_near_constant``
            True if near-constant by CV or dominance (but not constant).
        ``drop``
            True if the column should be removed (constant or near-constant).
    """
    cols = (
        columns if columns is not None
        else df.select_dtypes(include="number").columns.tolist()
    )

    records = []
    for col in cols:
        vals = df[col].dropna().to_numpy(dtype=float)
        n = len(vals)

        if n == 0:
            records.append(dict(
                column=col, variance=float("nan"), cv=float("nan"),
                pct_most_common=float("nan"),
                is_constant=False, is_near_constant=False, drop=False))
            continue

        var = float(np.var(vals, ddof=1)) if n > 1 else 0.0
        mean = float(np.mean(vals))
        std = float(np.std(vals, ddof=1)) if n > 1 else 0.0
        cv = std / abs(mean) if abs(mean) > 0 else float("nan")

        unique_vals, counts = np.unique(vals, return_counts=True)
        pct_most_common = float(counts.max() / n)

        is_constant = len(unique_vals) <= 1
        is_near_constant = (
            not is_constant
            and (
                (not np.isnan(cv) and cv < cv_threshold)
                or pct_most_common > dominant_pct_threshold
            )
        )

        records.append(dict(
            column=col, variance=var, cv=cv,
            pct_most_common=pct_most_common,
            is_constant=is_constant,
            is_near_constant=is_near_constant,
            drop=is_constant or is_near_constant))

    return pd.DataFrame(records).set_index("column")


def fix_variance(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        cv_threshold: float = 0.01,
        dominant_pct_threshold: float = 0.95,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop constant and near-constant columns.

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns, cv_threshold, dominant_pct_threshold
        Forwarded to :func:`check_variance`.

    Returns
    -------
    fixed_df : pd.DataFrame
        Copy of *df* with flagged columns removed.
    report : pd.DataFrame
        Full output of :func:`check_variance`.
    """
    report = check_variance(
        df, columns=columns,
        cv_threshold=cv_threshold,
        dominant_pct_threshold=dominant_pct_threshold)

    to_drop = report.index[report["drop"]].tolist()
    if to_drop:
        logger.info(
            "fix_variance: dropping %d column(s): %s", len(to_drop), to_drop)

    return df.drop(columns=to_drop), report


# ---------------------------------------------------------------------------
# Skewness
# ---------------------------------------------------------------------------

def check_skewness(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        skew_threshold: float = 1.0,
        max_median_ratio_threshold: float = 10.0,
) -> pd.DataFrame:
    """Report skewness statistics and log-transform recommendations per column.

    A log transform is recommended when a column is non-negative, has at
    least one positive value, and is right-skewed by either criterion:
    skewness coefficient > *skew_threshold*, or
    ``max / median`` > *max_median_ratio_threshold*.

    The pseudo-count ``c`` in ``log(c + x)`` follows the half-minimum
    rule: ``c = min_positive / 2`` when zeros are present, keeping zeros
    one unit below the detection floor in log space.  When there are no
    zeros, ``c = 0`` (plain ``log(x)``).

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns : list[str] or None
        Columns to inspect.  ``None`` checks all numeric columns.
    skew_threshold : float, default 1.0
        Fisher skewness above which a log transform is recommended.
    max_median_ratio_threshold : float, default 10.0
        ``max / median`` above which a log transform is recommended
        regardless of the skewness coefficient.

    Returns
    -------
    pd.DataFrame
        Index = column name.  Fields:

        ``skewness``
            Fisher skewness coefficient.
        ``max_median_ratio``
            ``max(x) / median(x)``; ``inf`` when median is zero.
        ``all_nonneg``
            True if no value is strictly negative.
        ``apply_log``
            True if a log transform is recommended.
        ``pseudo_count``
            ``c`` in ``log(c + x)``.  Zero means plain ``log(x)``
            (no zeros present).
    """
    from scipy.stats import skew as _skew

    cols = (
        columns if columns is not None
        else df.select_dtypes(include="number").columns.tolist()
    )

    records = []
    for col in cols:
        vals = df[col].dropna().to_numpy(dtype=float)

        if len(vals) == 0:
            records.append(dict(
                column=col, skewness=float("nan"),
                max_median_ratio=float("nan"),
                all_nonneg=False, apply_log=False, pseudo_count=0.0))
            continue

        all_nonneg = bool((vals >= 0).all())
        has_positives = bool((vals > 0).any())
        s = float(_skew(vals))
        median = float(np.median(vals))
        mmr = float(vals.max() / median) if median > 0 else float("inf")

        apply = (
            all_nonneg
            and has_positives
            and (s > skew_threshold or mmr > max_median_ratio_threshold)
        )

        pseudo_count = (
            float(vals[vals > 0].min()) / 2.0
            if apply and (vals == 0).any()
            else 0.0
        )

        records.append(dict(
            column=col, skewness=s, max_median_ratio=mmr,
            all_nonneg=all_nonneg, apply_log=apply,
            pseudo_count=pseudo_count))

    return pd.DataFrame(records).set_index("column")


def fix_skewness(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        skew_threshold: float = 1.0,
        max_median_ratio_threshold: float = 10.0,
        iterative: bool = True,
        max_iter: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply log transforms to right-skewed non-negative columns.

    Calls :func:`check_skewness` and applies ``log(c + x)`` to each
    flagged column, where ``c`` is determined by the half-minimum rule.

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns : list[str] or None
        Columns to inspect.  ``None`` checks all numeric columns.
    skew_threshold, max_median_ratio_threshold
        Forwarded to :func:`check_skewness`.
    iterative : bool, default True
        If True, repeat the check-and-transform cycle until no column
        is flagged or *max_iter* passes are exhausted.  Handles rare
        cases where a single log pass leaves residual skew (e.g.\
        data spanning many orders of magnitude).
    max_iter : int, default 5
        Maximum number of passes when *iterative* is True.

    Returns
    -------
    fixed_df : pd.DataFrame
        Copy of *df* with recommended transforms applied.
        Non-flagged and non-numeric columns are unchanged.
    report : pd.DataFrame
        Output of :func:`check_skewness` from the final pass.
    """
    import warnings

    fixed = df.copy()
    report = None
    for i in range(max_iter if iterative else 1):
        report = check_skewness(
            fixed, columns=columns,
            skew_threshold=skew_threshold,
            max_median_ratio_threshold=max_median_ratio_threshold)

        n_flagged = int(report["apply_log"].sum())
        if n_flagged == 0:
            logger.info("fix_skewness: converged after %d pass(es)", i)
            break

        for col, row in report.iterrows():
            if row["apply_log"]:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    fixed[col] = np.log(row["pseudo_count"] + fixed[col])

        logger.info(
            "fix_skewness pass %d: transformed %d / %d columns",
            i + 1, n_flagged, len(report))
    else:
        logger.warning(
            "fix_skewness: reached max_iter=%d without full convergence",
            max_iter)

    return fixed, report


# ---------------------------------------------------------------------------
# Missingness
# ---------------------------------------------------------------------------

def check_missingness(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        missing_threshold: float = 0.5,
) -> pd.DataFrame:
    """Report the fraction of missing values per column.

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns : list[str] or None
        Columns to inspect.  ``None`` checks all columns.
    missing_threshold : float, default 0.5
        Columns with a higher missing fraction are flagged.

    Returns
    -------
    pd.DataFrame
        Index = column name.  Fields:

        ``n_missing``
            Count of NaN entries.
        ``pct_missing``
            Fraction of entries that are NaN.
        ``flag``
            True if ``pct_missing > missing_threshold``.
    """
    cols = columns if columns is not None else df.columns.tolist()
    n_rows = len(df)

    records = []
    for col in cols:
        n_missing = int(df[col].isna().sum())
        pct = n_missing / n_rows if n_rows > 0 else float("nan")
        records.append(dict(
            column=col, n_missing=n_missing, pct_missing=pct,
            flag=pct > missing_threshold))

    return pd.DataFrame(records).set_index("column")


# ---------------------------------------------------------------------------
# Collinearity
# ---------------------------------------------------------------------------

def check_collinearity(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        corr_threshold: float = 0.95,
) -> pd.DataFrame:
    """Flag pairs of columns with near-perfect Pearson correlation.

    High correlation indicates redundant features.  PCA absorbs
    redundancy, but non-linear dimensionality reduction (e.g. Riemannian
    stats) can be dominated by groups of near-identical columns.

    Note: for large numbers of ATAC-seq samples from the same tissue,
    many pairs may be flagged — this is expected and informative rather
    than actionable in all cases.

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns : list[str] or None
        Columns to include.  ``None`` uses all numeric columns.
    corr_threshold : float, default 0.95
        Absolute Pearson r above which a pair is flagged.

    Returns
    -------
    pd.DataFrame
        One row per flagged pair.  Columns:

        ``col_a``, ``col_b``
            The two column names.
        ``correlation``
            Pearson r between the pair.

        Empty DataFrame if no pairs exceed the threshold.
    """
    cols = (
        columns if columns is not None
        else df.select_dtypes(include="number").columns.tolist()
    )

    corr = df[cols].corr(method="pearson")
    cols_list = list(cols)
    records = []
    for i, a in enumerate(cols_list):
        for b in cols_list[i + 1:]:
            r = float(corr.loc[a, b])
            if abs(r) > corr_threshold:
                records.append(dict(col_a=a, col_b=b, correlation=r))

    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=["col_a", "col_b", "correlation"])


# ---------------------------------------------------------------------------
# Combined interface
# ---------------------------------------------------------------------------

def check_all(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        cv_threshold: float = 0.01,
        dominant_pct_threshold: float = 0.95,
        skew_threshold: float = 1.0,
        max_median_ratio_threshold: float = 10.0,
        missing_threshold: float = 0.5,
        corr_threshold: float = 0.95,
) -> dict[str, pd.DataFrame]:
    """Run all covariate quality checks and return a report dictionary.

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns : list[str] or None
        Columns to inspect.  ``None`` checks all numeric columns.
    cv_threshold, dominant_pct_threshold
        Forwarded to :func:`check_variance`.
    skew_threshold, max_median_ratio_threshold
        Forwarded to :func:`check_skewness`.
    missing_threshold
        Forwarded to :func:`check_missingness`.
    corr_threshold
        Forwarded to :func:`check_collinearity`.

    Returns
    -------
    dict[str, pd.DataFrame]
        Keys: ``"variance"``, ``"skewness"``, ``"missingness"``,
        ``"collinearity"``.
    """
    return {
        "variance": check_variance(
            df, columns=columns,
            cv_threshold=cv_threshold,
            dominant_pct_threshold=dominant_pct_threshold),
        "skewness": check_skewness(
            df, columns=columns,
            skew_threshold=skew_threshold,
            max_median_ratio_threshold=max_median_ratio_threshold),
        "missingness": check_missingness(
            df, columns=columns,
            missing_threshold=missing_threshold),
        "collinearity": check_collinearity(
            df, columns=columns,
            corr_threshold=corr_threshold),
    }


def fix_all(
        df: pd.DataFrame,
        columns: list[str] | None = None,
        cv_threshold: float = 0.01,
        dominant_pct_threshold: float = 0.95,
        skew_threshold: float = 1.0,
        max_median_ratio_threshold: float = 10.0,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Drop low-variance columns then apply log transforms.

    Applies fixes in order:

    1. :func:`fix_variance` — drop constant / near-constant columns.
    2. :func:`fix_skewness` — log-transform right-skewed columns.

    Missingness and collinearity are reported by :func:`check_all` but
    not automatically fixed, as the right action is context-dependent.

    Parameters
    ----------
    df : pd.DataFrame
        Covariate matrix (genes × features).
    columns : list[str] or None
        Columns to inspect.  ``None`` processes all numeric columns.
    cv_threshold, dominant_pct_threshold
        Forwarded to :func:`fix_variance`.
    skew_threshold, max_median_ratio_threshold
        Forwarded to :func:`fix_skewness`.

    Returns
    -------
    fixed_df : pd.DataFrame
        Cleaned covariate matrix.
    reports : dict[str, pd.DataFrame]
        ``"variance"`` and ``"skewness"`` reports from each fix step.
    """
    fixed, var_report = fix_variance(
        df, columns=columns,
        cv_threshold=cv_threshold,
        dominant_pct_threshold=dominant_pct_threshold)

    fixed, skew_report = fix_skewness(
        fixed,
        skew_threshold=skew_threshold,
        max_median_ratio_threshold=max_median_ratio_threshold)

    return fixed, {"variance": var_report, "skewness": skew_report}
