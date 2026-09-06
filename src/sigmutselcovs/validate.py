"""Sanity and plausibility checks for built covariate matrices.

`validate_covariates` runs two tiers of checks on the **raw**
(pre-`fix_all`) full covariate matrix:

- **Data sanity**: NaN fractions, constant columns, value ranges,
  expected column counts against the registry, index integrity.
  These catch truncated downloads, empty joins, and assembly
  mismatches.
- **Biological plausibility**: direction-only Spearman correlations
  that any real tissue should show (expression anticorrelates with
  late replication, active marks and open chromatin correlate with
  expression, repressive marks anticorrelate, replicate epigenomes
  agree), plus a housekeeping-gene panel that must be expressed,
  early-replicating and accessible.

Each check yields a row (check, scope, status, value, threshold,
note) with status ``pass``, ``warn`` or ``fail``.  Direction misses
are warnings (a tissue can be legitimately unusual); a batch of them
is a failure.
"""

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .covariates_utilities import sanitize_feature_label
from .registry import ProjectSpec, get_project

logger = logging.getLogger(__name__)

# Housekeeping panel (Ensembl gene IDs): GAPDH, ACTB, RPL13A, B2M.
HOUSEKEEPING_GENES = {
    "GAPDH": "ENSG00000111640",
    "ACTB": "ENSG00000075624",
    "RPL13A": "ENSG00000142541",
    "B2M": "ENSG00000166710",
}

_ACTIVE_MARKS = (
    "h3k4me3",
    "h3k27ac",
    "h3k9ac",
    # ENCODE-sourced marks (encode_chromatin) beyond Roadmap's
    # default 7 -- see mutation_rates' covariates.md research on
    # dNdScv's 10 original marks. All three are activation-associated
    # acetylation marks, same expected direction as h3k27ac/h3k9ac.
    "h3k23ac",
    "h3k14ac",
    "h2ak9ac",
    # DNase-seq (open chromatin): same "_fc_signal_" naming
    # convention as the marks above (see download.py), and the same
    # expected direction as ATAC's open-chromatin check.
    "dnase",
)
_REPRESSIVE_MARKS = ("h3k9me3", "h3k27me3")
_DIRECTION_RHO = 0.05  # |rho| below this counts as "no clear signal"


class _Report:
    def __init__(self):
        self.rows: list[dict] = []

    def add(
        self,
        check: str,
        scope: str,
        status: str,
        value=None,
        threshold=None,
        note: str = "",
    ) -> None:
        self.rows.append(
            {
                "check": check,
                "scope": scope,
                "status": status,
                "value": value,
                "threshold": threshold,
                "note": note,
            }
        )

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            self.rows,
            columns=[
                "check",
                "scope",
                "status",
                "value",
                "threshold",
                "note",
            ],
        )


def _spearman(df: pd.DataFrame, a: pd.Series, b: pd.Series) -> float:
    joined = pd.concat([a, b], axis=1, keys=["a", "b"]).dropna()
    if len(joined) < 10:
        return float("nan")
    return float(joined["a"].corr(joined["b"], method="spearman"))


def _expression_series(
    df: pd.DataFrame, spec: ProjectSpec
) -> pd.Series | None:
    for col in (spec.gtex.representative_column, "tpm_unstranded"):
        if col in df.columns:
            return df[col]
    return None


def _mark_mean(
    df: pd.DataFrame, mark: str, where: str
) -> pd.Series | None:
    cols = [
        c
        for c in df.columns
        if mark in c and "fc_signal" in c and c.endswith(where)
    ]
    if not cols:
        return None
    return df[cols].mean(axis=1)


def _check_direction(
    report: _Report,
    df: pd.DataFrame,
    check: str,
    a: pd.Series | None,
    b: pd.Series | None,
    expected_sign: int,
    scope: str,
) -> None:
    if a is None or b is None:
        report.add(
            check,
            scope,
            "pass",
            note="not applicable " "(columns absent)",
        )
        return
    rho = _spearman(df, a, b)
    if np.isnan(rho):
        report.add(
            check,
            scope,
            "warn",
            value=rho,
            note="too few complete pairs",
        )
        return
    signed = rho * expected_sign
    if signed >= _DIRECTION_RHO:
        report.add(check, scope, "pass", value=round(rho, 3))
    elif signed > -_DIRECTION_RHO:
        report.add(
            check,
            scope,
            "warn",
            value=round(rho, 3),
            threshold=_DIRECTION_RHO * expected_sign,
            note="no clear signal",
        )
    else:
        report.add(
            check,
            scope,
            "warn",
            value=round(rho, 3),
            threshold=_DIRECTION_RHO * expected_sign,
            note="direction opposite to expectation",
        )


def _sanity_tier(
    report: _Report,
    df: pd.DataFrame,
    spec: ProjectSpec,
    *,
    nan_warn: float,
    nan_fail: float,
) -> None:
    # NaN fractions
    nan_fraction = df.isna().mean()
    worst = float(nan_fraction.max()) if len(nan_fraction) else 0.0
    n_fail = int((nan_fraction > nan_fail).sum())
    n_warn = int((nan_fraction > nan_warn).sum()) - n_fail
    if n_fail:
        report.add(
            "nan_fraction",
            "matrix",
            "fail",
            value=worst,
            threshold=nan_fail,
            note=f"{n_fail} columns above {nan_fail:.0%} NaN: "
            + ", ".join(
                nan_fraction[nan_fraction > nan_fail].index[:5]
            ),
        )
    elif n_warn:
        report.add(
            "nan_fraction",
            "matrix",
            "warn",
            value=worst,
            threshold=nan_warn,
            note=f"{n_warn} columns above {nan_warn:.0%} NaN",
        )
    else:
        report.add(
            "nan_fraction",
            "matrix",
            "pass",
            value=worst,
            threshold=nan_warn,
        )

    # constant columns (fix_variance would drop them silently)
    numeric = df.select_dtypes("number")
    constant = [
        c
        for c in numeric.columns
        if numeric[c].nunique(dropna=True) <= 1
    ]
    report.add(
        "constant_columns",
        "matrix",
        "warn" if constant else "pass",
        value=len(constant),
        note=", ".join(constant[:5]),
    )

    # value ranges
    if "mrt" in df.columns:
        mrt = df["mrt"].dropna()
        ok = bool(((mrt >= 0) & (mrt <= 1)).all())
        report.add(
            "mrt_range",
            "repliseq",
            "pass" if ok else "fail",
            value=[
                round(float(mrt.min()), 3),
                round(float(mrt.max()), 3),
            ],
            threshold="[0, 1]",
        )
    for label, prefix in (
        ("gtex_nonnegative", "gtex_"),
        ("roadmap_nonnegative", "fc_signal"),
    ):
        cols = [
            c
            for c in df.columns
            if (
                c.startswith(prefix)
                if prefix.endswith("_")
                else prefix in c
            )
        ]
        if not cols:
            continue
        minimum = float(df[cols].min().min())
        report.add(
            label,
            "matrix",
            "pass" if minimum >= 0 else "fail",
            value=minimum,
            threshold=">= 0",
        )
    if "tpm_unstranded" in df.columns:
        minimum = float(df["tpm_unstranded"].min())
        report.add(
            "tpm_nonnegative",
            "gexp",
            "pass" if minimum >= 0 else "fail",
            value=minimum,
            threshold=">= 0",
        )

    # column counts vs registry
    if spec.repliseq is not None:
        clr_cols = [c for c in df.columns if c.startswith("clr_rt_s")]
        n = spec.repliseq.n_fractions
        if clr_cols and n is not None:
            report.add(
                "clr_fraction_count",
                "repliseq",
                "pass" if len(clr_cols) == n else "fail",
                value=len(clr_cols),
                threshold=n,
            )
    if spec.atac is not None:
        atac_cols = [
            c
            for c in df.columns
            if c.startswith(spec.atac.column_prefix.lower() + "_")
        ]
        if atac_cols:
            report.add(
                "atac_body_promoter_pairs",
                "atac",
                "pass" if len(atac_cols) % 2 == 0 else "fail",
                value=len(atac_cols),
                threshold="even",
                note="each sample contributes body+promoter",
            )
    if spec.roadmap is not None:
        roadmap_cols = [c for c in df.columns if "fc_signal" in c]
        maximum = 2 * len(spec.roadmap.eids) * len(spec.roadmap.marks)
        if roadmap_cols:
            report.add(
                "roadmap_column_count",
                "roadmap",
                "pass" if len(roadmap_cols) <= maximum else "fail",
                value=len(roadmap_cols),
                threshold=f"<= {maximum}",
            )
    if spec.encode_chromatin is not None:
        encode_cols = [
            c
            for c in df.columns
            if "fc_signal" in c
            and any(
                sanitize_feature_label(t.accession) in c
                for t in spec.encode_chromatin.tracks
            )
        ]
        maximum = 2 * len(spec.encode_chromatin.tracks)
        if encode_cols:
            report.add(
                "encode_chromatin_column_count",
                "encode_chromatin",
                "pass" if len(encode_cols) <= maximum else "fail",
                value=len(encode_cols),
                threshold=f"<= {maximum}",
            )
    if spec.average_by_assay:
        # roadmap_column_count/encode_chromatin_column_count above
        # never fire here (collapsed columns no longer carry
        # "fc_signal"/accession substrings) -- this is their
        # equivalent for the collapsed shape: at most 2 columns
        # (body+promoter) per distinct mark pooled across sources.
        marks = set()
        if spec.roadmap is not None:
            marks.update(m.lower() for m in spec.roadmap.marks)
        if spec.encode_chromatin is not None:
            marks.update(
                t.label.lower() for t in spec.encode_chromatin.tracks
            )
        if spec.atac is not None:
            marks.add("atac")
        collapsed_cols = [
            c
            for c in df.columns
            if c.endswith(("_body", "_promoter"))
        ]
        maximum = 2 * len(marks)
        if collapsed_cols and maximum:
            report.add(
                "chromatin_collapsed_column_count",
                "chromatin_collapsed",
                "pass" if len(collapsed_cols) <= maximum else "fail",
                value=len(collapsed_cols),
                threshold=f"<= {maximum}",
            )

    # index integrity
    duplicated = int(df.index.duplicated().sum())
    ensg = df.index.astype(str).str.match(r"^ENSG\d+$")
    report.add(
        "index_unique",
        "matrix",
        "pass" if duplicated == 0 else "fail",
        value=duplicated,
        threshold=0,
    )
    bad = int((~ensg).sum())
    report.add(
        "index_ensg_format",
        "matrix",
        "pass" if bad == 0 else "warn",
        value=bad,
        threshold=0,
        note="non-ENSG or versioned ids in the index" if bad else "",
    )


def _biology_tier(
    report: _Report, df: pd.DataFrame, spec: ProjectSpec
) -> None:
    expression = _expression_series(df, spec)
    mrt = df["mrt"] if "mrt" in df.columns else None

    _check_direction(
        report,
        df,
        "expression_vs_mrt",
        expression,
        mrt,
        -1,
        "cross-block",
    )
    for mark in _ACTIVE_MARKS:
        _check_direction(
            report,
            df,
            f"{mark}_vs_expression",
            _mark_mean(df, mark, "body"),
            expression,
            +1,
            "roadmap",
        )
    for mark in _REPRESSIVE_MARKS:
        _check_direction(
            report,
            df,
            f"{mark}_vs_expression",
            _mark_mean(df, mark, "body"),
            expression,
            -1,
            "roadmap",
        )
    if spec.atac is not None:
        prefix = spec.atac.column_prefix.lower() + "_"
        atac_body = [
            c
            for c in df.columns
            if c.startswith(prefix) and c.endswith("body")
        ]
        atac = df[atac_body].mean(axis=1) if atac_body else None
        _check_direction(
            report,
            df,
            "atac_vs_expression",
            atac,
            expression,
            +1,
            "atac",
        )

    # replicate coherence: the same mark across Roadmap epigenomes
    if spec.roadmap is not None and len(spec.roadmap.eids) >= 2:
        for mark in set(_ACTIVE_MARKS) | set(_REPRESSIVE_MARKS):
            cols = [
                c
                for c in df.columns
                if mark in c and c.endswith("body")
            ]
            eids = {
                re.match(r"(e\d+)_", c).group(1)
                for c in cols
                if re.match(r"(e\d+)_", c)
            }
            if len(eids) < 2:
                continue
            pair = sorted(cols)[:2]
            _check_direction(
                report,
                df,
                f"{mark}_eid_coherence",
                df[pair[0]],
                df[pair[1]],
                +1,
                "roadmap",
            )

    # housekeeping panel
    present = {
        name: gid
        for name, gid in HOUSEKEEPING_GENES.items()
        if gid in df.index
    }
    if not present:
        report.add(
            "housekeeping_panel",
            "matrix",
            "warn",
            note="none of the panel genes in the index",
        )
    elif df.index.has_duplicates:
        report.add(
            "housekeeping_panel",
            "matrix",
            "warn",
            note="index has duplicates; panel skipped",
        )
    elif expression is not None:
        misses = []
        expression_rank = expression.rank(pct=True)
        mrt_rank = mrt.rank(pct=True) if mrt is not None else None
        for name, gid in present.items():
            rank = expression_rank.get(gid, np.nan)
            if not rank > 0.75:
                misses.append(f"{name} expression pct {rank:.2f}")
            if mrt_rank is not None:
                rank = mrt_rank.get(gid, np.nan)
                if not rank < 0.5:
                    misses.append(f"{name} mrt pct {rank:.2f}")
        report.add(
            "housekeeping_panel",
            "matrix",
            "warn" if misses else "pass",
            value=len(misses),
            note="; ".join(misses[:6]),
        )


def validate_covariates(
    project: str,
    data_dir: str | Path | None = None,
    *,
    cov_matrix_raw: pd.DataFrame | None = None,
    registry_path: str | Path | None = None,
    nan_warn: float = 0.5,
    nan_fail: float = 0.95,
    strict: bool = False,
    gencode_gtfs=None,
) -> pd.DataFrame:
    """Validate a project's covariate matrix.

    Parameters
    ----------
    project : str
        Registered TCGA study code.
    data_dir : str | Path | None
        Project data directory; used to (re)build the raw matrix
        from caches when ``cov_matrix_raw`` is not given.
    cov_matrix_raw : pd.DataFrame | None
        The raw (pre-fix_all) full matrix, if already in hand.
    strict : bool
        Raise ValueError when any check fails.

    Returns
    -------
    pd.DataFrame
        One row per check: check, scope, status (pass|warn|fail),
        value, threshold, note.
    """
    spec = get_project(project, registry_path)
    if cov_matrix_raw is None:
        if data_dir is None:
            raise ValueError("Give data_dir or cov_matrix_raw")
        from .builder import build_covariate_matrix

        cov_matrix_raw = build_covariate_matrix(
            project,
            data_dir,
            registry_path=registry_path,
            gencode_gtfs=gencode_gtfs,
            apply_fixes=False,
        ).full

    report = _Report()
    if cov_matrix_raw.columns.duplicated().any():
        duplicated = (
            cov_matrix_raw.columns[
                cov_matrix_raw.columns.duplicated()
            ]
            .unique()
            .tolist()
        )
        report.add(
            "duplicate_columns",
            "matrix",
            "warn",
            value=len(duplicated),
            note=f"checking first occurrence only: "
            f"{duplicated[:5]}",
        )
        cov_matrix_raw = cov_matrix_raw.loc[
            :, ~cov_matrix_raw.columns.duplicated()
        ]
    _sanity_tier(
        report,
        cov_matrix_raw,
        spec,
        nan_warn=nan_warn,
        nan_fail=nan_fail,
    )
    _biology_tier(report, cov_matrix_raw, spec)

    # a batch of direction misses is a failure even if each alone
    # is only a warning
    frame = report.frame()
    direction_misses = int(
        (
            (frame["status"] == "warn")
            & frame["note"].str.contains("opposite", na=False)
        ).sum()
    )
    if direction_misses >= 3:
        extra = pd.DataFrame(
            [
                {
                    "check": "biology_batch",
                    "scope": "matrix",
                    "status": "fail",
                    "value": direction_misses,
                    "threshold": "< 3",
                    "note": "multiple biological directions inverted — "
                    "check assemblies and column mapping",
                }
            ]
        )
        frame = pd.concat([frame, extra], ignore_index=True)

    for _, row in frame.iterrows():
        level = {
            "pass": logging.INFO,
            "warn": logging.WARNING,
            "fail": logging.ERROR,
        }[row["status"]]
        logger.log(
            level,
            "validate %s [%s]: %s value=%s %s",
            spec.code,
            row["status"],
            row["check"],
            row["value"],
            row["note"],
        )

    if strict and (frame["status"] == "fail").any():
        failures = frame[frame["status"] == "fail"]["check"].tolist()
        raise ValueError(f"Covariate validation failed: {failures}")
    return frame


def print_validation_report(frame: pd.DataFrame) -> None:
    """Print the validation report in a compact fixed-width form."""
    for _, row in frame.iterrows():
        symbol = {"pass": "ok", "warn": "WARN", "fail": "FAIL"}[
            row["status"]
        ]
        note = f"  {row['note']}" if row["note"] else ""
        print(
            f"[{symbol:4s}] {row['check']:28s} "
            f"value={row['value']}{note}"
        )
