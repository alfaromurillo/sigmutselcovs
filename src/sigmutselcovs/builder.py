"""Generalized covariate-matrix builder.

`build_covariate_matrix` assembles the three covariate matrices
(full, simple, tcga) for any project registered in
``data/projects.json``, from data already present in the project's
data directory (see `sigmutselcovs.download.download_covariates`).

The assembly reproduces the historical COAD pipeline exactly (same
loaders, same caches, same concatenation order), so building COAD
from an existing ``coad_data/`` tree is bit-for-bit identical to the
old ``coad_analysis/code/covariates.py`` ``build()``.
"""

import hashlib
import json
import logging
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import NamedTuple

import pandas as pd

from .covariates_checks import fix_all
from .covariates_chromatin import load_or_generate_chromatin_covariates
from .covariates_gene_expression import (
    import_gtex,
    load_or_generate_mean_tcga_gexp,
    load_or_generate_tcga_gexp_per_sample,
)
from .covariates_locations import location_gtex_tcga_mapping
from .covariates_replication_timing import (
    load_or_generate_mrt,
    load_or_generate_rt_fractions,
    load_or_generate_rt_wavelet,
)
from .covariates_utilities import clr_transform
from .paths import ProjectPaths, project_paths
from .registry import ProjectSpec, get_project, location_projects_registry

logger = logging.getLogger(__name__)

_SOURCES = ("gtex", "gexp", "repliseq", "roadmap", "atac")


class CovariateMatrices(NamedTuple):
    """The three covariate matrices; unpacks as the historical 3-tuple."""

    full: pd.DataFrame
    simple: pd.DataFrame
    tcga: pd.DataFrame


def _default_gencode_gtfs() -> dict[str, Path]:
    """GENCODE GTFs from sigmutsel (imported lazily: pulls in pymc)."""
    from sigmutsel.locations import (
        location_gencode19_annotation,
        location_gencode38_annotation,
    )
    return {"hg19": location_gencode19_annotation,
            "hg38": location_gencode38_annotation}


def _selected(include: Sequence[str] | None,
              exclude: Sequence[str]) -> set[str]:
    selected = set(include) if include is not None else set(_SOURCES)
    selected -= set(exclude)
    unknown = selected - set(_SOURCES)
    if unknown:
        raise ValueError(f"Unknown sources: {sorted(unknown)}; "
                         f"valid: {_SOURCES}")
    return selected


def _skip(source: str, reason: str) -> None:
    logger.warning("Skipping %s covariates: %s", source, reason)


def _load_repliseq(spec, paths: ProjectPaths, gtf: str | Path,
                   force_generation: bool) -> list[pd.DataFrame]:
    """Return the replication-timing frames for the full matrix.

    For fraction-based sources: [mrt, clr fractions] (CLR only when
    the fraction count supports it).  Wavelet sources yield a single
    rt_wavelet frame instead.
    """
    if spec.type == "mat":
        source = paths.rt_dir / spec.filename
        source_type = "mat"
        if not source.exists() and not paths.mrt_csv.exists():
            _skip("repliseq", f"neither {source} nor cache "
                  f"{paths.mrt_csv} exists")
            return []
    elif spec.type in ("fraction_bigwigs", "wavelet"):
        tracks = [paths.rt_encode_dir / f"{t.accession}.bigWig"
                  for t in spec.tracks]
        missing = [p for p in tracks if not p.exists()]
        cache = (paths.rt_wavelet_csv if spec.type == "wavelet"
                 else paths.mrt_csv)
        if missing and not cache.exists():
            _skip("repliseq", f"missing bigWigs {missing} and no "
                  f"cache at {cache}")
            return []
        if spec.type == "wavelet":
            rt_wavelet = load_or_generate_rt_wavelet(
                paths.rt_wavelet_csv,
                tracks[0],
                gtf,
                bin_size=spec.bin_size,
                force_generation=force_generation)
            return [rt_wavelet.to_frame()]
        source = tracks
        source_type = "fraction_bigwigs"
    else:
        raise ValueError(f"Unknown repliseq type {spec.type!r}")

    mrt_per_gene = load_or_generate_mrt(
        paths.mrt_csv,
        source,
        gtf,
        source_type=source_type,
        bin_size=spec.bin_size,
        mrt_fraction_cols=(list(spec.mrt_fraction_cols)
                           if spec.mrt_fraction_cols is not None
                           else None),
        force_generation=force_generation)

    rt_fractions = load_or_generate_rt_fractions(
        paths.rt_fractions_csv,
        source,
        gtf,
        source_type=source_type,
        bin_size=spec.bin_size,
        force_generation=force_generation)

    frames = [mrt_per_gene.to_frame()]
    n = rt_fractions.shape[1]
    use_clr = (spec.include_clr_fractions
               if spec.include_clr_fractions is not None
               else n >= 3)
    if use_clr:
        frames.append(clr_transform(rt_fractions).add_prefix("clr_"))
    else:
        logger.info(
            "Skipping CLR fractions: %d fractions give %s "
            "perfectly anticorrelated CLR columns; keeping mrt only",
            n, n)
    return frames


def _load_chromatin(source: str, covs_csv: Path,
                    bigwigs: list[Path], gtf: str | Path,
                    force_generation: bool) -> pd.DataFrame | None:
    """Chromatin covariates for one source (roadmap or TCGA ATAC)."""
    if not bigwigs and not covs_csv.exists():
        _skip(source, f"no bigWig files and no cache at {covs_csv}")
        return None
    if not bigwigs and force_generation:
        _skip(source, "force_generation without bigWig files")
        return None
    return load_or_generate_chromatin_covariates(
        covs_csv,
        bigwigs,
        gtf,
        biotypes=("protein_coding",),
        include_promoter=True,
        promoter_upstream=2000,
        promoter_downstream=200,
        force_generation=force_generation,
        average_by_assay=False)


def _atac_prefix(spec, atac_covs: pd.DataFrame) -> str:
    """Column prefix for the ATAC block of cov_matrix_tcga."""
    prefix = spec.column_prefix.lower()
    if any(c.startswith(prefix + "_") for c in atac_covs.columns):
        return prefix
    stems = [c.split("_")[0] for c in atac_covs.columns]
    fallback = stems[0] if stems and all(
        s == stems[0] for s in stems) else prefix
    logger.warning(
        "No ATAC column starts with %r; falling back to prefix %r "
        "derived from the bigWig names", prefix + "_", fallback)
    return fallback


def _describe_block(source: str, spec: ProjectSpec,
                    frame: pd.DataFrame) -> dict:
    """Block-level constants for the column dictionary."""
    if source == "gtex":
        return {"description": "GTEx median gene expression",
                "cell_line_or_epigenome": "",
                "assembly": "", "units": "median TPM"}
    if source == "gexp_mean":
        return {"description": "Mean TCGA gene expression over "
                "samples of one tissue type",
                "cell_line_or_epigenome": "",
                "assembly": "", "units": "TPM"}
    if source == "gexp_per_sample":
        return {"description": "Per-sample TCGA STAR expression "
                "metric ({barcode}_{metric})",
                "cell_line_or_epigenome": "",
                "assembly": "", "units": "STAR metric"}
    if source in ("mrt", "clr", "wavelet"):
        rt = spec.repliseq
        desc = {
            "mrt": "Mean replication time over gene body "
                   "(0 = early, 1 = late)",
            "clr": "CLR-transformed S-phase fraction over gene body",
            "wavelet": "Wavelet-smoothed early/late replication "
                       "signal over gene body",
        }[source]
        return {"description": desc,
                "cell_line_or_epigenome": rt.cell_line if rt else "",
                "assembly": rt.assembly if rt else "",
                "units": {"mrt": "MRT (0..1)",
                          "clr": "CLR(fraction)",
                          "wavelet": "signal"}[source]}
    if source == "atac":
        return {"description": "TCGA ATAC-seq insertions, mean over "
                "gene body or promoter window",
                "cell_line_or_epigenome": "TCGA tumor samples",
                "assembly": spec.atac.assembly if spec.atac else "",
                "units": "normalized insertions"}
    if source == "roadmap":
        return {"description": "Roadmap histone-mark ChIP fold-change "
                "signal, mean over gene body or promoter window",
                "cell_line_or_epigenome": "Roadmap epigenomes "
                + (", ".join(spec.roadmap.eids) if spec.roadmap else ""),
                "assembly": spec.roadmap.assembly if spec.roadmap else "",
                "units": "fold-change signal"}
    return {"description": source, "cell_line_or_epigenome": "",
            "assembly": "", "units": ""}


def _column_detail(source: str, column: str) -> str:
    if source in ("atac", "roadmap"):
        if column.endswith("_body"):
            return "gene body [TSS+200, TES]"
        if column.endswith("_promoter"):
            return "promoter [TSS-2000, TSS+200]"
    if source == "clr":
        return f"S-phase fraction {column.removeprefix('clr_rt_s')}"
    return ""


def _write_column_dictionary(
        paths: ProjectPaths,
        spec: ProjectSpec,
        blocks: list[tuple[str, list[str]]],
        final_columns: pd.Index,
        skew_report: pd.DataFrame | None) -> None:
    """One row per column of the built matrix (partial builds get a
    correct partial dictionary; fix-dropped columns are omitted)."""
    applied = {}
    if skew_report is not None:
        applied = skew_report.attrs.get("applied", {})
    final = set(final_columns)
    rows = []
    for source, columns in blocks:
        meta = _describe_block(source, spec, None)
        for column in columns:
            if column not in final:
                continue  # dropped by fix_variance
            transforms = applied.get(column)
            transform = ("; ".join(f"log({c:.6g} + x)"
                                   for c in transforms)
                         if transforms else "none")
            rows.append({
                "column": column,
                "source": source,
                "description": meta["description"],
                "detail": _column_detail(source, column),
                "cell_line_or_epigenome":
                    meta["cell_line_or_epigenome"],
                "assembly": meta["assembly"],
                "units": meta["units"],
                "transform": transform,
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(paths.column_dictionary_csv, index=False)
    logger.info("Wrote column dictionary (%d rows) to %s",
                len(frame), paths.column_dictionary_csv)


def _registry_hash(registry_path: str | Path | None) -> str:
    path = (Path(registry_path) if registry_path is not None
            else location_projects_registry)
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def _write_manifest(paths: ProjectPaths, spec: ProjectSpec,
                    included: set[str],
                    matrices: "CovariateMatrices",
                    registry_path: str | Path | None) -> None:
    try:
        from importlib.metadata import version
        pkg_version = version("sigmutselcovs")
    except Exception:  # noqa: BLE001 - version is informational only
        pkg_version = "unknown"
    manifest = {
        "project": spec.code,
        "sigmutselcovs_version": pkg_version,
        "registry_sha256_16": _registry_hash(registry_path),
        "sources_included": sorted(included),
        "shapes": {name: list(getattr(matrices, name).shape)
                   for name in ("full", "simple", "tcga")},
        "built": pd.Timestamp.now(tz="UTC").isoformat(),
    }
    paths.build_manifest_json.write_text(
        json.dumps(manifest, indent=2) + "\n")


def build_covariate_matrix(
        project: str,
        data_dir: str | Path,
        *,
        registry_path: str | Path | None = None,
        gencode_gtfs: Mapping[str, str | Path] | None = None,
        gtex_gct: str | Path | None = None,
        force_generation: bool = False,
        include: Sequence[str] | None = None,
        exclude: Sequence[str] = (),
        apply_fixes: bool = True,
        cache_matrices: bool = False,
        return_reports: bool = False,
        validate_output: bool = True,
        strict_validation: bool = False,
) -> CovariateMatrices | tuple[CovariateMatrices, dict]:
    """Build the covariate matrices for a registered project.

    Parameters
    ----------
    project : str
        TCGA study code registered in ``data/projects.json``
        (e.g. 'COAD', 'BRCA').
    data_dir : str | Path
        Project data directory laid out as `ProjectPaths` expects
        (the historical ``coad_data/`` structure).
    registry_path : str | Path | None
        Alternative registry JSON; defaults to the packaged one.
    gencode_gtfs : Mapping[str, str | Path] | None
        Assembly -> GENCODE GTF path, e.g. ``{"hg19": ..., "hg38":
        ...}``.  Defaults to the GTFs bundled with sigmutsel.  Each
        source block picks its GTF by its registry ``assembly``.
    gtex_gct : str | Path | None
        GTEx gene-median-TPM GCT; defaults to the packaged location.
    force_generation : bool
        Recompute the per-gene caches instead of loading them.
    include : Sequence[str] | None
        Sources to include (subset of gtex, gexp, repliseq, roadmap,
        atac).  None means all sources the registry defines.
    exclude : Sequence[str]
        Sources to drop from the selection.
    apply_fixes : bool
        Run `fix_all` (variance + skewness) on the concatenated
        matrix, as the COAD pipeline always did.
    cache_matrices : bool
        Write the three matrices and a build manifest under
        ``<data_dir>/covariate_matrices/``.
    return_reports : bool
        Also return the `fix_all` reports dict.
    validate_output : bool
        Run `sigmutselcovs.validate.validate_covariates` on the raw
        matrix at the end of the build (results are logged).
    strict_validation : bool
        Raise when a validation check fails.

    Returns
    -------
    CovariateMatrices, or (CovariateMatrices, reports) when
    ``return_reports`` — a NamedTuple that unpacks as the historical
    ``(cov_matrix_full, cov_matrix_simple, cov_matrix_tcga)``.
    """
    spec = get_project(project, registry_path)
    paths = project_paths(data_dir)
    selected = _selected(include, exclude)
    gtfs = dict(gencode_gtfs) if gencode_gtfs is not None else None

    def gtf_for(assembly: str) -> str | Path:
        nonlocal gtfs
        if gtfs is None:
            gtfs = _default_gencode_gtfs()
        if assembly not in gtfs:
            raise ValueError(
                f"No GENCODE GTF configured for assembly {assembly!r}")
        return gtfs[assembly]

    # Concatenation order must stay identical to the historical
    # covariates.py: gtex, mean tumor, mean normal, per-sample, mrt,
    # clr fractions, ATAC, Roadmap.
    frames: list[pd.DataFrame] = []
    included: set[str] = set()
    blocks: list[tuple[str, list[str]]] = []

    if "gtex" in selected:
        from .download import resolve_gtex_gct
        gtex = import_gtex(
            resolve_gtex_gct(gtex_gct),
            mapping_path=location_gtex_tcga_mapping,
            columns=spec.gtex.mapping_key)
        frames.append(gtex)
        blocks.append(("gtex", list(gtex.columns)))
        included.add("gtex")

    if "gexp" in selected:
        if spec.gexp is None:
            _skip("gexp", "not defined in the registry for "
                  f"{spec.code}")
        else:
            metrics = (list(spec.gexp.metrics)
                       if spec.gexp.metrics is not None else None)
            for tissue_type in spec.gexp.tissue_types:
                location_csv = (paths.mean_gexp_csv
                                if tissue_type == "Tumor"
                                else paths.mean_gexp_normal_csv)
                mean_gexp = load_or_generate_mean_tcga_gexp(
                    location_csv=location_csv,
                    tcga_dir=paths.gexp_tcga_dir,
                    tissue_type=tissue_type,
                    force_generation=force_generation)
                if tissue_type != "Tumor":
                    mean_gexp = mean_gexp.rename(
                        f"{mean_gexp.name}_{tissue_type.lower()}")
                frames.append(mean_gexp.to_frame())
                blocks.append(("gexp_mean", [mean_gexp.name]))
            per_sample = load_or_generate_tcga_gexp_per_sample(
                location_parquet=paths.gexp_per_sample_parquet,
                tcga_dir=paths.gexp_tcga_dir,
                metrics=metrics,
                force_generation=force_generation)
            frames.append(per_sample)
            blocks.append(("gexp_per_sample",
                           list(per_sample.columns)))
            included.add("gexp")

    if "repliseq" in selected:
        if spec.repliseq is None:
            _skip("repliseq", "not defined in the registry for "
                  f"{spec.code}")
        else:
            rt_frames = _load_repliseq(
                spec.repliseq, paths,
                gtf_for(spec.repliseq.assembly),
                force_generation)
            if rt_frames:
                frames.extend(rt_frames)
                for frame in rt_frames:
                    first = frame.columns[0]
                    kind = ("clr" if first.startswith("clr_")
                            else "wavelet" if first == "rt_wavelet"
                            else "mrt")
                    blocks.append((kind, list(frame.columns)))
                included.add("repliseq")

    atac_covs = None
    if "atac" in selected:
        if spec.atac is None:
            _skip("atac", f"not defined in the registry for {spec.code}")
        else:
            atac_covs = _load_chromatin(
                "atac", paths.atac_covs_csv, paths.atac_bigwig_files(),
                gtf_for(spec.atac.assembly), force_generation)
            if atac_covs is not None:
                frames.append(atac_covs)
                blocks.append(("atac", list(atac_covs.columns)))
                included.add("atac")

    if "roadmap" in selected:
        if spec.roadmap is None:
            _skip("roadmap", "not defined in the registry for "
                  f"{spec.code}")
        else:
            roadmap_covs = _load_chromatin(
                "roadmap", paths.roadmap_covs_csv,
                paths.roadmap_bigwig_files(),
                gtf_for(spec.roadmap.assembly), force_generation)
            if roadmap_covs is not None:
                frames.append(roadmap_covs)
                blocks.append(("roadmap",
                               list(roadmap_covs.columns)))
                included.add("roadmap")

    if not frames:
        raise ValueError(f"No covariate sources available for "
                         f"{spec.code} in {paths.root}")

    cov_matrix_full_raw = pd.concat(frames, axis=1, join="outer")

    duplicated = cov_matrix_full_raw.columns[
        cov_matrix_full_raw.columns.duplicated()].unique()
    if len(duplicated):
        logger.warning("Duplicated covariate columns after concat: %s",
                       list(duplicated))

    reports: dict = {}
    if apply_fixes:
        cov_matrix_full, reports = fix_all(cov_matrix_full_raw)
    else:
        cov_matrix_full = cov_matrix_full_raw

    # Simple interpretable model: one variable per major covariate
    # type. H3K4me3 marks active transcription (lower mutation rate);
    # H3K9me3 marks heterochromatin/late replication (higher mutation
    # rate). Each of these is a mean across several individual raw
    # tracks (tissues/samples) that are themselves columns of
    # cov_matrix_full_raw, so fix_all above already ran the skewness
    # test and applied log(pseudo_count + x) to each one individually,
    # before this averaging step. Do not log-transform again here:
    # these columns can be negative post-transform, and applying
    # log(1 + x) a second time is undefined once x < -1, which is
    # exactly what produced NaN for TTN and MUC16's gene expression.
    simple_cols: dict[str, pd.Series] = {}
    sm = spec.simple_matrix
    if sm.gtex_column in cov_matrix_full.columns:
        simple_cols[sm.gtex_column] = cov_matrix_full[sm.gtex_column]
    else:
        logger.warning("cov_matrix_simple: %r not in the full matrix "
                       "(dropped by fixes or gtex not built); omitting",
                       sm.gtex_column)
    if sm.include_mrt:
        if "mrt" in cov_matrix_full.columns:
            simple_cols["mrt"] = cov_matrix_full["mrt"]
        else:
            logger.warning("cov_matrix_simple: no mrt column; omitting")
    body_cols = [c for c in cov_matrix_full.columns
                 if sm.body_pattern in c]
    if body_cols:
        simple_cols["h3k4me3_body"] = (
            cov_matrix_full[body_cols].mean(axis=1))
    else:
        logger.warning("cov_matrix_simple: no columns match %r; omitting",
                       sm.body_pattern)
    promoter_cols = [c for c in cov_matrix_full.columns
                     if sm.promoter_pattern in c]
    if promoter_cols:
        simple_cols["h3k9me3_prom"] = (
            cov_matrix_full[promoter_cols].mean(axis=1))
    else:
        logger.warning("cov_matrix_simple: no columns match %r; omitting",
                       sm.promoter_pattern)
    cov_matrix_simple = pd.DataFrame(simple_cols)

    # TCGA-only matrix: gene expression + ATAC-seq
    tcga_cols: dict[str, pd.Series] = {}
    if "tpm_unstranded" in cov_matrix_full.columns:
        tcga_cols["gexp"] = cov_matrix_full["tpm_unstranded"]
    else:
        logger.warning("cov_matrix_tcga: no tpm_unstranded column; "
                       "omitting gexp")
    if atac_covs is not None and spec.atac is not None:
        prefix = _atac_prefix(spec.atac, atac_covs) + "_"
        for suffix, name in (("body", "atac_body"),
                             ("promoter", "atac_promoter")):
            cols = [c for c in cov_matrix_full.columns
                    if c.startswith(prefix) and c.endswith(suffix)]
            if cols:
                tcga_cols[name] = cov_matrix_full[cols].mean(axis=1)
            else:
                logger.warning("cov_matrix_tcga: no %s* columns ending "
                               "in %s; omitting %s", prefix, suffix, name)
    cov_matrix_tcga = pd.DataFrame(tcga_cols)

    matrices = CovariateMatrices(full=cov_matrix_full,
                                 simple=cov_matrix_simple,
                                 tcga=cov_matrix_tcga)

    if cache_matrices:
        paths.matrices_dir.mkdir(parents=True, exist_ok=True)
        matrices.full.to_parquet(paths.matrix_full_parquet)
        matrices.simple.to_csv(paths.matrix_simple_csv)
        matrices.tcga.to_csv(paths.matrix_tcga_csv)
        _write_manifest(paths, spec, included, matrices, registry_path)
        _write_column_dictionary(paths, spec, blocks,
                                 cov_matrix_full.columns,
                                 reports.get("skewness"))
        logger.info("Cached covariate matrices under %s",
                    paths.matrices_dir)

    if validate_output:
        from .validate import validate_covariates
        validate_covariates(spec.code,
                            cov_matrix_raw=cov_matrix_full_raw,
                            registry_path=registry_path,
                            strict=strict_validation)

    if return_reports:
        return matrices, reports
    return matrices
