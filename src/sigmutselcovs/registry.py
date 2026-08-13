"""Per-cancer-type registry of covariate data sources.

The registry maps a TCGA study code (e.g. ``COAD``, ``BRCA``) to the
data sources needed to download and build its covariate matrix: GTEx
tissue columns, TCGA gene expression (GDC query parameters), TCGA
ATAC-seq (GDC tarball UUID), Roadmap Epigenomics chromatin tracks,
and replication timing (one of three source types).

The registry lives in ``data/projects.json``, which has a
``defaults`` block (merged under every project's per-source entries)
and a ``projects`` block with one row per study code.  Any of the
sources ``gexp``, ``atac``, ``roadmap``, ``repliseq`` may be ``null``
for a project, in which case the builder skips that block; ``gtex``
is required.

Replication timing source types:

- ``mat``: a single transposed bins x fractions text table (the
  16-fraction Gaussian-smoothed Repli-seq MAT files from GEO
  GSE137764, e.g. HCT116 for COAD).
- ``fraction_bigwigs``: N >= 2 per-fraction bigWig tracks ordered
  early to late (the UW 6-fraction Repli-seq on ENCODE, e.g. MCF-7
  for BRCA).
- ``wavelet``: a single wavelet-smoothed early/late signal bigWig
  (the Gilbert-lab ENCODE series, e.g. LNCaP, Caki2, A549); yields
  one per-gene column, no MRT or CLR fractions.
"""

import json
import logging
from dataclasses import dataclass, fields
from pathlib import Path

from .covariates_locations import (
    location_covariates_data,
    location_gtex_tcga_mapping,
)

logger = logging.getLogger(__name__)

location_projects_registry = (
    location_covariates_data / "projects.json"
)

_REPLISEQ_TYPES = ("mat", "fraction_bigwigs", "wavelet")
_OPTIONAL_SOURCES = ("gexp", "atac", "roadmap", "repliseq")


@dataclass(frozen=True, slots=True)
class GtexSpec:
    """GTEx gene expression source (always required)."""

    mapping_key: str
    representative_column: str


@dataclass(frozen=True, slots=True)
class GexpSpec:
    """TCGA gene expression source (GDC files-API query)."""

    tcga_project_id: str
    data_type: str = "Gene Expression Quantification"
    workflow_type: str = "STAR - Counts"
    metrics: tuple[str, ...] | None = None
    tissue_types: tuple[str, ...] = ("Tumor", "Normal")


@dataclass(frozen=True, slots=True)
class AtacSpec:
    """TCGA ATAC-seq source (Corces et al. bigWig tarball on GDC)."""

    gdc_uuid: str
    column_prefix: str
    assembly: str = "hg38"


@dataclass(frozen=True, slots=True)
class RoadmapSpec:
    """Roadmap Epigenomics chromatin fold-change signal tracks."""

    eids: tuple[str, ...]
    marks: tuple[str, ...]
    required_marks: tuple[str, ...] = ()
    assembly: str = "hg19"
    url_template: str = (
        "https://egg2.wustl.edu/roadmap/data/byFileType/signal/"
        "consolidated/macs2signal/foldChange/"
        "{eid}-{mark}.fc.signal.bigwig"
    )


@dataclass(frozen=True, slots=True)
class TrackRef:
    """One replication-timing fraction track (ENCODE accession)."""

    label: str
    accession: str
    source: str = "encode"
    url: str | None = None


@dataclass(frozen=True, slots=True)
class RepliseqSpec:
    """Replication timing source."""

    type: str
    assembly: str
    cell_line: str
    filename: str | None = None
    url: str | None = None
    tracks: tuple[TrackRef, ...] = ()
    n_fractions: int | None = None
    mrt_fraction_cols: tuple[str, ...] | None = None
    include_clr_fractions: bool | None = None
    bin_size: int = 50_000


@dataclass(frozen=True, slots=True)
class SimpleMatrixSpec:
    """Recipe for the reduced covariate matrix (cov_matrix_simple)."""

    gtex_column: str
    body_pattern: str = "h3k4me3_fc_signal_body"
    promoter_pattern: str = "h3k9me3_fc_signal_promoter"
    include_mrt: bool = True


@dataclass(frozen=True, slots=True)
class ProjectSpec:
    """All covariate sources for one TCGA study code."""

    code: str
    description: str
    gtex: GtexSpec
    simple_matrix: SimpleMatrixSpec
    gexp: GexpSpec | None = None
    atac: AtacSpec | None = None
    roadmap: RoadmapSpec | None = None
    repliseq: RepliseqSpec | None = None


def _tupled(value):
    """Convert lists (from JSON) to tuples for the frozen dataclasses."""
    if isinstance(value, list):
        return tuple(_tupled(v) for v in value)
    return value


def _build_spec(cls, raw: dict):
    """Instantiate a spec dataclass from a raw dict, tupling lists."""
    known = {f.name for f in fields(cls)}
    unknown = set(raw) - known
    if unknown:
        raise ValueError(
            f"Unknown keys for {cls.__name__}: {sorted(unknown)}"
        )
    kwargs = {k: _tupled(v) for k, v in raw.items()}
    if cls is RepliseqSpec and kwargs.get("tracks"):
        kwargs["tracks"] = tuple(
            TrackRef(**dict(t)) if isinstance(t, dict) else t
            for t in raw["tracks"]
        )
    return cls(**kwargs)


def _merge_defaults(
    entry: dict | None, defaults: dict | None
) -> dict | None:
    """Merge a project's per-source entry over the defaults block.

    ``entry`` is None when the project explicitly disables the source
    (JSON null); an absent source key in the project row also yields
    None.  Defaults only apply when the project provides an entry.
    """
    if entry is None:
        return None
    merged = dict(defaults or {})
    merged.update(entry)
    return merged


def validate_registry(raw: dict) -> None:
    """Validate the raw registry dict; raise ValueError on problems."""
    if raw.get("schema_version") != 1:
        raise ValueError(
            "Unsupported registry schema_version: "
            f"{raw.get('schema_version')!r} (expected 1)"
        )
    projects = raw.get("projects")
    if not projects:
        raise ValueError("Registry has no projects")

    with open(location_gtex_tcga_mapping) as fh:
        gtex_mapping = {str(k).upper() for k in json.load(fh)}

    for code, row in projects.items():
        if row.get("gtex") is None:
            raise ValueError(f"{code}: gtex source is required")
        mapping_key = row["gtex"].get("mapping_key", "").upper()
        if mapping_key not in gtex_mapping:
            raise ValueError(
                f"{code}: gtex.mapping_key {mapping_key!r} not in "
                f"{location_gtex_tcga_mapping.name}"
            )
        repliseq = row.get("repliseq")
        if repliseq is not None:
            rtype = repliseq.get("type")
            if rtype not in _REPLISEQ_TYPES:
                raise ValueError(
                    f"{code}: repliseq.type {rtype!r} not one of "
                    f"{_REPLISEQ_TYPES}"
                )
            if rtype == "mat" and not repliseq.get("filename"):
                raise ValueError(
                    f"{code}: repliseq type 'mat' needs a filename"
                )
            if rtype == "fraction_bigwigs" and not repliseq.get(
                "tracks"
            ):
                raise ValueError(
                    f"{code}: repliseq type 'fraction_bigwigs' "
                    "needs tracks"
                )
            if rtype == "wavelet" and not repliseq.get("tracks"):
                raise ValueError(
                    f"{code}: repliseq type 'wavelet' needs one track"
                )


def load_registry(
    path: str | Path | None = None,
) -> dict[str, ProjectSpec]:
    """Load and validate the project registry.

    Parameters
    ----------
    path : str | Path | None
        Registry JSON to load; defaults to the packaged
        ``data/projects.json``.

    Returns
    -------
    dict[str, ProjectSpec]
        Study code -> project specification.
    """
    path = (
        Path(path) if path is not None else location_projects_registry
    )
    with open(path) as fh:
        raw = json.load(fh)
    validate_registry(raw)

    defaults = raw.get("defaults", {})
    registry = {}
    for code, row in raw["projects"].items():
        code = code.upper()

        def source(key, cls, *, row=row):
            merged = _merge_defaults(row.get(key), defaults.get(key))
            return (
                None if merged is None else _build_spec(cls, merged)
            )

        gtex = _build_spec(GtexSpec, row["gtex"])
        simple_raw = _merge_defaults(
            row.get("simple_matrix", {}),
            defaults.get("simple_matrix"),
        )
        simple_raw.setdefault(
            "gtex_column", gtex.representative_column
        )
        registry[code] = ProjectSpec(
            code=code,
            description=row.get("description", code),
            gtex=gtex,
            simple_matrix=_build_spec(SimpleMatrixSpec, simple_raw),
            gexp=source("gexp", GexpSpec),
            atac=source("atac", AtacSpec),
            roadmap=source("roadmap", RoadmapSpec),
            repliseq=source("repliseq", RepliseqSpec),
        )
    return registry


def available_projects(path: str | Path | None = None) -> list[str]:
    """Return the study codes available in the registry, sorted."""
    return sorted(load_registry(path))


def get_project(
    code: str, path: str | Path | None = None
) -> ProjectSpec:
    """Return the spec for one study code.

    Raises ValueError (not KeyError) for unknown codes, listing the
    available ones.
    """
    registry = load_registry(path)
    code_up = code.upper()
    if code_up not in registry:
        raise ValueError(
            f"Unknown project {code!r}; available: "
            f"{', '.join(sorted(registry))}"
        )
    return registry[code_up]
