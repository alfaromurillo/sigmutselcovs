"""Per-project data-directory layout.

`ProjectPaths` maps a project data directory (e.g. ``coad_data/``,
``brca_data/``) to the standard locations used by the download layer
and the covariate-matrix builder.  The layout mirrors the historical
``coad_analysis/coad_data/`` tree exactly, so existing COAD caches
are reused without moving any file.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

# Order matters: it determines bigWig track order and therefore the
# column order of the chromatin blocks in the covariate matrix.
_BIGWIG_PATTERNS = [
    "*.bigwig",
    "*.bigWig",
    "*.BigWig",
    "*.bw",
    "*.BW",
]


def glob_patterns_dedup(
    root: Path, patterns: Sequence[str]
) -> list[Path]:
    """Return a de-duplicated list of files matching patterns.

    Preserves the order produced by iterating over patterns, while
    removing duplicates when a file matches multiple patterns.

    Parameters
    ----------
    root : pathlib.Path
        Directory to search in.
    patterns : sequence of str
        Glob patterns, e.g., ["*.bigwig", "*.bw"].

    Returns
    -------
    list[pathlib.Path]
        Ordered, de-duplicated list of matching files.
    """
    files: list[Path] = []
    for pat in patterns:
        files.extend(sorted(root.glob(pat)))
    out: list[Path] = []
    seen: set[Path] = set()
    for p in files:
        if p in seen:
            continue
        out.append(p)
        seen.add(p)
    return out


def bigwig_files(directory: str | Path) -> list[Path]:
    """Return the bigWig files in a directory, in canonical order.

    The order (pattern list, then sorted within each pattern) is the
    one the historical COAD pipeline used, and covariate column order
    depends on it — do not change.
    """
    directory = Path(directory)
    if not directory.exists():
        return []
    return glob_patterns_dedup(directory, _BIGWIG_PATTERNS)


@dataclass(frozen=True, slots=True)
class ProjectPaths:
    """Standard locations inside a project data directory."""

    root: Path

    # --- gene expression ---
    @property
    def gexp_dir(self) -> Path:
        return self.root / "gene_expression"

    @property
    def gexp_tcga_dir(self) -> Path:
        """GDC manifests + sample sheets live at this level."""
        return self.gexp_dir / "tcga"

    @property
    def gexp_star_dir(self) -> Path:
        """STAR count files, one UUID subdirectory per file."""
        return self.gexp_tcga_dir / "star_gene_counts"

    @property
    def mean_gexp_csv(self) -> Path:
        return self.gexp_dir / "tcga_mean_gene_expression.csv"

    @property
    def mean_gexp_normal_csv(self) -> Path:
        return self.gexp_dir / "tcga_mean_gene_expression_normal.csv"

    @property
    def gexp_per_sample_parquet(self) -> Path:
        return self.gexp_dir / "tcga_gexp_per_sample.parquet"

    # --- replication timing ---
    @property
    def rt_dir(self) -> Path:
        return self.root / "replication_timing"

    @property
    def rt_encode_dir(self) -> Path:
        """Downloaded ENCODE fraction/wavelet bigWigs."""
        return self.rt_dir / "encode"

    @property
    def mrt_csv(self) -> Path:
        return self.rt_dir / "mean_replication_time_per_gene.csv"

    @property
    def rt_fractions_csv(self) -> Path:
        return self.rt_dir / "rt_fractions_per_gene.csv"

    @property
    def rt_wavelet_csv(self) -> Path:
        return self.rt_dir / "rt_wavelet_per_gene.csv"

    # --- chromatin ---
    @property
    def roadmap_dir(self) -> Path:
        return self.root / "chromatin" / "roadmap"

    @property
    def roadmap_covs_csv(self) -> Path:
        return self.roadmap_dir / "chromatin_covs_roadmap.csv"

    @property
    def atac_dir(self) -> Path:
        return self.root / "chromatin" / "tcga"

    @property
    def atac_covs_csv(self) -> Path:
        return self.atac_dir / "chromatin_covs_tcga.csv"

    @property
    def encode_chromatin_dir(self) -> Path:
        """Downloaded ENCODE histone/DNase bigWigs (not repli-seq --
        see rt_encode_dir for those, which live under a separate
        replication_timing/encode/ tree)."""
        return self.root / "chromatin" / "encode"

    @property
    def encode_chromatin_covs_csv(self) -> Path:
        return self.encode_chromatin_dir / "chromatin_covs_encode.csv"

    # --- MAF files (downloaded by sigmutsel, not this package) ---
    @property
    def maf_dir(self) -> Path:
        return self.root / "tcga" / "all_maf_files"

    # --- built covariate matrices ---
    @property
    def matrices_dir(self) -> Path:
        return self.root / "covariate_matrices"

    @property
    def matrix_full_parquet(self) -> Path:
        return self.matrices_dir / "cov_matrix_full.parquet"

    @property
    def matrix_simple_csv(self) -> Path:
        return self.matrices_dir / "cov_matrix_simple.csv"

    @property
    def matrix_tcga_csv(self) -> Path:
        return self.matrices_dir / "cov_matrix_tcga.csv"

    @property
    def build_manifest_json(self) -> Path:
        return self.matrices_dir / "build_manifest.json"

    @property
    def matrix_pca_parquet(self) -> Path:
        return self.matrices_dir / "cov_matrix_pca.parquet"

    @property
    def pca_manifest_json(self) -> Path:
        return self.matrices_dir / "pca_manifest.json"

    @property
    def column_dictionary_csv(self) -> Path:
        return self.matrices_dir / "cov_matrix_columns.csv"

    def roadmap_bigwig_files(self) -> list[Path]:
        return bigwig_files(self.roadmap_dir)

    def atac_bigwig_files(self) -> list[Path]:
        return bigwig_files(self.atac_dir)

    def encode_chromatin_bigwig_files(self) -> list[Path]:
        return bigwig_files(self.encode_chromatin_dir)


def project_paths(data_dir: str | Path) -> ProjectPaths:
    """Return the standard layout for a project data directory."""
    return ProjectPaths(root=Path(data_dir).expanduser().resolve())
