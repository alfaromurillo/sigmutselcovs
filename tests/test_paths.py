"""Tests for the per-project data-directory layout."""

from pathlib import Path

from sigmutselcovs.paths import (
    ProjectPaths,
    bigwig_files,
    project_paths,
)


def test_layout_matches_historical_coad_tree(tmp_path):
    """Literal-string test: this is what protects COAD cache reuse."""
    p = project_paths(tmp_path)
    root = Path(tmp_path).resolve()
    assert p.gexp_tcga_dir == root / "gene_expression" / "tcga"
    assert p.gexp_star_dir == (
        root / "gene_expression" / "tcga" / "star_gene_counts"
    )
    assert p.mean_gexp_csv == (
        root / "gene_expression" / "tcga_mean_gene_expression.csv"
    )
    assert p.mean_gexp_normal_csv == (
        root
        / "gene_expression"
        / "tcga_mean_gene_expression_normal.csv"
    )
    assert p.gexp_per_sample_parquet == (
        root / "gene_expression" / "tcga_gexp_per_sample.parquet"
    )
    assert p.mrt_csv == (
        root
        / "replication_timing"
        / "mean_replication_time_per_gene.csv"
    )
    assert p.rt_fractions_csv == (
        root / "replication_timing" / "rt_fractions_per_gene.csv"
    )
    assert p.roadmap_covs_csv == (
        root / "chromatin" / "roadmap" / "chromatin_covs_roadmap.csv"
    )
    assert p.atac_covs_csv == (
        root / "chromatin" / "tcga" / "chromatin_covs_tcga.csv"
    )
    assert p.maf_dir == root / "tcga" / "all_maf_files"
    assert p.matrices_dir == root / "covariate_matrices"


def test_bigwig_files_order_and_dedup(tmp_path):
    """Pattern order first, sorted within each pattern, no dupes."""
    names = ["b.bigwig", "a.bigwig", "z.bw", "c.bigWig", "y.bw"]
    for name in names:
        (tmp_path / name).touch()
    got = [p.name for p in bigwig_files(tmp_path)]
    assert got == ["a.bigwig", "b.bigwig", "c.bigWig", "y.bw", "z.bw"]


def test_bigwig_files_missing_dir(tmp_path):
    assert bigwig_files(tmp_path / "nope") == []


def test_paths_frozen(tmp_path):
    p = project_paths(tmp_path)
    try:
        p.root = Path("/elsewhere")
        raised = False
    except AttributeError:
        raised = True
    assert raised


def test_project_paths_accepts_str(tmp_path):
    assert project_paths(str(tmp_path)) == ProjectPaths(
        root=Path(tmp_path).resolve()
    )
