"""Builder: chromatin windows, log-means and peaks end to end.

Real (synthetic) bigWigs, peak files and GTF from ``conftest.py``,
built through `build_covariate_matrix`; no network.
"""

import json

import numpy as np
import pandas as pd
import pytest
from conftest import (
    GENE_IDS,
    legacy_registry_raw,
    write_bigwig,
    write_gtf,
    write_peaks,
)

from sigmutselcovs import builder
from sigmutselcovs import covariates_chromatin as cc
from sigmutselcovs.builder import (
    build_covariate_matrix,
    covariate_column_blocks,
)
from sigmutselcovs.paths import project_paths

# in file-name order, which is column order
TRACKS = ("E075-H3K4me3", "E075-H3K9me3", "E106-H3K4me3")


@pytest.fixture
def project(tmp_path):
    """COAD's data dir with synthetic Roadmap bigWigs and peaks."""
    paths = project_paths(tmp_path / "data")
    paths.roadmap_dir.mkdir(parents=True)
    paths.roadmap_peaks_dir.mkdir()
    for seed, track in enumerate(TRACKS):
        write_bigwig(
            paths.roadmap_dir / f"{track}.fc.signal.bigwig", seed
        )
        write_peaks(
            paths.roadmap_peaks_dir / f"{track}.narrowPeak.gz",
            gzipped=True,
        )
    gtf = write_gtf(tmp_path / "synthetic.gtf")
    return paths, {"hg19": gtf, "hg38": gtf}


def _registry(tmp_path, **roadmap):
    raw = legacy_registry_raw()
    raw["projects"]["COAD"]["roadmap"].update(roadmap)
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    return path


def _build(paths, gtfs, registry_path=None, **kwargs):
    kwargs = {
        "include": ["roadmap"],
        "apply_fixes": False,
        "validate_output": False,
        **kwargs,
    }
    return build_covariate_matrix(
        "COAD",
        paths.root,
        registry_path=registry_path,
        gencode_gtfs=gtfs,
        **kwargs,
    )


def _caches(paths):
    return sorted(
        p.name for p in paths.roadmap_dir.glob("chromatin_covs_*.csv")
    )


def test_default_build_is_unchanged(project):
    paths, gtfs = project
    full = _build(paths, gtfs).full
    # what the builder computed before the options existed
    expected = cc.load_or_generate_chromatin_covariates(
        paths.root / "reference.csv",
        paths.roadmap_bigwig_files(),
        gtfs["hg19"],
        biotypes=("protein_coding",),
        include_promoter=True,
        promoter_upstream=2000,
        promoter_downstream=200,
    )
    pd.testing.assert_frame_equal(full, expected, check_exact=True)
    assert list(full.columns) == [
        f"{t.lower().replace('-', '_')}_fc_signal_{w}"
        for t in TRACKS
        for w in ("body", "promoter")
    ]
    assert _caches(paths) == ["chromatin_covs_roadmap.csv"]
    assert set(full.attrs["column_sources"].values()) == {"roadmap"}


def test_explicit_defaults_hit_the_default_cache(project, tmp_path):
    paths, gtfs = project
    _build(paths, gtfs)
    # read back from the cache, as the explicit build below is (a
    # fresh build can differ in the last ulp: pandas' CSV parser)
    default = _build(paths, gtfs).full
    registry = _registry(
        tmp_path,
        windows=["body", "promoter"],
        transform=None,
        peaks=[],
    )
    explicit = _build(paths, gtfs, registry_path=registry).full
    pd.testing.assert_frame_equal(default, explicit, check_exact=True)
    assert _caches(paths) == ["chromatin_covs_roadmap.csv"]


def test_windows_log_and_peaks(project, tmp_path):
    paths, gtfs = project
    registry = _registry(
        tmp_path,
        windows=["cds_span", "cds_span2k"],
        transform="log1p",
        peaks=["narrowPeak"],
    )
    full = _build(paths, gtfs, registry_path=registry).full
    signal = [
        f"{t.lower().replace('-', '_')}_fc_signal_{w}_log"
        for t in TRACKS
        for w in ("cds_span", "cds_span2k")
    ]
    peaks = [
        f"{t.lower().replace('-', '_')}_narrowpeak_{w}"
        for t in TRACKS
        for w in ("cds_span", "cds_span2k")
    ]
    assert list(full.columns) == signal + peaks
    assert list(full.index) == GENE_IDS
    sources = full.attrs["column_sources"]
    assert {sources[c] for c in signal} == {"roadmap"}
    assert {sources[c] for c in peaks} == {"roadmap_peaks"}
    assert set(covariate_column_blocks(full)) == {"chromatin"}

    # values are the summaries' own
    bigwig = paths.roadmap_dir / "E106-H3K4me3.fc.signal.bigwig"
    genes = cc._load_gene_bodies(
        gtfs["hg19"], None, with_cds_span=True
    )
    genes = genes.loc[GENE_IDS]
    direct = cc.summarize_bigwig_windows(
        bigwig,
        genes,
        windows=["cds_span", "cds_span2k"],
        transform="log1p",
        label="E106-H3K4me3.fc.signal",
    )
    pd.testing.assert_frame_equal(
        full[direct.columns], direct, rtol=1e-12
    )
    # gene 3 has no CDS
    assert full.loc["ENSG00000000003"].isna().all()
    assert full.loc["ENSG00000000001", peaks[0]] == 500 / 3200

    assert _caches(paths) == [
        "chromatin_covs_roadmap_cds_span-cds_span2k_log1p.csv",
        "chromatin_covs_roadmap_peaks_narrowpeak_cds_span-cds_span2k.csv",
    ]


def test_other_options_do_not_read_the_default_cache(
    project, tmp_path, monkeypatch
):
    paths, gtfs = project
    _build(paths, gtfs)
    # mark the default cache, so reading it shows
    cache = paths.roadmap_covs_csv
    marked = pd.read_csv(cache, index_col=0) * 0 - 1.0
    marked.to_csv(cache)

    log_registry = _registry(tmp_path, transform="log1p")
    log = _build(paths, gtfs, registry_path=log_registry).full
    assert all(c.endswith("_log") for c in log.columns)
    assert (log.fillna(0) >= 0).all().all()

    # defaults still read their (marked) cache
    again = _build(paths, gtfs).full
    assert (again.dropna() == -1.0).all().all()

    # and the log build its own, without recomputing
    def fail(*args, **kwargs):
        raise AssertionError("recomputed a cached summary")

    monkeypatch.setattr(cc, "summarise_tracks_to_genes", fail)
    cached = _build(paths, gtfs, registry_path=log_registry).full
    pd.testing.assert_frame_equal(log, cached, rtol=1e-12)


def test_peak_cache_follows_its_files(project, tmp_path):
    paths, gtfs = project
    registry = _registry(tmp_path, peaks=["narrowPeak"])
    first = _build(paths, gtfs, registry_path=registry).full
    gone = paths.roadmap_peaks_dir / "E075-H3K9me3.narrowPeak.gz"
    gone.unlink()
    second = _build(paths, gtfs, registry_path=registry).full
    dropped = set(first.columns) - set(second.columns)
    assert dropped == {
        "e075_h3k9me3_narrowpeak_body",
        "e075_h3k9me3_narrowpeak_promoter",
    }


def test_peaks_without_bigwigs(project, tmp_path, caplog):
    paths, gtfs = project
    for path in paths.roadmap_bigwig_files():
        path.unlink()
    registry = _registry(tmp_path, peaks=["narrowPeak"])
    matrices = _build(paths, gtfs, registry_path=registry)
    assert all("narrowpeak" in c for c in matrices.full.columns)
    assert any(
        "Skipping roadmap covariates" in m for m in caplog.messages
    )


def test_average_by_assay_with_options(project, tmp_path):
    paths, gtfs = project
    raw = legacy_registry_raw()
    raw["projects"]["COAD"]["roadmap"].update(
        windows=["body", "cds_span"],
        transform="log1p",
        peaks=["narrowPeak"],
    )
    raw["projects"]["COAD"]["average_by_assay"] = True
    registry = tmp_path / "projects.json"
    registry.write_text(json.dumps(raw))
    full = _build(paths, gtfs, registry_path=registry).full
    assert list(full.columns) == [
        "h3k4me3_body_log",
        "h3k4me3_cds_span_log",
        "h3k9me3_body_log",
        "h3k9me3_cds_span_log",
        "h3k4me3_narrowpeak_body",
        "h3k4me3_narrowpeak_cds_span",
        "h3k9me3_narrowpeak_body",
        "h3k9me3_narrowpeak_cds_span",
    ]
    signal = cc.summarise_tracks_to_genes(
        [
            p
            for p in paths.roadmap_bigwig_files()
            if "H3K4me3" in p.name
        ],
        cc._load_gene_bodies(
            gtfs["hg19"], None, with_cds_span=True
        ).loc[GENE_IDS],
        windows=["body"],
        transform="log1p",
    )
    np.testing.assert_allclose(
        full["h3k4me3_body_log"], signal.mean(axis=1), rtol=1e-12
    )


def test_column_dictionary_and_validation(project, tmp_path):
    paths, gtfs = project
    registry = _registry(
        tmp_path,
        windows=["promoter", "cds_span2k"],
        transform="log1p",
        peaks=["narrowPeak"],
    )
    _build(
        paths,
        gtfs,
        registry_path=registry,
        cache_matrices=True,
        validate_output=True,
        strict_validation=True,
    )
    dictionary = pd.read_csv(paths.column_dictionary_csv).set_index(
        "column"
    )
    log_row = dictionary.loc["e075_h3k4me3_fc_signal_cds_span2k_log"]
    assert log_row["source"] == "roadmap"
    assert log_row["detail"] == (
        "CDS span +-2 kb, mean of log1p(signal) per base"
    )
    assert "log1p" in log_row["description"]
    peak_row = dictionary.loc["e075_h3k4me3_narrowpeak_promoter"]
    assert peak_row["source"] == "roadmap_peaks"
    assert peak_row["detail"] == (
        "promoter [TSS-2000, TSS+200], fraction covered by "
        "narrowPeak peaks"
    )
    assert peak_row["units"] == "fraction of bases"


def test_simple_matrix_ignores_log_columns(project, tmp_path):
    paths, gtfs = project
    registry = _registry(tmp_path, transform="log1p")
    matrices = _build(paths, gtfs, registry_path=registry)
    assert "h3k4me3_body" not in matrices.simple.columns
    default = _build(paths, gtfs)
    assert "h3k4me3_body" in default.simple.columns


def test_load_chromatin_signature_kept(tmp_path, monkeypatch):
    """The default call passes no new options to the loader."""
    seen = {}

    def fake(covs_csv, bigwigs, gtf, **kwargs):
        seen.update(kwargs)
        return pd.DataFrame(index=["ENSG_A"])

    monkeypatch.setattr(
        builder, "load_or_generate_chromatin_covariates", fake
    )
    bigwig = tmp_path / "a.bigwig"
    bigwig.touch()
    builder._load_chromatin(
        "x", tmp_path / "c.csv", [bigwig], "g", False
    )
    assert "windows" not in seen
    assert "transform" not in seen
