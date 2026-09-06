"""Builder tests with stubbed loaders (no real data, no network)."""

import logging

import numpy as np
import pandas as pd
import pytest

from sigmutselcovs import builder
from sigmutselcovs.builder import (
    CovariateMatrices,
    build_covariate_matrix,
)
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import RepliseqSpec

GENES = [f"ENSG{i:011d}" for i in range(1, 6)]
GTFS = {"hg19": "gencode19.gtf", "hg38": "gencode38.gtf"}


def _series(name, offset=0.0):
    return pd.Series(
        np.arange(len(GENES), dtype=float) + offset,
        index=pd.Index(GENES, name="ensembl_gene_id"),
        name=name,
    )


def _frame(columns, offset=0.0):
    return pd.DataFrame(
        {
            c: np.arange(len(GENES), dtype=float) + offset + i
            for i, c in enumerate(columns)
        },
        index=pd.Index(GENES, name="ensembl_gene_id"),
    )


@pytest.fixture
def stubbed(tmp_path, monkeypatch):
    """Stub every loader; create the files the skip logic checks."""
    paths = project_paths(tmp_path)
    paths.rt_dir.mkdir(parents=True)
    (
        paths.rt_dir
        / "GSE137764_HCT_GaussiansGSE137764_mooth_scaled_autosome.mat"
    ).touch()
    paths.roadmap_dir.mkdir(parents=True)
    (paths.roadmap_dir / "E075-H3K4me3.fc.signal.bigwig").touch()
    (paths.roadmap_dir / "E075-H3K9me3.fc.signal.bigwig").touch()
    paths.atac_dir.mkdir(parents=True)
    (paths.atac_dir / "COAD_ABC_T1.insertions.bw").touch()

    gtex = _frame(
        ["gtex_colon_sigmoid", "gtex_colon_transverse_mucosa"]
    )
    per_sample = _frame(
        [
            "TCGA-AA-0001-01A_tpm_unstranded",
            "TCGA-AA-0002-01A_tpm_unstranded",
        ],
        offset=10,
    )
    rt_fractions = _frame(
        [f"rt_s{i}" for i in range(1, 17)], offset=0.1
    )
    atac = _frame(
        [
            "coad_abc_t1_insertions_body",
            "coad_abc_t1_insertions_promoter",
        ],
        offset=20,
    )
    roadmap = _frame(
        [
            "e075_h3k4me3_fc_signal_body",
            "e075_h3k9me3_fc_signal_promoter",
        ],
        offset=30,
    )

    def fake_mean_gexp(
        location_csv,
        tcga_dir,
        *,
        tissue_type="Tumor",
        force_generation=False,
        **kw,
    ):
        return _series(
            "tpm_unstranded",
            offset=100 if tissue_type == "Tumor" else 200,
        )

    monkeypatch.setattr(
        builder, "import_gtex", lambda *a, **k: gtex.copy()
    )
    monkeypatch.setattr(
        builder, "load_or_generate_mean_tcga_gexp", fake_mean_gexp
    )
    monkeypatch.setattr(
        builder,
        "load_or_generate_tcga_gexp_per_sample",
        lambda **k: per_sample.copy(),
    )
    monkeypatch.setattr(
        builder,
        "load_or_generate_mrt",
        lambda *a, **k: _series("mrt", offset=0.2),
    )
    monkeypatch.setattr(
        builder,
        "load_or_generate_rt_fractions",
        lambda *a, **k: rt_fractions.copy(),
    )

    def fake_chromatin(location_df, tracks, gtf_path, **kw):
        if "roadmap" in str(location_df):
            return roadmap.copy()
        return atac.copy()

    monkeypatch.setattr(
        builder,
        "load_or_generate_chromatin_covariates",
        fake_chromatin,
    )
    return tmp_path


def test_full_matrix_column_order(stubbed):
    """gtex, tumor, normal, per-sample, mrt, clr, ATAC, Roadmap."""
    matrices = build_covariate_matrix(
        "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
    )
    cols = list(matrices.full.columns)
    expected = (
        [
            "gtex_colon_sigmoid",
            "gtex_colon_transverse_mucosa",
            "tpm_unstranded",
            "tpm_unstranded_normal",
            "TCGA-AA-0001-01A_tpm_unstranded",
            "TCGA-AA-0002-01A_tpm_unstranded",
            "mrt",
        ]
        + [f"clr_rt_s{i}" for i in range(1, 17)]
        + [
            "coad_abc_t1_insertions_body",
            "coad_abc_t1_insertions_promoter",
            "e075_h3k4me3_fc_signal_body",
            "e075_h3k9me3_fc_signal_promoter",
        ]
    )
    assert cols == expected


def test_unpacks_as_three_tuple(stubbed):
    full, simple, tcga = build_covariate_matrix(
        "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
    )
    assert isinstance(full, pd.DataFrame)
    assert list(simple.columns) == [
        "gtex_colon_transverse_mucosa",
        "mrt",
        "h3k4me3_body",
        "h3k9me3_prom",
    ]
    assert list(tcga.columns) == [
        "gexp",
        "atac_body",
        "atac_promoter",
    ]


def test_simple_and_tcga_values(stubbed):
    matrices = build_covariate_matrix(
        "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
    )
    full = matrices.full
    pd.testing.assert_series_equal(
        matrices.simple["h3k4me3_body"],
        full[["e075_h3k4me3_fc_signal_body"]].mean(axis=1),
        check_names=False,
    )
    pd.testing.assert_series_equal(
        matrices.tcga["atac_body"],
        full[["coad_abc_t1_insertions_body"]].mean(axis=1),
        check_names=False,
    )
    pd.testing.assert_series_equal(
        matrices.tcga["gexp"],
        full["tpm_unstranded"],
        check_names=False,
    )


def test_exclude_sources(stubbed, caplog):
    with caplog.at_level(logging.WARNING):
        matrices = build_covariate_matrix(
            "COAD",
            stubbed,
            gencode_gtfs=GTFS,
            apply_fixes=False,
            exclude=("atac", "repliseq"),
        )
    cols = matrices.full.columns
    assert not any(c.startswith("coad_") for c in cols)
    assert "mrt" not in cols
    # simple/tcga degrade with warnings instead of KeyError
    assert "mrt" not in matrices.simple.columns
    assert list(matrices.tcga.columns) == ["gexp"]
    assert any("no mrt column" in m for m in caplog.messages)


def test_include_subset(stubbed):
    matrices = build_covariate_matrix(
        "COAD",
        stubbed,
        gencode_gtfs=GTFS,
        apply_fixes=False,
        include=("gtex",),
    )
    assert list(matrices.full.columns) == [
        "gtex_colon_sigmoid",
        "gtex_colon_transverse_mucosa",
    ]


def test_unknown_source_rejected(stubbed):
    with pytest.raises(ValueError, match="Unknown sources"):
        build_covariate_matrix(
            "COAD", stubbed, gencode_gtfs=GTFS, include=("bogus",)
        )


def test_missing_repliseq_data_skips(stubbed):
    paths = project_paths(stubbed)
    (
        paths.rt_dir
        / "GSE137764_HCT_GaussiansGSE137764_mooth_scaled_autosome.mat"
    ).unlink()
    matrices = build_covariate_matrix(
        "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
    )
    assert "mrt" not in matrices.full.columns


def test_clr_gate_two_fractions(stubbed, monkeypatch, caplog):
    two = _frame(["rt_s1", "rt_s2"], offset=0.1)
    monkeypatch.setattr(
        builder,
        "load_or_generate_rt_fractions",
        lambda *a, **k: two.copy(),
    )
    spec = RepliseqSpec(
        type="mat",
        assembly="hg38",
        cell_line="X",
        filename=(
            "GSE137764_HCT_GaussiansGSE137764_mooth_"
            "scaled_autosome.mat"
        ),
        include_clr_fractions=None,
    )
    with caplog.at_level(logging.INFO):
        frames = builder._load_repliseq(
            spec, project_paths(stubbed), "gencode38.gtf", False
        )
    assert len(frames) == 1  # mrt only, CLR gated off
    assert any("Skipping CLR" in m for m in caplog.messages)


def test_duplicate_columns_warned(stubbed, monkeypatch, caplog):
    # per-sample block re-uses a gtex column name
    dup = _frame(["gtex_colon_sigmoid"], offset=10)
    monkeypatch.setattr(
        builder,
        "load_or_generate_tcga_gexp_per_sample",
        lambda **k: dup.copy(),
    )
    with caplog.at_level(logging.WARNING):
        build_covariate_matrix(
            "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
        )
    assert any(
        "Duplicated covariate columns" in m for m in caplog.messages
    )


def test_cache_matrices_and_manifest(stubbed):
    matrices = build_covariate_matrix(
        "COAD",
        stubbed,
        gencode_gtfs=GTFS,
        apply_fixes=False,
        cache_matrices=True,
    )
    paths = project_paths(stubbed)
    assert paths.matrix_full_parquet.exists()
    assert paths.matrix_simple_csv.exists()
    assert paths.matrix_tcga_csv.exists()
    import json

    manifest = json.loads(paths.build_manifest_json.read_text())
    assert manifest["project"] == "COAD"
    assert manifest["shapes"]["full"] == list(matrices.full.shape)
    assert set(manifest["sources_included"]) == {
        "gtex",
        "gexp",
        "repliseq",
        "atac",
        "roadmap",
    }


def test_column_dictionary(stubbed):
    matrices = build_covariate_matrix(
        "COAD",
        stubbed,
        gencode_gtfs=GTFS,
        apply_fixes=False,
        cache_matrices=True,
    )
    paths = project_paths(stubbed)
    dictionary = pd.read_csv(paths.column_dictionary_csv)
    # one row per column of the built matrix, same order
    assert dictionary["column"].tolist() == list(
        matrices.full.columns
    )
    by_col = dictionary.set_index("column")
    assert by_col.loc["gtex_colon_sigmoid", "source"] == "gtex"
    assert by_col.loc["tpm_unstranded", "source"] == "gexp_mean"
    assert (
        by_col.loc["TCGA-AA-0001-01A_tpm_unstranded", "source"]
        == "gexp_per_sample"
    )
    assert by_col.loc["mrt", "source"] == "mrt"
    assert by_col.loc["clr_rt_s1", "source"] == "clr"
    assert by_col.loc["clr_rt_s7", "detail"] == "S-phase fraction 7"
    atac_row = by_col.loc["coad_abc_t1_insertions_body"]
    assert atac_row["source"] == "atac"
    assert atac_row["detail"] == "gene body [TSS+200, TES]"
    assert atac_row["assembly"] == "hg38"
    roadmap_row = by_col.loc["e075_h3k9me3_fc_signal_promoter"]
    assert roadmap_row["source"] == "roadmap"
    assert roadmap_row["detail"] == "promoter [TSS-2000, TSS+200]"
    assert roadmap_row["assembly"] == "hg19"
    assert by_col.loc["mrt", "cell_line_or_epigenome"] == "HCT116"
    # no fixes applied -> every transform is none
    assert (dictionary["transform"] == "none").all()


def test_column_dictionary_partial_build(stubbed):
    build_covariate_matrix(
        "COAD",
        stubbed,
        gencode_gtfs=GTFS,
        apply_fixes=False,
        cache_matrices=True,
        include=("gtex", "repliseq"),
    )
    dictionary = pd.read_csv(
        project_paths(stubbed).column_dictionary_csv
    )
    assert set(dictionary["source"]) == {"gtex", "mrt", "clr"}


def test_column_dictionary_records_transforms(stubbed, monkeypatch):
    # a strongly right-skewed gtex column triggers fix_skewness
    import numpy as np

    skewed = pd.DataFrame(
        {
            "gtex_colon_sigmoid": [0.0, 0.1, 0.2, 0.1, 1000.0],
            "gtex_colon_transverse_mucosa": np.linspace(1, 2, 5),
        },
        index=pd.Index(GENES, name="ensembl_gene_id"),
    )
    monkeypatch.setattr(
        builder, "import_gtex", lambda *a, **k: skewed.copy()
    )
    build_covariate_matrix(
        "COAD",
        stubbed,
        gencode_gtfs=GTFS,
        apply_fixes=True,
        cache_matrices=True,
        include=("gtex",),
    )
    dictionary = pd.read_csv(
        project_paths(stubbed).column_dictionary_csv
    )
    by_col = dictionary.set_index("column")
    assert by_col.loc["gtex_colon_sigmoid", "transform"].startswith(
        "log("
    )


def test_named_tuple_type(stubbed):
    matrices = build_covariate_matrix(
        "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
    )
    assert isinstance(matrices, CovariateMatrices)


def test_fraction_bigwigs_skips_when_missing(stubbed, caplog):
    from sigmutselcovs.registry import TrackRef

    spec = RepliseqSpec(
        type="fraction_bigwigs",
        assembly="hg19",
        cell_line="MCF-7",
        tracks=(
            TrackRef(label="s1", accession="ENCFF000AAA"),
            TrackRef(label="s2", accession="ENCFF000BBB"),
        ),
    )
    with caplog.at_level(logging.WARNING):
        frames = builder._load_repliseq(
            spec, project_paths(stubbed), "gencode19.gtf", False
        )
    assert frames == []
    assert any("missing bigWigs" in m for m in caplog.messages)


def test_wavelet_uses_cached_csv(stubbed, monkeypatch):
    from sigmutselcovs.registry import TrackRef

    called = {}

    def fake_wavelet(location_csv, bigwig, gtf, **kw):
        called["args"] = (location_csv, bigwig)
        return _series("rt_wavelet")

    monkeypatch.setattr(
        builder, "load_or_generate_rt_wavelet", fake_wavelet
    )
    paths = project_paths(stubbed)
    paths.rt_wavelet_csv.parent.mkdir(parents=True, exist_ok=True)
    paths.rt_wavelet_csv.touch()  # cache present, bigWig absent
    spec = RepliseqSpec(
        type="wavelet",
        assembly="hg19",
        cell_line="LNCaP",
        tracks=(TrackRef(label="el", accession="ENCFF000CCC"),),
    )
    frames = builder._load_repliseq(
        spec, paths, "gencode19.gtf", False
    )
    assert len(frames) == 1
    assert list(frames[0].columns) == ["rt_wavelet"]
    assert called["args"][0] == paths.rt_wavelet_csv


# --- encode_chromatin (new source, expanded-covariate-sources) ---


def _registry_with_encode_chromatin(tmp_path):
    """A COAD row extended with encode_chromatin tracks, written to
    its own registry JSON (mirrors test_registry.py's tmp_path
    pattern) -- no shipped cohort has this source yet."""
    import json

    from sigmutselcovs.registry import location_projects_registry

    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["encode_chromatin"] = {
        "tracks": [
            {"label": "H3K23ac", "accession": "ENCFF000DDD"},
            {"label": "DNase", "accession": "ENCFF000EEE"},
        ]
    }
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    return path


def test_encode_chromatin_block_appends_after_roadmap(
    stubbed, monkeypatch
):
    """encode_chromatin is skipped for the shipped COAD row (no
    tracks registered); with a synthetic row that has tracks, its
    columns land last, after Roadmap -- same append-only ordering
    documented in build_covariate_matrix's concatenation comment."""
    registry_path = _registry_with_encode_chromatin(stubbed)
    paths = project_paths(stubbed)
    paths.encode_chromatin_dir.mkdir(parents=True)
    (
        paths.encode_chromatin_dir
        / "H3K23ac_fc_signal_ENCFF000DDD.bigWig"
    ).touch()
    (
        paths.encode_chromatin_dir
        / "DNase_fc_signal_ENCFF000EEE.bigWig"
    ).touch()

    encode_chromatin = _frame(
        [
            "h3k23ac_fc_signal_encff000ddd_body",
            "dnase_fc_signal_encff000eee_body",
        ],
        offset=40,
    )

    def fake_chromatin(location_df, tracks, gtf_path, **kw):
        if "roadmap" in str(location_df):
            return _frame(
                [
                    "e075_h3k4me3_fc_signal_body",
                    "e075_h3k9me3_fc_signal_promoter",
                ],
                offset=30,
            )
        if "encode" in str(location_df):
            return encode_chromatin.copy()
        return _frame(
            [
                "coad_abc_t1_insertions_body",
                "coad_abc_t1_insertions_promoter",
            ],
            offset=20,
        )

    monkeypatch.setattr(
        builder,
        "load_or_generate_chromatin_covariates",
        fake_chromatin,
    )

    matrices = build_covariate_matrix(
        "COAD",
        stubbed,
        gencode_gtfs=GTFS,
        apply_fixes=False,
        registry_path=registry_path,
    )
    assert list(matrices.full.columns)[-2:] == [
        "h3k23ac_fc_signal_encff000ddd_body",
        "dnase_fc_signal_encff000eee_body",
    ]


def test_encode_chromatin_skipped_when_not_in_registry(stubbed):
    """The shipped COAD row has no encode_chromatin -- unchanged
    behavior, no new columns, no error."""
    matrices = build_covariate_matrix(
        "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
    )
    assert not any(
        "fc_signal" in c and "encff" in c.lower()
        for c in matrices.full.columns
    )


def test_encode_chromatin_skipped_when_files_missing(stubbed, caplog):
    """Registry entry present but no bigWigs downloaded yet and no
    cache -- skip with a warning, same as any other chromatin source
    (mirrors test_missing_repliseq_data_skips)."""
    registry_path = _registry_with_encode_chromatin(stubbed)
    with caplog.at_level(logging.WARNING):
        matrices = build_covariate_matrix(
            "COAD",
            stubbed,
            gencode_gtfs=GTFS,
            apply_fixes=False,
            registry_path=registry_path,
        )
    assert not any(
        "encff" in c.lower() for c in matrices.full.columns
    )
    assert any(
        "Skipping encode_chromatin" in m for m in caplog.messages
    )


# --- average_by_assay (pooled collapsing across chromatin sources) ---


def test_average_by_assay_pools_roadmap_and_encode_together(
    stubbed, monkeypatch
):
    """average_by_assay=True must combine roadmap+encode_chromatin's
    raw per-track pools BEFORE collapsing, not collapse each source
    separately -- otherwise two sources sharing a mark (e.g. both
    having a "dnase" track) would each produce their own
    "dnase_body" column, colliding on the same name after concat
    instead of averaging into one shared column (the exact bug
    average_by_assay=True was built to avoid for GENERIC's native-
    Roadmap DNase + ENCODE's 83 DNase tracks)."""
    import json

    from sigmutselcovs.registry import location_projects_registry

    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["average_by_assay"] = True
    raw["projects"]["COAD"]["encode_chromatin"] = {
        "tracks": [{"label": "DNase", "accession": "ENCFF000EEE"}]
    }
    registry_path = stubbed / "projects.json"
    registry_path.write_text(json.dumps(raw))

    paths = project_paths(stubbed)
    paths.encode_chromatin_dir.mkdir(parents=True)
    (
        paths.encode_chromatin_dir
        / "DNase_fc_signal_ENCFF000EEE.bigWig"
    ).touch()

    roadmap_frame = _frame(
        ["e075_h3k4me3_fc_signal_body", "e075_dnase_body"], offset=30
    )
    encode_frame = _frame(
        ["dnase_fc_signal_encff000eee_body"], offset=40
    )

    def fake_chromatin(location_df, tracks, gtf_path, **kw):
        if "roadmap" in str(location_df):
            return roadmap_frame.copy()
        if "encode" in str(location_df):
            return encode_frame.copy()
        return _frame(
            [
                "coad_abc_t1_insertions_body",
                "coad_abc_t1_insertions_promoter",
            ],
            offset=20,
        )

    monkeypatch.setattr(
        builder,
        "load_or_generate_chromatin_covariates",
        fake_chromatin,
    )

    matrices = build_covariate_matrix(
        "COAD",
        stubbed,
        gencode_gtfs=GTFS,
        apply_fixes=False,
        exclude=("atac",),
        registry_path=registry_path,
    )
    full = matrices.full

    assert list(full.columns).count("dnase_body") == 1
    assert "h3k4me3_body" in full.columns
    assert "e075_dnase_body" not in full.columns
    assert "dnase_fc_signal_encff000eee_body" not in full.columns

    expected_dnase = (
        roadmap_frame["e075_dnase_body"]
        + encode_frame["dnase_fc_signal_encff000eee_body"]
    ) / 2
    pd.testing.assert_series_equal(
        full["dnase_body"], expected_dnase.rename("dnase_body")
    )


def test_average_by_assay_false_keeps_per_track_columns(stubbed):
    """Default behavior (the shipped COAD row) is unaffected -- no
    collapsing, per-track column names survive untouched."""
    matrices = build_covariate_matrix(
        "COAD", stubbed, gencode_gtfs=GTFS, apply_fixes=False
    )
    assert "e075_h3k4me3_fc_signal_body" in matrices.full.columns
    assert "dnase_body" not in matrices.full.columns


# --- combine_with_generic ---


def test_combine_with_generic_outer_joins_on_gene_index():
    cohort = pd.DataFrame(
        {"cohort_col": [1.0, 2.0]},
        index=pd.Index(["ENSG_A", "ENSG_B"], name="ensembl_gene_id"),
    )
    generic = pd.DataFrame(
        {"generic_col": [10.0, 20.0]},
        index=pd.Index(["ENSG_B", "ENSG_C"], name="ensembl_gene_id"),
    )
    combined = builder.combine_with_generic(cohort, generic)

    assert list(combined.columns) == ["cohort_col", "generic_col"]
    assert set(combined.index) == {"ENSG_A", "ENSG_B", "ENSG_C"}
    # ENSG_A has no generic data, ENSG_C has no cohort data -- outer
    # join means NaN, not a dropped row.
    assert pd.isna(combined.loc["ENSG_A", "generic_col"])
    assert pd.isna(combined.loc["ENSG_C", "cohort_col"])
    assert combined.loc["ENSG_B", "cohort_col"] == 2.0
    assert combined.loc["ENSG_B", "generic_col"] == 10.0
