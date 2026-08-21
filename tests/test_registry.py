"""Tests for the project registry."""

import json

import pytest

from sigmutselcovs.registry import (
    GexpSpec,
    ProjectSpec,
    RepliseqSpec,
    available_projects,
    get_project,
    load_registry,
    location_projects_registry,
    validate_registry,
)


def test_available_projects():
    assert available_projects() == [
        "ACC",
        "BLCA",
        "BRCA",
        "CESC",
        "CHOL",
        "COAD",
        "DLBC",
        "OV",
        "SKCM",
        "STAD",
        "TGCT",
        "UCEC",
    ]


def test_get_project_unknown_lists_available():
    with pytest.raises(ValueError, match="BRCA, CESC, CHOL, COAD"):
        get_project("LUAD")


def test_get_project_case_insensitive():
    assert get_project("coad").code == "COAD"


def test_coad_golden_row():
    """The COAD row must reproduce the historical setup exactly."""
    coad = get_project("COAD")
    assert isinstance(coad, ProjectSpec)
    assert coad.gtex.mapping_key == "COAD"
    assert (
        coad.gtex.representative_column
        == "gtex_colon_transverse_mucosa"
    )
    assert coad.gexp.tcga_project_id == "TCGA-COAD"
    assert coad.gexp.workflow_type == "STAR - Counts"
    assert coad.gexp.metrics is None  # all six STAR metrics
    assert (
        coad.atac.gdc_uuid == "26b96cd9-dce0-4340-b15f-9e0afbb6312c"
    )
    assert coad.atac.column_prefix == "coad"
    assert coad.roadmap.eids == ("E075", "E106", "E101", "E102")
    assert coad.roadmap.marks == (
        "H3K4me1",
        "H3K4me3",
        "H3K9ac",
        "H3K9me3",
        "H3K27ac",
        "H3K27me3",
        "H3K36me3",
    )
    assert coad.repliseq.type == "mat"
    assert coad.repliseq.cell_line == "HCT116"
    assert coad.repliseq.assembly == "hg38"
    assert coad.repliseq.filename == (
        "GSE137764_HCT_GaussiansGSE137764_mooth_scaled_autosome.mat"
    )
    assert coad.repliseq.n_fractions == 16
    assert (
        coad.simple_matrix.gtex_column
        == "gtex_colon_transverse_mucosa"
    )
    assert coad.simple_matrix.body_pattern == "h3k4me3_fc_signal_body"
    assert coad.simple_matrix.promoter_pattern == (
        "h3k9me3_fc_signal_promoter"
    )


def test_brca_row():
    brca = get_project("BRCA")
    assert (
        brca.gtex.representative_column
        == "gtex_breast_mammary_tissue"
    )
    assert (
        brca.atac.gdc_uuid == "f1c06cd3-cf35-41cc-bc75-6db273c94273"
    )
    assert brca.roadmap.eids == ("E027", "E028", "E119")
    assert brca.repliseq.type == "fraction_bigwigs"
    assert brca.repliseq.assembly == "hg19"
    assert brca.repliseq.cell_line == "MCF-7"
    assert [t.label for t in brca.repliseq.tracks] == [
        "g1b",
        "s1",
        "s2",
        "s3",
        "s4",
        "g2",
    ]
    assert brca.repliseq.tracks[0].accession == "ENCFF001GSV"
    assert brca.repliseq.mrt_fraction_cols is None  # all fractions
    # simple_matrix.gtex_column defaults to the representative column
    assert (
        brca.simple_matrix.gtex_column == "gtex_breast_mammary_tissue"
    )


def test_chol_row():
    """CHOL has no bile-duct-specific epigenome or cell line; both
    Roadmap and repliseq fall back to liver (organ-level proxy, not
    a cholangiocyte match -- see the row's description)."""
    chol = get_project("CHOL")
    assert chol.gtex.representative_column == "gtex_liver"
    assert (
        chol.atac.gdc_uuid == "78b6baeb-8b9b-486b-ae07-d9094cadaaa2"
    )
    assert chol.roadmap.eids == ("E066", "E118")
    assert chol.repliseq.type == "fraction_bigwigs"
    assert chol.repliseq.assembly == "hg19"
    assert chol.repliseq.cell_line == "HepG2"
    assert [t.label for t in chol.repliseq.tracks] == [
        "g1b",
        "s1",
        "s2",
        "s3",
        "s4",
        "g2",
    ]
    assert chol.repliseq.tracks[0].accession == "ENCFF001GPC"


def test_dlbc_row():
    """DLBC has no lymph-node GTEx column (falls back to
    Whole_Blood/Spleen) and no TCGA ATAC-seq coverage; Roadmap uses
    the direct cell-of-origin match (E032 primary B cells) plus a
    broader-coverage epigenome (E062 PBMCs), and repliseq uses
    GM12878, the standard lymphoblastoid line -- see the row's
    description."""
    dlbc = get_project("DLBC")
    assert dlbc.gtex.representative_column == "gtex_whole_blood"
    assert dlbc.atac is None
    assert dlbc.roadmap.eids == ("E032", "E062")
    assert dlbc.repliseq.type == "fraction_bigwigs"
    assert dlbc.repliseq.assembly == "hg19"
    assert dlbc.repliseq.cell_line == "GM12878"
    assert [t.label for t in dlbc.repliseq.tracks] == [
        "g1b",
        "s1",
        "s2",
        "s3",
        "s4",
        "g2",
    ]
    assert dlbc.repliseq.tracks[0].accession == "ENCFF001GNK"


def test_ucec_row():
    """UCEC has no matching Roadmap epigenome or repliseq cell line;
    both sources are disabled rather than using a mismatched proxy."""
    ucec = get_project("UCEC")
    assert ucec.gtex.mapping_key == "UCEC"
    assert ucec.gtex.representative_column == "gtex_uterus"
    assert ucec.gexp.tcga_project_id == "TCGA-UCEC"
    assert (
        ucec.atac.gdc_uuid == "e447195b-eeb9-459f-83db-9c748aafe395"
    )
    assert ucec.atac.column_prefix == "ucec"
    assert ucec.roadmap is None
    assert ucec.repliseq is None
    assert ucec.simple_matrix.gtex_column == "gtex_uterus"


def test_stad_row():
    """STAD has real Roadmap coverage (Gastric + Stomach Mucosa) but
    no matching repliseq cell line in either ENCODE cohort or GEO
    GSE137764."""
    stad = get_project("STAD")
    assert stad.gtex.mapping_key == "STAD"
    assert stad.gtex.representative_column == "gtex_stomach_mucosa"
    assert stad.gexp.tcga_project_id == "TCGA-STAD"
    assert (
        stad.atac.gdc_uuid == "a59a7128-62ff-44f9-8f35-497fda8b0beb"
    )
    assert stad.atac.column_prefix == "stad"
    assert stad.roadmap.eids == ("E094", "E110")
    assert stad.repliseq is None
    assert stad.simple_matrix.gtex_column == "gtex_stomach_mucosa"


def test_tgct_row():
    """TGCT has ATAC and GTEx coverage but no matching Roadmap
    epigenome (no testis anatomy in the standard 127-panel) and no
    matching repliseq cell line."""
    tgct = get_project("TGCT")
    assert tgct.gtex.mapping_key == "TGCT"
    assert tgct.gtex.representative_column == "gtex_testis"
    assert tgct.gexp.tcga_project_id == "TCGA-TGCT"
    assert (
        tgct.atac.gdc_uuid == "2f27b5c3-73ce-4ddc-9c04-eb6187875788"
    )
    assert tgct.atac.column_prefix == "tgct"
    assert tgct.roadmap is None
    assert tgct.repliseq is None
    assert tgct.simple_matrix.gtex_column == "gtex_testis"


def test_ov_row():
    """OV has no TCGA ATAC-seq coverage (Corces et al. 2018 didn't
    profile ovarian cancer) and no matching repliseq cell line (no
    ovarian line in ENCODE's Repli-seq experiment set or GEO
    GSE137764), but does have a real Roadmap match (E097 Ovary)."""
    ov = get_project("OV")
    assert ov.gtex.mapping_key == "OV"
    assert ov.gtex.representative_column == "gtex_ovary"
    assert ov.gexp.tcga_project_id == "TCGA-OV"
    assert ov.atac is None
    assert ov.roadmap.eids == ("E097",)
    assert ov.repliseq is None
    assert ov.simple_matrix.gtex_column == "gtex_ovary"


def test_skcm_row():
    """SKCM has real ATAC and Roadmap coverage (primary melanocyte
    cultures E059/E061, melanoma's cell of origin) but no matching
    repliseq cell line (no melanoma/melanocyte line in ENCODE's
    Repli-seq experiment set or GEO GSE137764, which remains
    H1/H9/HCT116-only)."""
    skcm = get_project("SKCM")
    assert skcm.gtex.mapping_key == "SKCM"
    assert (
        skcm.gtex.representative_column
        == "gtex_skin_sun_exposed_lower_leg"
    )
    assert skcm.gexp.tcga_project_id == "TCGA-SKCM"
    assert (
        skcm.atac.gdc_uuid == "9b87d207-8c8e-47da-b56b-333bcf85856d"
    )
    assert skcm.atac.column_prefix == "skcm"
    assert skcm.roadmap.eids == ("E059", "E061")
    assert skcm.repliseq is None
    assert (
        skcm.simple_matrix.gtex_column
        == "gtex_skin_sun_exposed_lower_leg"
    )


def test_defaults_merge():
    """Roadmap marks and gexp defaults come from the defaults block."""
    raw = json.loads(location_projects_registry.read_text())
    # no project row spells out roadmap marks (a project may disable
    # roadmap entirely with a null row, e.g. UCEC -- nothing to merge)
    for row in raw["projects"].values():
        if row["roadmap"] is None:
            continue
        assert "marks" not in row["roadmap"]
    brca = get_project("BRCA")
    assert brca.roadmap.marks == tuple(
        raw["defaults"]["roadmap"]["marks"]
    )
    assert brca.gexp.data_type == "Gene Expression Quantification"
    assert brca.gexp.tissue_types == ("Tumor", "Normal")


def test_null_source_disables_block(tmp_path):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["atac"] = None
    raw["projects"]["COAD"]["repliseq"] = None
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    coad = load_registry(path)["COAD"]
    assert coad.atac is None
    assert coad.repliseq is None
    assert coad.roadmap is not None


def test_validate_rejects_bad_schema_version():
    with pytest.raises(ValueError, match="schema_version"):
        validate_registry(
            {"schema_version": 2, "projects": {"X": {}}}
        )


def test_validate_rejects_unknown_mapping_key(tmp_path):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["gtex"]["mapping_key"] = "NOPE"
    with pytest.raises(ValueError, match="mapping_key"):
        validate_registry(raw)


def test_validate_rejects_missing_gtex():
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["gtex"] = None
    with pytest.raises(ValueError, match="gtex source is required"):
        validate_registry(raw)


def test_validate_rejects_bad_repliseq():
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["repliseq"]["type"] = "bogus"
    with pytest.raises(ValueError, match="repliseq.type"):
        validate_registry(raw)
    raw["projects"]["COAD"]["repliseq"] = {
        "type": "mat",
        "assembly": "hg38",
        "cell_line": "HCT116",
    }
    with pytest.raises(ValueError, match="needs a filename"):
        validate_registry(raw)


def test_unknown_spec_key_rejected(tmp_path):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["atac"]["bogus_key"] = 1
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="bogus_key"):
        load_registry(path)


def test_specs_are_frozen():
    coad = get_project("COAD")
    with pytest.raises(AttributeError):
        coad.gexp.tcga_project_id = "TCGA-LUAD"


def test_gexp_spec_defaults():
    spec = GexpSpec(tcga_project_id="TCGA-X")
    assert spec.metrics is None
    assert spec.tissue_types == ("Tumor", "Normal")


def test_repliseq_bin_size_default():
    assert (
        RepliseqSpec(
            type="mat", assembly="hg38", cell_line="X", filename="f"
        ).bin_size
        == 50_000
    )
