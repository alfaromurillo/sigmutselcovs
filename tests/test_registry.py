"""Tests for the project registry."""

import json

import pytest

from sigmutselcovs.registry import (
    EncodeChromatinSpec,
    GexpSpec,
    GtexSpec,
    ProjectSpec,
    RepliseqSpec,
    TrackRef,
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
        "ESCA",
        "GENERIC",
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


def test_esca_row():
    """ESCA has real ATAC and Roadmap coverage (E079 Esophagus) but
    no matching repliseq cell line anywhere (no esophageal/gastric/
    GI-tract line in ENCODE's Repli-seq set or GEO GSE137764) --
    see the row's description for the squamous/adenocarcinoma
    histology-mix caveat."""
    esca = get_project("ESCA")
    assert esca.gtex.representative_column == "gtex_esophagus_mucosa"
    assert (
        esca.atac.gdc_uuid == "fa18db6a-00fe-410f-982a-04d14d029812"
    )
    assert esca.roadmap.eids == ("E079",)
    assert esca.repliseq is None


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
    # No real-cohort row spells out roadmap marks (a project may
    # disable roadmap entirely with a null row, e.g. UCEC -- nothing
    # to merge). GENERIC is the deliberate exception: it adds DNase
    # to the standard 7 marks, since Roadmap has native DNase-seq
    # tracks for most EIDs via the same signal system.
    for code, row in raw["projects"].items():
        if row["roadmap"] is None or code == "GENERIC":
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


# --- encode_chromatin (new source, expanded-covariate-sources) ---


def test_gtex_spec_reduce_defaults_to_none():
    spec = GtexSpec(mapping_key="COAD", representative_column="x")
    assert spec.reduce is None


def test_encode_chromatin_spec_defaults():
    spec = EncodeChromatinSpec(
        tracks=(TrackRef(label="DNase", accession="ENCFF000AAA"),)
    )
    assert spec.assembly == "hg38"
    assert spec.tracks[0].label == "DNase"


def test_golden_row_with_encode_chromatin_and_generic_gtex(tmp_path):
    """A synthetic row exercising both new mechanisms together: a
    GENERIC-style tissue-agnostic gtex.reduce, plus encode_chromatin
    tracks -- no real cohort has encode_chromatin yet (Objective A/B
    of the covariate-expansion task adds those), so this is the
    registry-level contract test until then."""
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["gtex"]["reduce"] = "median"
    raw["projects"]["COAD"]["encode_chromatin"] = {
        "tracks": [
            {"label": "H3K23ac", "accession": "ENCFF000AAA"},
            {"label": "DNase", "accession": "ENCFF000BBB"},
        ]
    }
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    coad = load_registry(path)["COAD"]

    assert coad.gtex.reduce == "median"
    assert coad.encode_chromatin.assembly == "hg38"
    assert [t.label for t in coad.encode_chromatin.tracks] == [
        "H3K23ac",
        "DNase",
    ]
    assert coad.encode_chromatin.tracks[1].accession == "ENCFF000BBB"


def test_encode_chromatin_defaults_to_none():
    """SKCM/CESC are genuine negative results from the live
    ENCODE-portal search (Objective A.4, 2026-09-03): no adult
    melanocyte DNase-seq and no cervix DNase-seq at all exist there
    -- see each row's description. Every other cohort in the shipped
    registry now has at least a DNase encode_chromatin track."""
    assert get_project("SKCM").encode_chromatin is None
    assert get_project("CESC").encode_chromatin is None


def test_validate_rejects_bad_gtex_reduce(tmp_path):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["gtex"]["reduce"] = "mean"
    with pytest.raises(ValueError, match="gtex.reduce"):
        validate_registry(raw)


def test_validate_rejects_empty_encode_chromatin_tracks(tmp_path):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["encode_chromatin"] = {"tracks": []}
    with pytest.raises(ValueError, match="encode_chromatin"):
        validate_registry(raw)


def test_zero_roadmap_cohorts_have_real_encode_chromatin_rows():
    """ACC, BLCA, TGCT, UCEC had zero Roadmap epigenomes (Objective
    A.1 of the covariate-expansion task, 2026-09-03) -- ENCODE-native
    adult tissue coverage found via the portal search, real accessions
    (no placeholders). BLCA has no adult tissue DNase-seq available
    (only an embryonic experiment, rejected as a stage mismatch);
    the other three do."""
    expected_marks = {
        "ACC": {
            "H3K4me1",
            "H3K4me3",
            "H3K9me3",
            "H3K27me3",
            "H3K36me3",
            "DNase",
        },
        "BLCA": {"H3K4me1", "H3K27ac", "H3K4me3", "H3K36me3"},
        "TGCT": {
            "H3K4me1",
            "H3K4me3",
            "H3K27ac",
            "H3K27me3",
            "H3K9me3",
            "H3K36me3",
            "DNase",
        },
        "UCEC": {
            "H3K4me1",
            "H3K9me3",
            "H3K36me3",
            "H3K27ac",
            "H3K27me3",
            "DNase",
        },
    }
    for code, marks in expected_marks.items():
        spec = get_project(code)
        assert spec.roadmap is None, code  # still no Roadmap coverage
        assert spec.encode_chromatin is not None, code
        assert {
            t.label for t in spec.encode_chromatin.tracks
        } == marks
        assert all(
            t.accession.startswith("ENCFF")
            for t in spec.encode_chromatin.tracks
        ), code


def test_dnase_only_encode_chromatin_rows():
    """Objective A.3/A.4 (2026-09-03): DNase-seq added, on top of
    each cohort's existing Roadmap/ATAC coverage, for every cohort
    with a matched adult-tissue ENCODE DNase-seq experiment. OV/DLBC
    substitute it for their missing TCGA ATAC; the rest add it
    alongside their existing ATAC."""
    expected_dnase = {
        "COAD": "ENCFF013JSI",
        "BRCA": "ENCFF874CNE",
        "STAD": "ENCFF493HHP",
        "CHOL": "ENCFF972AOH",
        "ESCA": "ENCFF293DZU",
        "OV": "ENCFF596SZA",
        "DLBC": "ENCFF749AVF",
    }
    for code, accession in expected_dnase.items():
        spec = get_project(code)
        assert spec.encode_chromatin is not None, code
        tracks = {
            t.label: t.accession for t in spec.encode_chromatin.tracks
        }
        assert tracks == {"DNase": accession}, code


def test_generic_row():
    """GENERIC (Objective B, 2026-09-03, rebuilt 2026-09-04 with the
    disk-bounded streaming path) is a tissue-agnostic pool, not a
    real TCGA code -- gtex.reduce collapses every GTEx tissue into
    one column instead of picking a representative, atac/gexp/
    repliseq are null (TCGA-specific, no tissue-agnostic
    equivalent), and simple_matrix.gtex_column defaults correctly
    off representative_column without an explicit override.

    roadmap.eids is the full 127-epigenome Roadmap consolidated
    panel (E001-E129 excluding E060/E064, which genuinely don't
    exist -- verified live); roadmap.marks adds DNase to the
    standard 7 histone marks, since Roadmap has its own native
    DNase-seq tracks via the same signal system for most EIDs.
    encode_chromatin is 83 adult primary-tissue DNase-seq
    experiments from the ENCODE portal, one per unique biosample
    term, complementing Roadmap's own DNase coverage."""
    generic = get_project("GENERIC")
    assert generic.gtex.mapping_key == "GENERIC"
    assert generic.gtex.reduce == "median"
    assert generic.gtex.representative_column == (
        "gtex_pantissue_median"
    )
    assert generic.atac is None
    assert generic.gexp is None
    assert generic.repliseq is None
    assert generic.roadmap is not None
    assert len(generic.roadmap.eids) == 127
    assert len(set(generic.roadmap.eids)) == 127  # no duplicates
    assert "E060" not in generic.roadmap.eids
    assert "E064" not in generic.roadmap.eids
    assert set(generic.roadmap.marks) == {
        "H3K4me1",
        "H3K4me3",
        "H3K9ac",
        "H3K9me3",
        "H3K27ac",
        "H3K27me3",
        "H3K36me3",
        "DNase",
    }
    assert generic.encode_chromatin is not None
    assert len(generic.encode_chromatin.tracks) == 83
    assert (
        len({t.accession for t in generic.encode_chromatin.tracks})
        == 83
    )  # no duplicate accessions
    assert all(
        t.label == "DNase" for t in generic.encode_chromatin.tracks
    )
    assert generic.simple_matrix.gtex_column == (
        "gtex_pantissue_median"
    )


def test_generic_not_split_by_registry_project_convention():
    """tcga_analysis's REGISTRY_PROJECT = PROJECT.split('-')[0]
    subcohort logic (not sigmutselcovs code, but a constraint this
    registry entry must survive) is a no-op on a hyphen-free code."""
    assert "GENERIC".split("-")[0] == "GENERIC"  # noqa: SIM905
