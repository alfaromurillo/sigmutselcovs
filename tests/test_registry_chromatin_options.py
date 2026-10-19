"""Registry: the optional chromatin windows, transform and peaks."""

import json

import pytest

from sigmutselcovs.registry import (
    EncodeChromatinSpec,
    PeakRef,
    RoadmapSpec,
    TrackRef,
    load_registry,
    location_projects_registry,
)


def _registry(tmp_path, edit):
    raw = json.loads(location_projects_registry.read_text())
    edit(raw)
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    return path


def test_new_fields_are_optional():
    spec = RoadmapSpec(eids=("E075",), marks=("H3K4me3",))
    assert spec.windows == ("body", "promoter")
    assert spec.transform is None
    assert spec.peaks == ()
    assert spec.peaks_url_template.format(
        eid="E075", mark="H3K4me3", peak_type="narrowPeak"
    ) == (
        "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/"
        "consolidated/narrowPeak/E075-H3K4me3.narrowPeak.gz"
    )
    encode = EncodeChromatinSpec(
        tracks=(TrackRef(label="DNase", accession="ENCFF000AAA"),)
    )
    assert encode.windows == ("body", "promoter")
    assert encode.transform is None
    assert encode.peak_tracks == ()


def test_shipped_rows_keep_the_defaults():
    for code, spec in load_registry().items():
        for source in (spec.roadmap, spec.encode_chromatin):
            if source is None:
                continue
            assert source.windows == ("body", "promoter"), code
            assert source.transform is None, code
        if spec.roadmap is not None:
            assert spec.roadmap.peaks == (), code
        if spec.encode_chromatin is not None:
            assert spec.encode_chromatin.peak_tracks == (), code


def test_new_fields_parse(tmp_path):
    def edit(raw):
        coad = raw["projects"]["COAD"]
        coad["roadmap"].update(
            windows=["cds_span", "cds_span2k"],
            transform="log1p",
            peaks=["narrowPeak", "broadPeak"],
        )
        coad["encode_chromatin"]["windows"] = ["body"]
        coad["encode_chromatin"]["peak_tracks"] = [
            {"label": "H3K4me3", "accession": "ENCFF000PKA"},
            {
                "label": "H3K27ac",
                "accession": "ENCFF000PKB",
                "peak_type": "broadPeak",
            },
        ]

    coad = load_registry(_registry(tmp_path, edit))["COAD"]
    assert coad.roadmap.windows == ("cds_span", "cds_span2k")
    assert coad.roadmap.transform == "log1p"
    assert coad.roadmap.peaks == ("narrowPeak", "broadPeak")
    assert coad.roadmap.marks  # defaults still merged under
    assert coad.encode_chromatin.windows == ("body",)
    assert coad.encode_chromatin.transform is None
    assert coad.encode_chromatin.peak_tracks == (
        PeakRef(label="H3K4me3", accession="ENCFF000PKA"),
        PeakRef(
            label="H3K27ac",
            accession="ENCFF000PKB",
            peak_type="broadPeak",
        ),
    )


def test_defaults_block_options_apply(tmp_path):
    def edit(raw):
        raw["defaults"]["roadmap"]["windows"] = ["cds_span"]

    registry = load_registry(_registry(tmp_path, edit))
    assert registry["BRCA"].roadmap.windows == ("cds_span",)


@pytest.mark.parametrize(
    "source, field, value, match",
    [
        ("roadmap", "windows", ["body", "exon"], "unknown names"),
        ("roadmap", "windows", [], "must not be empty"),
        ("roadmap", "windows", ["body", "body"], "repeats"),
        ("roadmap", "windows", "body", "must be a list"),
        ("encode_chromatin", "windows", ["gene"], "unknown names"),
        ("roadmap", "transform", "log2", "transform"),
        ("encode_chromatin", "transform", "sqrt", "transform"),
        ("roadmap", "peaks", ["gappedPeak"], "unknown names"),
        (
            "encode_chromatin",
            "peak_tracks",
            [{"label": "H3K4me3"}],
            "label and an accession",
        ),
        (
            "encode_chromatin",
            "peak_tracks",
            [
                {
                    "label": "H3K4me3",
                    "accession": "ENCFF1",
                    "peak_type": "bed",
                }
            ],
            "peak_type",
        ),
    ],
)
def test_bad_options_rejected(tmp_path, source, field, value, match):
    def edit(raw):
        raw["projects"]["COAD"][source][field] = value

    with pytest.raises(ValueError, match=match):
        load_registry(_registry(tmp_path, edit))


def test_bad_default_window_rejected(tmp_path):
    def edit(raw):
        raw["defaults"]["roadmap"]["windows"] = ["utr"]

    with pytest.raises(ValueError, match="roadmap.windows"):
        load_registry(_registry(tmp_path, edit))


def test_null_transform_is_the_default(tmp_path):
    def edit(raw):
        raw["projects"]["COAD"]["roadmap"]["transform"] = None

    coad = load_registry(_registry(tmp_path, edit))["COAD"]
    assert coad.roadmap.transform is None
