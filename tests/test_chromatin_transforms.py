"""Linear and log columns from one source (`transform` as a list)."""

import json

import pandas as pd
import pytest
from conftest import write_bigwig, write_gtf

from sigmutselcovs import covariates_chromatin as cc
from sigmutselcovs.builder import build_covariate_matrix
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import (
    chromatin_transforms,
    load_registry,
    location_projects_registry,
)


def _registry(tmp_path, **roadmap):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["roadmap"].update(roadmap)
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    return path


def test_chromatin_transforms():
    assert chromatin_transforms(None) == (None,)
    assert chromatin_transforms("log1p") == ("log1p",)
    assert chromatin_transforms([None, "log1p"]) == (None, "log1p")
    assert chromatin_transforms(("log1p", None)) == ("log1p", None)


def test_registry_takes_a_list(tmp_path):
    path = _registry(tmp_path, transform=[None, "log1p"])
    coad = load_registry(path)["COAD"]
    assert coad.roadmap.transform == (None, "log1p")
    scalar = load_registry(_registry(tmp_path, transform="log1p"))
    assert scalar["COAD"].roadmap.transform == "log1p"


@pytest.mark.parametrize(
    "value, match",
    [
        ([], "must not be empty"),
        ([None, "log2"], "not one of"),
        ([None, None], "repeats"),
        (["log1p", "log1p"], "repeats"),
    ],
)
def test_registry_rejects_bad_lists(tmp_path, value, match):
    with pytest.raises(ValueError, match=match):
        load_registry(_registry(tmp_path, transform=value))


def test_check_transforms():
    assert cc.check_transforms([None, "log1p"]) == (None, "log1p")
    for bad in ([], ["sqrt"], [None, None]):
        with pytest.raises(ValueError):
            cc.check_transforms(bad)


def test_cache_path_covers_the_list(tmp_path):
    base = tmp_path / "chromatin_covs_roadmap.csv"
    path = cc.chromatin_cache_path
    # one-element lists are their scalars
    assert path(base, transform=[None]) == base
    assert path(base, transform=["log1p"]) == path(
        base, transform="log1p"
    )
    both = path(base, transform=[None, "log1p"])
    assert (
        both.name
        == "chromatin_covs_roadmap_body-promoter_linear-log1p.csv"
    )
    assert path(base, transform=["log1p", None]) != both
    assert len({both, base, path(base, transform="log1p")}) == 3


def test_summaries_are_the_scalar_ones_side_by_side(
    synthetic_bigwig, synthetic_gtf
):
    path, _ = synthetic_bigwig
    genes = cc._load_gene_bodies(
        synthetic_gtf, ("protein_coding",), with_cds_span=True
    )
    windows = ["body", "cds_span"]
    both = cc.summarise_tracks_to_genes(
        [path], genes, windows=windows, transform=[None, "log1p"]
    )
    linear = cc.summarise_tracks_to_genes(
        [path], genes, windows=windows
    )
    log = cc.summarise_tracks_to_genes(
        [path], genes, windows=windows, transform="log1p"
    )
    pd.testing.assert_frame_equal(
        both, pd.concat([linear, log], axis=1), check_exact=True
    )
    streamed = cc.summarize_tracks_to_genes_streaming(
        [(path.stem, path)],
        genes,
        windows=windows,
        transform=[None, "log1p"],
    )
    pd.testing.assert_frame_equal(both, streamed, check_exact=True)


def test_build_emits_both(tmp_path):
    paths = project_paths(tmp_path / "data")
    paths.roadmap_dir.mkdir(parents=True)
    for seed, track in enumerate(("E075-H3K4me3", "E106-H3K4me3")):
        write_bigwig(
            paths.roadmap_dir / f"{track}.fc.signal.bigwig", seed
        )
    gtf = write_gtf(tmp_path / "g.gtf")
    gtfs = {"hg19": gtf, "hg38": gtf}

    def build(registry=None, **kwargs):
        return build_covariate_matrix(
            "COAD",
            paths.root,
            registry_path=registry,
            gencode_gtfs=gtfs,
            include=["roadmap"],
            apply_fixes=False,
            **kwargs,
        )

    reg = _registry(tmp_path, transform=[None, "log1p"])
    full = build(
        reg, cache_matrices=True, strict_validation=True
    ).full
    assert list(full.columns) == [
        f"{eid}_h3k4me3_fc_signal_{w}{log}"
        for eid in ("e075", "e106")
        for log in ("", "_log")
        for w in ("body", "promoter")
    ]
    caches = sorted(
        p.name for p in paths.roadmap_dir.glob("chromatin_covs_*.csv")
    )
    assert caches == [
        "chromatin_covs_roadmap_body-promoter_linear-log1p.csv"
    ]
    # the linear half is the default build's, the log half the
    # scalar log build's (each read back from its own cache file)
    default = build(validate_output=False).full
    log = build(
        _registry(tmp_path, transform="log1p"), validate_output=False
    ).full
    pd.testing.assert_frame_equal(
        full[default.columns], default, rtol=1e-15
    )
    pd.testing.assert_frame_equal(full[log.columns], log, rtol=1e-15)
    dictionary = pd.read_csv(paths.column_dictionary_csv)
    assert dictionary["description"].str.contains("log1p").all()
