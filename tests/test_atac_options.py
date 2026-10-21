"""ATAC sources take windows and transforms, as the other chromatin
sources do."""

import json

import pandas as pd
import pytest
from conftest import write_bigwig, write_gtf

from sigmutselcovs.builder import build_covariate_matrix
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import (
    load_registry,
    location_projects_registry,
)


def _registry(tmp_path, **atac):
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["atac"].update(atac)
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    return path


def test_registry_reads_atac_options(tmp_path):
    reg = _registry(
        tmp_path, windows=["cds_span"], transform="log1p_mean"
    )
    atac = load_registry(reg)["COAD"].atac
    assert atac.windows == ("cds_span",)
    assert atac.transform == "log1p_mean"


def test_registry_rejects_bad_atac_options(tmp_path):
    with pytest.raises(ValueError, match="atac.windows"):
        load_registry(_registry(tmp_path, windows=["exon"]))


def test_build_uses_atac_options(tmp_path):
    paths = project_paths(tmp_path / "data")
    paths.atac_dir.mkdir(parents=True)
    for seed, sample in enumerate(("COAD_ABC_T1", "COAD_DEF_T1")):
        write_bigwig(paths.atac_dir / f"{sample}.insertions.bw", seed)
    gtf = write_gtf(tmp_path / "g.gtf")
    reg = _registry(
        tmp_path, windows=["cds_span"], transform="log1p_mean"
    )
    matrices = build_covariate_matrix(
        "COAD",
        paths.root,
        registry_path=reg,
        gencode_gtfs={"hg19": gtf, "hg38": gtf},
        include=["atac"],
        validate_output=False,
    )
    assert list(matrices.full.columns) == [
        "coad_abc_t1_insertions_cds_span_logmean",
        "coad_def_t1_insertions_cds_span_logmean",
    ]
    assert list(matrices.tcga.columns) == ["atac_cds_span_logmean"]
    pd.testing.assert_series_equal(
        matrices.tcga["atac_cds_span_logmean"],
        matrices.full.mean(axis=1),
        check_names=False,
    )
