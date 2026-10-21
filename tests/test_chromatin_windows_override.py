"""build_covariate_matrix(chromatin_windows=...) overrides every
chromatin source's windows."""

import json

import pytest
from conftest import write_bigwig, write_gtf

from sigmutselcovs.builder import build_covariate_matrix
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import location_projects_registry


def test_override_adds_the_promoter(tmp_path):
    paths = project_paths(tmp_path / "data")
    paths.roadmap_dir.mkdir(parents=True)
    write_bigwig(paths.roadmap_dir / "E075-H3K4me3.fc.signal.bigwig")
    gtf = write_gtf(tmp_path / "g.gtf")
    raw = json.loads(location_projects_registry.read_text())
    raw["projects"]["COAD"]["roadmap"].update(
        windows=["cds_span"], transform="log1p_mean"
    )
    reg = tmp_path / "projects.json"
    reg.write_text(json.dumps(raw))

    def build(**kwargs):
        return build_covariate_matrix(
            "COAD",
            paths.root,
            registry_path=reg,
            gencode_gtfs={"hg19": gtf, "hg38": gtf},
            include=["roadmap"],
            apply_fixes=False,
            validate_output=False,
            **kwargs,
        ).full

    stem = "e075_h3k4me3_fc_signal"
    assert list(build().columns) == [f"{stem}_cds_span_logmean"]
    assert list(
        build(chromatin_windows=["cds_span", "promoter"]).columns
    ) == [f"{stem}_cds_span_logmean", f"{stem}_promoter_logmean"]
    with pytest.raises(ValueError, match="Unknown chromatin windows"):
        build(chromatin_windows=["exon"])
