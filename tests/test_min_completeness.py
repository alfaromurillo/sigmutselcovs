"""``min_completeness`` drops chromatin columns missing too many genes."""

import json

import pytest
from conftest import (
    CHROM_SIZES,
    legacy_registry_raw,
    write_bigwig,
    write_gtf,
)

from sigmutselcovs.builder import build_covariate_matrix
from sigmutselcovs.paths import project_paths
from sigmutselcovs.registry import (
    load_registry,
)


def _chr1_only_bigwig(path):
    import pyBigWig

    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(list(CHROM_SIZES.items()))
    # varying along chr1, so the variance fix keeps the column
    starts = list(range(0, CHROM_SIZES["chr1"], 1000))
    ends = [s + 1000 for s in starts]
    values = [float(1 + k % 7) for k in range(len(starts))]
    bw.addEntries(
        ["chr1"] * len(starts), starts, ends=ends, values=values
    )
    bw.close()


def _registry(tmp_path, **roadmap):
    raw = legacy_registry_raw()
    raw["projects"]["COAD"]["roadmap"].update(roadmap)
    path = tmp_path / "projects.json"
    path.write_text(json.dumps(raw))
    return path


@pytest.mark.parametrize("bad", [0, 1.5, "0.9", True])
def test_registry_rejects_bad_fractions(tmp_path, bad):
    with pytest.raises(ValueError, match="min_completeness"):
        load_registry(_registry(tmp_path, min_completeness=bad))


def test_build_drops_incomplete_tracks(tmp_path):
    paths = project_paths(tmp_path / "data")
    paths.roadmap_dir.mkdir(parents=True)
    write_bigwig(paths.roadmap_dir / "E075-H3K4me3.fc.signal.bigwig")
    _chr1_only_bigwig(
        paths.roadmap_dir / "E106-H3K4me3.fc.signal.bigwig"
    )
    gtf = write_gtf(tmp_path / "g.gtf")

    def build(**roadmap):
        return build_covariate_matrix(
            "COAD",
            paths.root,
            registry_path=_registry(
                tmp_path,
                windows=["body"],
                transform="log1p_mean",
                **roadmap,
            ),
            gencode_gtfs={"hg19": gtf, "hg38": gtf},
            include=["roadmap"],
            validate_output=False,
        ).full

    both = [
        "e075_h3k4me3_fc_signal_body_logmean",
        "e106_h3k4me3_fc_signal_body_logmean",
    ]
    assert list(build().columns) == both
    # the chr1-only track has values for 3 of the 5 genes
    assert list(build(min_completeness=0.9).columns) == both[:1]
    assert list(build(min_completeness=0.5).columns) == both
