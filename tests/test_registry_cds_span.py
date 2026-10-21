"""The shipped registry: one CDS-span method for every chromatin row."""

import pytest

from sigmutselcovs.registry import load_registry

pytestmark = pytest.mark.shipped_registry


def test_every_chromatin_row_is_on_the_cds_span():
    for code, spec in load_registry().items():
        for source in (
            spec.roadmap,
            spec.encode_chromatin,
            spec.atac,
        ):
            if source is None:
                continue
            assert source.windows == ("cds_span",), code
            assert source.transform == "log1p_mean", code
            assert source.min_completeness == 0.9, code
        assert not spec.average_by_assay, code


def test_roadmap_rows_take_every_consolidated_mark():
    registry = load_registry()
    marks = registry["GENERIC"].roadmap.marks
    assert len(marks) == 32
    assert {"DNase", "H3K27ac", "H3K79me2", "H2A.Z"} <= set(marks)
    assert registry["COAD"].roadmap.marks == marks
