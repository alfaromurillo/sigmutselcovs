"""combine_with_generic(drop_own_tracks=True) matches tracks, not names."""

import pandas as pd

from sigmutselcovs.builder import combine_with_generic
from sigmutselcovs.covariates_chromatin import (
    chromatin_track_identity,
)


def _frame(sources):
    df = pd.DataFrame(
        {c: [1.0, 2.0] for c in sources}, index=["ENSG1", "ENSG2"]
    )
    df.attrs["column_sources"] = dict(sources)
    return df


OWN = _frame(
    {
        "e075_h3k4me3_fc_signal_body": "roadmap",
        "h3k27ac_fc_signal_encff001aaa_promoter": "encode_chromatin",
        "gtex_colon": "gtex",
    }
)
POOL = _frame(
    {
        "e075_h3k4me3_fc_signal_cds_span_logmean": "roadmap",
        "e075_h3k4me3_broadpeak_cds_span": "roadmap_peaks",
        "e106_h3k4me3_fc_signal_cds_span_logmean": "roadmap",
        "h3k27ac_fc_signal_encff001aaa_cds_span_logmean": "encode_chromatin",
        "dnase_fc_signal_encff002bbb_cds_span_logmean": "encode_chromatin",
        "rt_pool_z": "rt_pool",
    }
)


def test_identity_ignores_window_transform_and_peaks():
    for column in (
        "e075_h3k4me3_fc_signal_body",
        "e075_h3k4me3_fc_signal_cds_span_logmean",
        "e075_h3k4me3_fc_signal_cds_span2k_log",
        "e075_h3k4me3_broadpeak_cds_span",
    ):
        assert chromatin_track_identity(column) == "e075_h3k4me3"
    assert (
        chromatin_track_identity("h3k27ac_fc_signal_encff001aaa_body")
        == "h3k27ac_encff001aaa"
    )


def test_drops_the_cohorts_own_tracks():
    combined = combine_with_generic(OWN, POOL, drop_own_tracks=True)
    assert list(combined.columns) == [
        *OWN.columns,
        "e106_h3k4me3_fc_signal_cds_span_logmean",
        "dnase_fc_signal_encff002bbb_cds_span_logmean",
        "rt_pool_z",
    ]
    sources = combined.attrs["column_sources"]
    assert set(sources) == set(combined.columns)


def test_default_keeps_them():
    combined = combine_with_generic(OWN, POOL)
    assert combined.shape[1] == OWN.shape[1] + POOL.shape[1]
