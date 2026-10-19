"""Chromatin summary windows, the log-mean and peak coverage.

Every summary is checked against a per-base numpy computation on the
synthetic genome of ``conftest.py`` (and the linear means against
``pyBigWig.stats(..., exact=True)``), and the default options against
the original body/promoter path, value for value.
"""

import numpy as np
import pandas as pd
import pyBigWig
import pytest
from conftest import CHROM_SIZES, GENE_IDS, write_gtf, write_peaks

from sigmutselcovs import covariates_chromatin as cc
from sigmutselcovs.covariates_utilities import (
    load_gene_bodies_from_gtf,
)

WINDOWS = ("body", "promoter", "cds_span", "cds_span2k")


@pytest.fixture
def genes(synthetic_gtf):
    return cc._load_gene_bodies(
        synthetic_gtf, ("protein_coding",), with_cds_span=True
    )


def _reference(dense, chrom, start, end, transform):
    """Per-base mean over the covered bases of [start, end), clipped
    to the chromosome."""
    if not (np.isfinite(start) and np.isfinite(end)):
        return np.nan
    start = int(max(start, 0))
    end = int(min(end, CHROM_SIZES[chrom]))
    values = dense[chrom][start:end]
    values = values[~np.isnan(values)]
    if not values.size:
        return np.nan
    if transform == "log1p":
        values = np.log1p(values)
    return float(values.mean())


def _assert_close(got, expected):
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    assert np.array_equal(np.isnan(got), np.isnan(expected))
    ok = ~np.isnan(expected)
    rel = np.abs(got[ok] - expected[ok]) / np.maximum(
        np.abs(expected[ok]), 1e-12
    )
    assert rel.max(initial=0.0) < 1e-6


# --- gene table ---


def test_cds_span_columns(synthetic_gtf):
    genes = load_gene_bodies_from_gtf(
        synthetic_gtf,
        biotypes=["protein_coding"],
        autosomes_only=False,
        with_cds_span=True,
    )
    assert list(genes.index) == GENE_IDS
    spans = genes[["cds_start", "cds_end"]]
    assert str(spans["cds_start"].dtype) == "Int64"
    # merged across transcripts, 0-based half-open
    assert tuple(spans.loc["ENSG00000000001"]) == (2400, 5600)
    assert tuple(spans.loc["ENSG00000000002"]) == (10000, 12500)
    assert spans.loc["ENSG00000000003"].isna().all()
    assert tuple(spans.loc["ENSG00000000004"]) == (800, 2000)
    assert tuple(spans.loc["ENSG00000000005"]) == (9000, 11000)


def test_cds_span_leaves_default_table_unchanged(synthetic_gtf):
    kwargs = {"biotypes": ["protein_coding"], "autosomes_only": False}
    plain = load_gene_bodies_from_gtf(synthetic_gtf, **kwargs)
    with_cds = load_gene_bodies_from_gtf(
        synthetic_gtf, with_cds_span=True, **kwargs
    )
    assert list(plain.columns) == [
        "Chromosome",
        "start",
        "end",
        "strand",
    ]
    pd.testing.assert_frame_equal(plain, with_cds[plain.columns])


def test_cds_span_keeps_par_copies_apart(tmp_path):
    gtf = write_gtf(
        tmp_path / "par.gtf",
        [
            (
                "ENSG00000000009.2",
                "chrX",
                1001,
                5000,
                "+",
                "protein_coding",
                {"TX": [(1201, 1300)]},
            ),
            (
                "ENSG00000000009.2",
                "chrY",
                2001,
                6000,
                "+",
                "protein_coding",
                {"TY": [(4001, 4500)]},
            ),
        ],
    )
    genes = load_gene_bodies_from_gtf(
        gtf, autosomes_only=False, with_cds_span=True
    )
    by_chrom = genes.set_index("Chromosome")
    assert tuple(by_chrom.loc["chrX", ["cds_start", "cds_end"]]) == (
        1200,
        1300,
    )
    assert tuple(by_chrom.loc["chrY", ["cds_start", "cds_end"]]) == (
        4000,
        4500,
    )


# --- windows ---


def test_gene_windows(genes):
    def coords(window, gene):
        row = cc.gene_windows(genes, window).loc[gene]
        return (row["start"], row["end"])

    g1, g2, g3, g4, g5 = GENE_IDS
    assert coords("body", g1) == (2000, 7000)
    assert coords("promoter", g1) == (0, 2200)
    assert coords("promoter", g2) == (14800, 17000)  # minus strand
    assert coords("cds_span", g1) == (2400, 5600)
    assert coords("cds_span2k", g1) == (400, 7600)
    assert coords("cds_span2k", g4) == (0, 4000)  # clipped at 0
    assert coords("cds_span2k", g5) == (7000, 13000)  # past the end
    assert np.isnan(coords("cds_span", g3)).all()
    assert np.isnan(coords("cds_span2k", g3)).all()


def test_cds_windows_need_cds_columns(genes):
    plain = genes.drop(columns=["cds_start", "cds_end"])
    with pytest.raises(ValueError, match="with_cds_span"):
        cc.gene_windows(plain, "cds_span")


@pytest.mark.parametrize(
    "windows, error",
    [
        (["body", "exon"], ValueError),
        ([], ValueError),
        (["body", "body"], ValueError),
        ("body", TypeError),
    ],
)
def test_check_windows_rejects(windows, error):
    with pytest.raises(error):
        cc.check_windows(windows)


def test_check_transform_rejects():
    assert cc.check_transform(None) is None
    assert cc.check_transform("log1p") == "log1p"
    with pytest.raises(ValueError, match="transform"):
        cc.check_transform("log2")


# --- bigWig summaries ---


@pytest.mark.parametrize("transform", [None, "log1p"])
def test_window_summaries_match_per_base_numpy(
    synthetic_bigwig, genes, transform
):
    path, dense = synthetic_bigwig
    out = cc.summarize_bigwig_windows(
        path, genes, windows=WINDOWS, transform=transform, label="trk"
    )
    suffix = "_log" if transform else ""
    assert list(out.columns) == [f"trk_{w}{suffix}" for w in WINDOWS]
    with pyBigWig.open(str(path)) as bw:
        for window in WINDOWS:
            coords = cc.gene_windows(genes, window)
            expected = []
            for gene, row in coords.iterrows():
                legacy = transform is None and window in (
                    "body",
                    "promoter",
                )
                if (
                    legacy
                    and row["end"] > CHROM_SIZES[row["Chromosome"]]
                ):
                    # the original path: NaN past the chromosome end
                    expected.append(np.nan)
                    continue
                expected.append(
                    _reference(
                        dense,
                        row["Chromosome"],
                        row["start"],
                        row["end"],
                        transform,
                    )
                )
                if transform is None and np.isfinite(row["start"]):
                    end = min(
                        row["end"], CHROM_SIZES[row["Chromosome"]]
                    )
                    stat = bw.stats(
                        row["Chromosome"],
                        int(row["start"]),
                        int(end),
                        type="mean",
                        exact=True,
                    )[0]
                    got = out.loc[gene, f"trk_{window}"]
                    if stat is None:  # no covered base
                        assert np.isnan(got)
                    else:
                        assert got == stat
            _assert_close(out[f"trk_{window}{suffix}"], expected)


def test_log_mean_is_per_base_not_log_of_mean(
    synthetic_bigwig, genes
):
    path, _ = synthetic_bigwig
    linear = cc.summarize_bigwig_windows(
        path, genes, windows=["body"]
    )
    log = cc.summarize_bigwig_windows(
        path, genes, windows=["body"], transform="log1p"
    )
    stem = cc.sanitize_feature_label(path.stem)
    log_of_mean = np.log1p(linear[f"{stem}_body"])
    # Jensen: the mean of the logs is below the log of the mean
    assert (log[f"{stem}_body_log"] < log_of_mean).all()


def test_log_mean_chunking(synthetic_bigwig, genes, monkeypatch):
    path, _ = synthetic_bigwig
    whole = cc.summarize_bigwig_windows(
        path, genes, windows=WINDOWS, transform="log1p"
    )
    # chunks smaller than most windows: each long window is read alone
    monkeypatch.setattr(cc, "_LOG1P_CHUNK_BASES", 1_000)
    chunked = cc.summarize_bigwig_windows(
        path, genes, windows=WINDOWS, transform="log1p"
    )
    pd.testing.assert_frame_equal(whole, chunked, rtol=1e-12)


def test_absent_chromosome_is_nan(synthetic_bigwig, genes, caplog):
    path, _ = synthetic_bigwig
    moved = genes.copy()
    moved.loc["ENSG00000000001", "Chromosome"] = "chr7"
    out = cc.summarize_bigwig_windows(
        path,
        moved,
        windows=["cds_span"],
        transform="log1p",
        label="t",
    )
    assert np.isnan(out.loc["ENSG00000000001", "t_cds_span_log"])
    assert out.drop(index="ENSG00000000001").notna().sum().item() == 3
    assert any("chr7" in m for m in caplog.messages)


def test_default_windows_equal_original_path(synthetic_bigwig, genes):
    path, _ = synthetic_bigwig
    original = cc.summarise_tracks_to_genes([path], genes)
    explicit = cc.summarise_tracks_to_genes(
        [path], genes, windows=["body", "promoter"]
    )
    pd.testing.assert_frame_equal(
        original, explicit, check_exact=True
    )
    assert list(original.columns) == [
        "e900_h3k4me3_fc_signal_body",
        "e900_h3k4me3_fc_signal_promoter",
    ]


def test_streaming_windows_match_batch(synthetic_bigwig, genes):
    path, _ = synthetic_bigwig
    kwargs = {"windows": ["cds_span", "body"], "transform": "log1p"}
    batch = cc.summarise_tracks_to_genes([path], genes, **kwargs)
    streamed = cc.summarize_tracks_to_genes_streaming(
        [(path.stem, path)], genes, **kwargs
    )
    pd.testing.assert_frame_equal(batch, streamed, check_exact=True)


def test_load_or_generate_with_windows(
    synthetic_bigwig, synthetic_gtf, tmp_path
):
    path, dense = synthetic_bigwig
    cache = tmp_path / "covs.csv"
    out = cc.load_or_generate_chromatin_covariates(
        cache,
        [path],
        synthetic_gtf,
        windows=["cds_span2k"],
        transform="log1p",
    )
    assert list(out.columns) == [
        "e900_h3k4me3_fc_signal_cds_span2k_log"
    ]
    assert list(out.index) == GENE_IDS
    expected = _reference(dense, "chr2", 0, 4000, "log1p")
    _assert_close([out.loc["ENSG00000000004"].item()], [expected])
    reread = cc.load_or_generate_chromatin_covariates(
        cache, [], synthetic_gtf
    )
    pd.testing.assert_frame_equal(out, reread, rtol=1e-12)


# --- peaks ---

# Hand count over the merged peaks of conftest.PEAKS
# (chr1: [2300, 2800) [3000, 3100) [14900, 15100);
#  chr2: [0, 100) [900, 1200)).
EXPECTED_COVERAGE = {
    "body": [600 / 5000, 100 / 5500, 0.0, 400 / 3000, 0.0],
    "promoter": [0.0, 200 / 2200, 200 / 2200, 100 / 200, 0.0],
    "cds_span": [500 / 3200, 0.0, np.nan, 300 / 1200, 0.0],
    "cds_span2k": [600 / 7200, 0.0, np.nan, 400 / 4000, 0.0],
}


@pytest.mark.parametrize("gzipped", [False, True])
def test_peak_coverage_hand_count(tmp_path, genes, gzipped):
    name = "E900-H3K4me3.narrowPeak" + (".gz" if gzipped else "")
    peaks = write_peaks(tmp_path / name, gzipped=gzipped)
    out = cc.summarize_peak_coverage(peaks, genes, windows=WINDOWS)
    assert list(out.columns) == [
        f"e900_h3k4me3_narrowpeak_{w}" for w in WINDOWS
    ]
    for window, expected in EXPECTED_COVERAGE.items():
        got = out[f"e900_h3k4me3_narrowpeak_{window}"].to_numpy()
        np.testing.assert_allclose(got, expected, rtol=0, atol=1e-15)


def test_read_peak_intervals_merges(tmp_path):
    peaks = write_peaks(tmp_path / "x.broadPeak")
    merged = cc.read_peak_intervals(peaks)
    assert merged["chr1"][0].tolist() == [2300, 3000, 14900]
    assert merged["chr1"][1].tolist() == [2800, 3100, 15100]
    assert merged["chr2"][0].tolist() == [0, 900]
    assert merged["chr2"][1].tolist() == [100, 1200]


def test_peak_file_parts():
    assert cc.peak_file_parts("E075-H3K4me3.narrowPeak.gz") == (
        "E075-H3K4me3",
        "narrowPeak",
    )
    assert cc.peak_file_parts("/d/H3K27ac_ENCFF1.broadPeak") == (
        "H3K27ac_ENCFF1",
        "broadPeak",
    )
    with pytest.raises(ValueError, match="not a peak file"):
        cc.peak_file_parts("E075-H3K4me3.bed.gz")


def test_peak_chromosome_without_peaks_is_zero(
    tmp_path, genes, caplog
):
    peaks = write_peaks(
        tmp_path / "a.narrowPeak",
        peaks=[("chr1", 2_300, 2_800)],
    )
    out = cc.summarize_peak_coverage(peaks, genes, windows=["body"])
    assert out["a_narrowpeak_body"].tolist() == [
        0.1,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    assert any("chr2" in m for m in caplog.messages)


def test_load_or_generate_peak_covariates(tmp_path, synthetic_gtf):
    narrow = write_peaks(
        tmp_path / "E900-H3K4me3.narrowPeak.gz", gzipped=True
    )
    broad = write_peaks(tmp_path / "E900-H3K4me3.broadPeak")
    cache = tmp_path / "peaks.csv"
    out = cc.load_or_generate_peak_covariates(
        cache, [narrow, broad], synthetic_gtf, windows=["cds_span"]
    )
    assert list(out.columns) == [
        "e900_h3k4me3_narrowpeak_cds_span",
        "e900_h3k4me3_broadpeak_cds_span",
    ]
    np.testing.assert_allclose(
        out.iloc[:, 0], EXPECTED_COVERAGE["cds_span"], atol=1e-15
    )
    assert cache.exists()


# --- average_by_assay ---


def test_collapse_by_assay_keeps_windows_transforms_and_peaks_apart():
    idx = pd.Index(["ENSG_A"], name="ensembl_gene_id")
    cols = {
        "e075_h3k4me3_fc_signal_body": 1.0,
        "e106_h3k4me3_fc_signal_body": 3.0,
        "e075_h3k4me3_fc_signal_cds_span": 5.0,
        "e106_h3k4me3_fc_signal_cds_span2k": 7.0,
        "e075_h3k4me3_fc_signal_cds_span_log": 0.5,
        "h3k4me3_fc_signal_encff1_cds_span_log": 1.5,
        "e075_h3k4me3_narrowpeak_body": 0.2,
        "h3k4me3_encff2_narrowpeak_body": 0.4,
        "e075_h3k4me3_broadpeak_promoter": 0.6,
        "dnase_fc_signal_encff3_promoter_log": 2.0,
    }
    frame = pd.DataFrame({k: [v] for k, v in cols.items()}, index=idx)
    collapsed = cc._collapse_by_assay(frame)
    assert collapsed.iloc[0].to_dict() == {
        "h3k4me3_body": 2.0,
        "h3k4me3_cds_span": 5.0,
        "h3k4me3_cds_span2k": 7.0,
        "h3k4me3_cds_span_log": 1.0,
        "h3k4me3_narrowpeak_body": pytest.approx(0.3),
        "h3k4me3_broadpeak_promoter": 0.6,
        "dnase_promoter_log": 2.0,
    }


# --- cache names ---


def test_chromatin_cache_path(tmp_path):
    base = tmp_path / "chromatin_covs_roadmap.csv"
    assert cc.chromatin_cache_path(base) == base
    assert (
        cc.chromatin_cache_path(base, windows=["body", "promoter"])
        == base
    )
    assert cc.chromatin_cache_path(
        base, windows=["cds_span", "cds_span2k"], transform="log1p"
    ) == tmp_path / (
        "chromatin_covs_roadmap_cds_span-cds_span2k_log1p.csv"
    )
    assert cc.chromatin_cache_path(
        base, transform="log1p"
    ) == tmp_path / ("chromatin_covs_roadmap_body-promoter_log1p.csv")
    assert cc.chromatin_cache_path(
        base, peak_types=["narrowPeak"]
    ) == tmp_path / (
        "chromatin_covs_roadmap_peaks_narrowpeak_body-promoter.csv"
    )
    # window order is column order: a different cache
    assert cc.chromatin_cache_path(
        base, windows=["promoter", "body"]
    ) != cc.chromatin_cache_path(base)
