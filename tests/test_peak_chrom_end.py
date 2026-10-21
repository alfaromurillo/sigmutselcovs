"""Peak windows clipped at the chromosome end (packaged chrom.sizes)."""

import json

import numpy as np
import pandas as pd
import pytest
from conftest import (
    CHROM_SIZES,
    GENE_IDS,
    legacy_registry_raw,
    write_gtf,
    write_peaks,
)

from sigmutselcovs import builder
from sigmutselcovs import covariates_chromatin as cc
from sigmutselcovs.builder import build_covariate_matrix
from sigmutselcovs.covariates_utilities import load_chrom_sizes
from sigmutselcovs.paths import project_paths

# a peak at chr2's very end (chr2 is 12,000 bases long)
END_PEAKS = [("chr1", 2_300, 2_800), ("chr2", 11_800, 12_000)]
G5 = GENE_IDS[
    4
]  # minus strand, chr2 [6000, 11500), CDS [9000, 11000)


@pytest.fixture
def genes(synthetic_gtf):
    return cc._load_gene_bodies(
        synthetic_gtf, ("protein_coding",), with_cds_span=True
    )


@pytest.mark.parametrize("assembly", ["hg19", "hg38"])
def test_packaged_chrom_sizes(assembly):
    sizes = load_chrom_sizes(assembly)
    primary = [f"chr{i}" for i in range(1, 23)] + [
        "chrX",
        "chrY",
        "chrM",
    ]
    assert list(sizes) == primary
    expected = {
        "hg19": (249_250_621, 155_270_560, 16_571),
        "hg38": (248_956_422, 156_040_895, 16_569),
    }[assembly]
    assert (sizes["chr1"], sizes["chrX"], sizes["chrM"]) == expected


def test_chrom_sizes_from_a_file(tmp_path):
    path = tmp_path / "toy.chrom.sizes"
    path.write_text("chrA\t100\nchrB\t50\n")
    assert load_chrom_sizes(path) == {"chrA": 100, "chrB": 50}
    with pytest.raises(ValueError, match="No chromosome sizes"):
        load_chrom_sizes("mm10")


def test_window_past_the_end_is_clipped(tmp_path, genes):
    peaks = write_peaks(tmp_path / "t.narrowPeak", peaks=END_PEAKS)
    windows = ["body", "promoter", "cds_span", "cds_span2k"]
    clipped = cc.summarize_peak_coverage(
        peaks, genes, windows=windows, chrom_sizes=CHROM_SIZES
    )
    unclipped = cc.summarize_peak_coverage(
        peaks, genes, windows=windows
    )
    # G5's promoter [11300, 13500) and CDS span +-2 kb [7000, 13000)
    # run past chr2's end: 200 peak bases over the 700 and 5000 bases
    # that exist, not over 2200 and 6000
    assert clipped.loc[G5, "t_narrowpeak_promoter"] == 200 / 700
    assert unclipped.loc[G5, "t_narrowpeak_promoter"] == 200 / 2200
    assert clipped.loc[G5, "t_narrowpeak_cds_span2k"] == 200 / 5000
    assert unclipped.loc[G5, "t_narrowpeak_cds_span2k"] == 200 / 6000
    # windows inside their chromosome are untouched
    inside = clipped.drop(index=G5)
    pd.testing.assert_frame_equal(
        inside, unclipped.drop(index=G5), check_exact=True
    )
    assert clipped.loc[G5, "t_narrowpeak_body"] == 0.0


def test_window_starting_past_the_end_is_nan(tmp_path, genes):
    peaks = write_peaks(tmp_path / "t.narrowPeak", peaks=END_PEAKS)
    short = {"chr1": 20_000, "chr2": 5_000}  # G5 lies beyond chr2
    out = cc.summarize_peak_coverage(
        peaks, genes, windows=["body"], chrom_sizes=short
    )
    assert np.isnan(out.loc[G5, "t_narrowpeak_body"])


def test_unsized_chromosome_is_not_clipped(tmp_path, genes, caplog):
    peaks = write_peaks(tmp_path / "t.narrowPeak", peaks=END_PEAKS)
    out = cc.summarize_peak_coverage(
        peaks,
        genes,
        windows=["promoter"],
        chrom_sizes={"chr1": 20_000},
    )
    assert out.loc[G5, "t_narrowpeak_promoter"] == 200 / 2200
    assert any(
        "No chromosome size for chr2" in m for m in caplog.messages
    )


def test_builder_clips_with_the_assembly_table(tmp_path, monkeypatch):
    paths = project_paths(tmp_path / "data")
    paths.roadmap_peaks_dir.mkdir(parents=True)
    write_peaks(
        paths.roadmap_peaks_dir / "E075-H3K4me3.narrowPeak.gz",
        peaks=END_PEAKS,
        gzipped=True,
    )
    gtf = write_gtf(tmp_path / "g.gtf")
    asked = []

    def sizes(assembly):
        asked.append(assembly)
        return dict(CHROM_SIZES)

    # the synthetic genome's chromosomes are far shorter than hg19's
    monkeypatch.setattr(builder, "load_chrom_sizes", sizes)
    raw = legacy_registry_raw()
    raw["projects"]["COAD"]["roadmap"]["peaks"] = ["narrowPeak"]
    registry = tmp_path / "projects.json"
    registry.write_text(json.dumps(raw))
    full = build_covariate_matrix(
        "COAD",
        paths.root,
        registry_path=registry,
        gencode_gtfs={"hg19": gtf, "hg38": gtf},
        include=["roadmap"],
        apply_fixes=False,
        validate_output=False,
    ).full
    assert asked == ["hg19"]
    assert (
        full.loc[G5, "e075_h3k4me3_narrowpeak_promoter"] == 200 / 700
    )
