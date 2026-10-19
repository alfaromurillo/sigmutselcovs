"""Synthetic genome shared by the chromatin tests (no network).

Two chromosomes with stretches of no data, a bigWig with values from
a fixed RNG, a small GENCODE-style GTF and a narrowPeak file with
overlapping peaks. `dense_signal` returns the same signal per base
(NaN where the track has no data) for numpy reference computations.
"""

import gzip

import numpy as np
import pytest

CHROM_SIZES = {"chr1": 20_000, "chr2": 12_000}
# [start, end) stretches the bigWig leaves without data.
GAPS = {
    "chr1": [(0, 300), (8_000, 9_000), (15_500, 15_520)],
    "chr2": [(0, 500), (11_000, 11_400)],
}

# (gene_id, chrom, start, end, strand, gene_type, CDS segments by
# transcript), GTF 1-based inclusive.
GTF_GENES = [
    (
        "ENSG00000000001.5",
        "chr1",
        2001,
        7000,
        "+",
        "protein_coding",
        {
            "T1": [(2501, 2600), (3001, 3200), (4001, 4100)],
            "T2": [(2401, 2450), (5001, 5600)],
        },
    ),
    (
        "ENSG00000000002.3",
        "chr1",
        9501,
        15000,
        "-",
        "protein_coding",
        {"T3": [(10001, 10200), (12001, 12500)]},
    ),
    # protein-coding gene without CDS features: NaN CDS windows
    (
        "ENSG00000000003.1",
        "chr1",
        16001,
        18000,
        "+",
        "protein_coding",
        {},
    ),
    # CDS span +-2 kb runs off the chromosome start
    (
        "ENSG00000000004.2",
        "chr2",
        1,
        3000,
        "+",
        "protein_coding",
        {"T4": [(801, 1000), (1501, 2000)], "T5": [(1801, 1900)]},
    ),
    # minus strand near the chromosome end: the promoter and the
    # CDS span +-2 kb run past it
    (
        "ENSG00000000005.7",
        "chr2",
        6001,
        11500,
        "-",
        "protein_coding",
        {"T6": [(9001, 11000)]},
    ),
    # filtered out by biotypes=("protein_coding",)
    (
        "ENSG00000000006.1",
        "chr1",
        500,
        1500,
        "+",
        "lncRNA",
        {},
    ),
]

GENE_IDS = [
    g[0].split(".")[0] for g in GTF_GENES if g[5] == "protein_coding"
]


def _attrs(gene_id, gene_type, transcript=None):
    out = f'gene_id "{gene_id}"; gene_type "{gene_type}";'
    if transcript is not None:
        out += f' transcript_id "{transcript}";'
    return out


def write_gtf(path, genes=GTF_GENES):
    lines = ["##description: synthetic test annotation"]
    for gene_id, chrom, start, end, strand, gtype, cds in genes:
        lines.append(
            "\t".join(
                [
                    chrom,
                    "TEST",
                    "gene",
                    str(start),
                    str(end),
                    ".",
                    strand,
                    ".",
                    _attrs(gene_id, gtype),
                ]
            )
        )
        for transcript, segments in cds.items():
            for s, e in segments:
                lines.append(
                    "\t".join(
                        [
                            chrom,
                            "TEST",
                            "CDS",
                            str(s),
                            str(e),
                            ".",
                            strand,
                            "0",
                            _attrs(gene_id, gtype, transcript),
                        ]
                    )
                )
    path.write_text("\n".join(lines) + "\n")
    return path


def _entries(chrom, rng):
    """bedGraph-style runs of random length and value, skipping the
    chromosome's gaps (and a few short random holes)."""
    size = CHROM_SIZES[chrom]
    gaps = GAPS[chrom]
    starts, ends, values = [], [], []
    pos = 0
    while pos < size:
        gap = next((g for g in gaps if g[0] <= pos < g[1]), None)
        if gap is not None:
            pos = gap[1]
            continue
        run = int(rng.integers(1, 60))
        end = min(pos + run, size)
        nxt = next((g[0] for g in gaps if pos < g[0] < end), None)
        if nxt is not None:
            end = nxt
        if rng.random() < 0.05:  # a short hole with no data
            pos = end
            continue
        value = (
            0.0 if rng.random() < 0.1 else float(rng.gamma(1.5, 2.0))
        )
        starts.append(pos)
        ends.append(end)
        values.append(value)
        pos = end
    return starts, ends, values


@pytest.fixture(scope="session")
def synthetic_bigwig(tmp_path_factory):
    """Path to the synthetic bigWig and its dense per-base signal."""
    import pyBigWig

    rng = np.random.default_rng(20261019)
    path = (
        tmp_path_factory.mktemp("bigwig")
        / "E900-H3K4me3.fc.signal.bigwig"
    )
    dense = {}
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(list(CHROM_SIZES.items()))
    for chrom, size in CHROM_SIZES.items():
        starts, ends, values = _entries(chrom, rng)
        bw.addEntries(
            [chrom] * len(starts), starts, ends=ends, values=values
        )
        signal = np.full(size, np.nan)
        for s, e, v in zip(starts, ends, values):
            # bigWig stores float32
            signal[s:e] = np.float32(v)
        dense[chrom] = signal
    bw.close()
    return path, dense


@pytest.fixture
def synthetic_gtf(tmp_path):
    return write_gtf(tmp_path / "synthetic.gtf")


# Overlapping, abutting and contained peaks; GENE windows' coverage is
# counted by hand in test_chromatin_windows.py.
PEAKS = [
    ("chr1", 2_300, 2_500),
    ("chr1", 2_450, 2_700),  # overlaps the previous one
    ("chr1", 2_700, 2_800),  # abuts it
    ("chr1", 3_000, 3_100),
    ("chr1", 3_050, 3_080),  # contained in the previous one
    ("chr1", 14_900, 15_100),
    ("chr2", 0, 100),
    ("chr2", 900, 1_200),
]


def write_peaks(path, peaks=PEAKS, gzipped=False):
    """narrowPeak rows (10 columns), out of order, after a track
    line."""
    rows = ['track type=narrowPeak name="synthetic"']
    for i, (chrom, s, e) in reversed(list(enumerate(peaks))):
        rows.append(
            f"{chrom}\t{s}\t{e}\tpeak{i}\t100\t.\t5.0\t3.0\t2.0\t50"
        )
    text = "\n".join(rows) + "\n"
    if gzipped:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path
