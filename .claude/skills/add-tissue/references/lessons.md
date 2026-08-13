# Lessons from the BRCA bring-up

Detailed reference for steps 4 and 5 of `SKILL.md`, plus gotchas
that apply throughout. Read the relevant section when you reach
that step -- no need to read this file front to back.

## Roadmap Epigenomics

URL pattern (confirmed working, use exactly this):

```
https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/
macs2signal/foldChange/{EID}-{MARK}.fc.signal.bigwig
```

The standard 7 marks COAD and BRCA both use:
`H3K4me1, H3K4me3, H3K9ac, H3K9me3, H3K27ac, H3K27me3, H3K36me3`.

**404s on specific `(EID, mark)` pairs are the norm, not a bug.**
COAD's E106 is missing H3K9ac; BRCA's E027 is missing H3K27ac and
E028 is missing both H3K9ac and H3K27ac. The download layer
(`download_roadmap_tracks` in `download.py`) already tolerates this
-- it logs a warning and continues unless the mark is listed in
`required_marks` (which is empty by default; leave it empty unless
you have a specific reason a mark absolutely must exist). Don't
treat a 404 here as a sign something is configured wrong -- verify
each track with a live request (`curl -I` or equivalent) before
writing it into the registry, but expect some fraction to be
missing regardless of how careful the tissue-to-epigenome matching
was.

To find candidate EIDs for a primary site: the Roadmap web portal
(https://egg2.wustl.edu/roadmap/web_portal/) has the full
127-epigenome list with tissue/cell-type labels. There is no
GDC-style API for this -- it's a lookup against that static list by
eyeballing the tissue description.

## Replication timing: two incompatible ENCODE cohorts

This is the part that's easy to get wrong if you don't already know
it. Querying `ReplicationTimingSeries` on the ENCODE portal
(`https://www.encodeproject.org/search/?type=ReplicationTimingSeries&format=json`)
returns series from **two different labs with incompatible output
formats** -- the series accession alone doesn't tell you which kind
you got; you have to check the files.

### Cohort A: UW/Stamatoyannopoulos, 6-fraction Repli-seq (hg19)

One percentage-normalized bigWig per S-phase fraction: **G1b, S1,
S2, S3, S4, G2**. This is the rich format -- it gives both a
per-gene MRT scalar and 6 CLR-transformed fraction columns. Cell
lines seen in this cohort as of the BRCA research (2026-08):
MCF-7 (used for BRCA), HeLa-S3 (CESC), K562 (LAML), GM12878 + 3
other LCLs (DLBC), IMR-90, keratinocyte/NHEK, SK-N-SH, BJ, HUVEC,
BG02, HepG2 (LIHC).

If the tissue you're adding has a cell line in this list, use it --
this is the `fraction_bigwigs` registry type, same as BRCA's row.

### Cohort B: Gilbert lab, 2-fraction early/late (hg19)

Despite being called "2-fraction," **the only processed output is a
single wavelet-smoothed log2(early/late) bigWig** -- the separate
early and late fraction files exist only as raw fastq/bam, not
usable signal tracks. This is the `wavelet` registry type: one
per-gene column, no MRT, no CLR fractions. Cell lines seen: LNCaP
(PRAD), Caki2 (KIRC), A549 + NCI-H460 (LUAD/LUSC), G401, SK-N-MC,
**T47D**.

**T47D is a dead end** -- it's listed as a breast cell line and was
the first thing checked for BRCA, but it has zero usable signal
files (only fastq/bam for its two "fraction" experiments, and the
series-level file is empty of processed output). MCF-7 (Cohort A)
was used instead. Don't spend time re-verifying T47D if it comes up
for a future tissue -- it's confirmed unusable.

### Also check GEO

`GSE137764` has 16-fraction Gaussian-smoothed MAT files (the
richest format, used for COAD/HCT116), but as of the BRCA work it
only covers **H1, H9, HCT116** (plus mouse). Worth a quick check for
any new tissue, but don't expect a hit.

### If nothing fits

Set `"repliseq": null`. A matrix missing the replication-timing
block is meaningfully better than one built from a mismatched cell
line (wrong tissue's chromatin state), and the builder already
handles the null case cleanly.

## Other gotchas (apply across all steps)

- **GDC's `/data/<uuid>` endpoint rejects HEAD requests with a 400,
  and `GET /files/<uuid>` only works for UUIDs that are actually
  indexed in search.** The ATAC tarball UUIDs are *not* indexed --
  `GET /files/<uuid>` 404s for them even though the blob is real
  (verified 2026-08 against the confirmed-working COAD/BRCA UUIDs).
  The reliable way to get metadata for any `/data/<uuid>` blob,
  indexed or not, is a 1-byte ranged GET
  (`Range: bytes=0-0`), which returns a `Content-Range:
  0-0/<total>` header giving the exact size --
  `gdcfetch.get_data_size` does this (used by
  `check_updates()`'s `tcga_atac_bigwigs` entry, method
  `gdc_blob_size`). Note it only gives size, not md5 -- the
  `Content-MD5` header on a ranged response covers just the
  requested byte range, not the whole file.
- **`pyBigWig` cannot open bigWigs from a URL** -- neither the
  ENCODE portal's `@@download` link nor the S3 URL it redirects to
  worked when tested directly. Always download to disk first (the
  `download.py` layer already does this); don't try to stream a
  bigWig for a quick check.
- **GDC occasionally 500s on an otherwise-healthy batch.** The
  download layer retries each file a few times with backoff and
  only reports the files that never succeeded after retries -- a
  single transient failure won't lose an otherwise-complete
  download. If `download_gdc_files` does raise, rerunning is safe
  and cheap: already-downloaded files are skipped.
- **Assembly bookkeeping isn't optional bookkeeping.** BRCA mixes
  hg19 (Roadmap, MCF-7 replication timing) and hg38 (ATAC, GDC gene
  expression) -- genes annotated in only one GENCODE version get
  NaN in blocks built from the other. This looks alarming in
  `validate`'s NaN-fraction check but is expected, not a bug (COAD
  already has a smaller version of this same pattern from its
  Roadmap block). Know which blocks are which assembly *before*
  running validate, so a high NaN fraction doesn't read as a
  surprise.
