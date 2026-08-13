---
name: add-tissue
description: Use when adding a new TCGA cancer type to sigmutselcovs' covariate registry (data/projects.json) -- e.g. "add UCEC to sigmutselcovs", "let's do stomach cancer next", "register TGCT", "bring up ovarian covariates". Walks GTEx/GDC/Roadmap/ENCODE source discovery through to a downloaded, built, and validated covariate matrix for one TCGA project code. TEMPORARY: delete this skill once tissues_TODO.md's priority list (and any tissues added after it) are all checked off -- it exists only to carry the BRCA bring-up lessons through the next few tissues.
---

## Why this skill exists

Adding BRCA (this package's first tissue beyond COAD) took a full
research pass: two incompatible ENCODE replication-timing cohorts,
a GDC endpoint that 400s on HEAD but not GET, Roadmap tracks that
404 in ways that are normal rather than broken, and an hg19/hg38
assembly split that inflates NaN in ways that look alarming but
aren't. None of that is written down anywhere Claude would
otherwise see it. This skill carries it forward so the next few
tissues (starting with the priority list in `tissues_TODO.md`:
UCEC, STAD, TGCT, OV) don't require re-deriving it from scratch.

Read `references/lessons.md` before step 4 (Roadmap) and step 5
(replication timing) -- those are where BRCA's research took the
longest and where the gotchas live. The rest of this file is the
workflow.

## Prerequisites

- The registry schema lives in `src/sigmutselcovs/registry.py`
  (dataclasses: `GtexSpec`, `GexpSpec`, `AtacSpec`, `RoadmapSpec`,
  `RepliseqSpec`, `SimpleMatrixSpec`) -- read it once if unfamiliar.
- `data/projects.json` has COAD and BRCA rows already; both are
  templates. BRCA is the more relevant one since it's the last
  tissue actually added this way.
- `DEVELOPMENT.md` § "Adding a new cancer type" has a condensed
  version of this workflow; this skill goes deeper on sourcing.
- Downloads run on whichever machine has disk space (BRCA ran on
  `gauss`, reached via `ssh gauss` -- no GitHub auth there, so sync
  code with `scp`, not `git pull`; see project CLAUDE.md files).

## Workflow

### 1. Confirm the code and pull its GDC metadata

The TCGA code (e.g. `UCEC`) is not something to guess -- confirm it
against the live GDC project list:

```
GET https://api.gdc.cancer.gov/projects
    ?filters={"op":"=","content":{"field":"program.name","value":"TCGA"}}
    &fields=project_id,name,primary_site&size=100
```

Note the `primary_site` -- it's the anchor for steps 3-5 (matching
GDC ATAC cohorts, Roadmap epigenomes, and ENCODE cell lines to the
right tissue).

### 2. GTEx mapping

Check `src/sigmutselcovs/data/gene_expression/gtex_tcga_mapping.json`
for the code. As of the BRCA work, 31/33 TCGA codes are already
mapped (only MESO and UVM are missing) -- for any of the priority
tissues (UCEC, STAD, TGCT, OV) this step is just confirming the
entry exists, not creating one. If it's genuinely missing, pick the
matching GTEx tissue column(s) from the packaged GCT
(`sigmutselcovs download-gtex` fetches it) by name -- COAD's entry
(`Colon_Sigmoid`, `Colon_Transverse*`) shows the pattern of
including close-but-distinct GTEx subtissues.

### 3. TCGA ATAC-seq

Check `tissues_TODO.md`'s table first -- it already records which
23 of 33 TCGA types have ATAC-seq coverage (Corces et al. 2018). If
covered, the tarball UUID is on
https://gdc.cancer.gov/about-data/publications/ATACseq-AWG (one
`[TAR-GZ]` link per cancer type; the type name is in the
surrounding text, not the link itself). If not covered (OV is the
known example), set `"atac": null` in the registry row -- the
builder already handles a null ATAC block gracefully (COAD and
BRCA's tests cover this path).

### 4. Roadmap Epigenomics

Read `references/lessons.md` § Roadmap before doing this step --
the URL pattern and the "404s are normal" point matter here.

Find candidate EIDs for the primary site from the Roadmap
epigenome list, then **live-verify** each `(EID, mark)` combination
actually resolves -- don't assume a track exists just because the
epigenome does. BRCA's E027 turned out to be missing H3K27ac; this
is the rule, not the exception, so budget time for a handful of
404s per epigenome.

### 5. Replication timing -- the hard part

Read `references/lessons.md` § Replication timing in full before
starting -- this is where most of the BRCA research time went, and
the two ENCODE cohorts are easy to conflate if you don't already
know they're incompatible.

Search ENCODE's `ReplicationTimingSeries` for a cell line matching
the primary site, then also check GEO GSE137764 for a 16-fraction
MAT (currently only H1/H9/HCT116). Prefer sources in this order:
`mat` (richest) > `fraction_bigwigs` (good) > `wavelet` (single
column, no CLR fractions, but still useful) > nothing.

**It is fine to ship without replication timing.** Set
`"repliseq": null` if nothing reasonable turns up -- the builder
degrades gracefully (this is exactly the path `_skip("repliseq",
...)` in `builder.py` exists for), and a matrix missing one source
is far better than one built from a mismatched cell line.

### 6. Write the registry row

Add the row to `data/projects.json`, following the COAD/BRCA
pattern (a `defaults` block is already merged in for `roadmap`,
`gexp`, `atac`, and `simple_matrix` -- only specify what differs).
Note in a comment-equivalent place (the row's `description` field,
or your commit message) which blocks ended up hg19 vs hg38 -- this
determines expected NaN inflation later and is worth a mental note
now rather than a surprise during validation.

Validate the schema:

```bash
pytest tests/test_registry.py
```

`load_registry()` will raise a clear `ValueError` if the row is
malformed (unknown keys, bad `repliseq.type`, `gtex.mapping_key`
not in the mapping file, etc.) -- fix and retry rather than
guessing at the schema.

Optionally add a golden-row test for the new tissue mirroring
`test_coad_golden_row` / `test_brca_row` in `tests/test_registry.py`
-- this guards the row against accidental future edits, same as
the existing two.

### 7. Download, staged small to large

```bash
sigmutselcovs download <CODE> --data-dir <code>_data --which repliseq
sigmutselcovs download <CODE> --data-dir <code>_data --which roadmap
sigmutselcovs download <CODE> --data-dir <code>_data --which gexp
sigmutselcovs download <CODE> --data-dir <code>_data --which atac
```

Staging this way surfaces problems (a bad ENCODE accession, a wrong
Roadmap EID) before the largest, slowest source starts. All four
downloaders are idempotent and resumable -- rerunning after a
partial failure only refetches what's missing (see
`references/lessons.md` § GDC downloads for a specific failure mode
this protects against).

### 8. Build and validate

```bash
sigmutselcovs build <CODE> --data-dir <code>_data
sigmutselcovs validate <CODE> --data-dir <code>_data
```

`validate` runs two tiers automatically (also runs inside `build`
by default): data sanity (NaN fractions, constant columns, value
ranges, column counts against the registry) and biological
plausibility (expression vs. replication timing, active/repressive
histone marks, a housekeeping-gene panel). A single failed
direction warns -- some tissues are genuinely unusual -- but three
or more failing together usually means an assembly mismatch or a
wrong column mapping; don't dismiss that as noise.

Sanity-check column counts match what the registry predicts: ATAC
columns = 2 x (number of bigWigs downloaded); Roadmap columns =
2 x (EIDs x marks that actually resolved, i.e. minus the 404s from
step 4); CLR fraction columns = the fraction count from the
`repliseq` source (0 if `wavelet` or `null`).

### 9. Record and commit

Update `tissues_TODO.md` (check the row, fill in the registry /
repliseq / roadmap columns with what was actually found) and
`covs_TODO.md` if relevant. Commit the registry row and the
bring-up results as separate commits -- the row is a code change,
the bring-up is a verification result, and splitting them makes
`git log` legible later. Follow the `git` skill's 50/70 rule for
commit messages.

## When to stop and ask

- A tissue has no usable replication-timing source, no ATAC
  coverage, *and* thin Roadmap coverage -- at that point the matrix
  is thin enough that it's worth confirming with Jorge whether it's
  still useful before investing in the build/validate steps.
- Validation fails with 3+ inverted biological directions and the
  cause isn't an obvious assembly mismatch -- this has historically
  meant a real bug, not tissue quirk, and is worth a second pair of
  eyes before shipping the row.
