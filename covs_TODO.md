# sigmutselcovs generalization TODO

Tracks the plan to generalize sigmutselcovs from COAD-only to any
TCGA cancer type (BRCA first — collaborators need breast mutation
rates from sigmutsel), and to make it a proper installable package
like sigmutsel (both eventually on PyPI).

Check items off in the same commit that completes them.

Key design decisions (2026-08-12):
- Covariate-matrix builder moves INTO this package;
  `coad_analysis/code/covariates.py` becomes a thin wrapper.
- Replication timing supports three source types: `mat`
  (16-fraction transposed table, HCT116), `fraction_bigwigs`
  (N≥2 ENCODE fraction tracks; BRCA uses MCF-7 UW 6-fraction,
  hg19), `wavelet` (single smoothed early/late track, for
  Gilbert-lab-only cell lines: LNCaP/PRAD, Caki2/KIRC, A549,
  NCI-H460, G401, SK-N-MC, T47D).
- MCF-7 MRT uses all six fractions with midpoint weights (equals
  the canonical UW weighted average, Hansen et al. 2010 PNAS).
- BRCA per-sample expression: all six STAR metrics (parity w/ COAD).
- OSF: fetch API implemented now; OSF project created + populated
  later. Default artifact to fetch = PCA-reduced matrix (small);
  full matrix and per-gene intermediates opt-in.
- Update checks: `check_updates()` + packaged sources.json (no CI
  cron).
- BRCA analysis itself will live in a future `brca_analysis` repo.

## Phase 0 — setup and in-package fixes

- [x] Add covs_TODO.md checklist (this file); pointer in
      mutation_rates/TODO.md Objective 3
- [x] import_gtex: default mapping_path to packaged
      gtex_tcga_mapping.json; fix docstring (`tcga_study_code`
      param does not exist)
- [x] load_or_generate_tcga_gexp_per_sample: fix dead cols_to_load
      and dropped tissue_type pass-through
- [x] annotate_rt_twidth: honor bin_size kwarg (literal 50000 bug)
- [x] pyproject.toml: add pyarrow + requests; pytest config
      (markers: network, slow); requires-python >=3.12
      (console entry point added later with cli.py in Phase 5)

## Phase 1 — project registry

- [x] registry.py: frozen dataclasses (ProjectSpec: GtexSpec,
      GexpSpec, AtacSpec, RoadmapSpec, RepliseqSpec,
      SimpleMatrixSpec); load_registry / available_projects /
      get_project / validate_registry
- [x] data/projects.json: defaults block + COAD row (reproduces
      current setup exactly) + BRCA row (ATAC UUID f1c06cd3-…,
      Roadmap E027/E028/E119, MCF-7 fraction accessions
      ENCFF001GSV/GTD/GTF/GTH/GTK/GSX early→late,
      gtex_breast_mammary_tissue)
- [x] tests: validation, COAD golden row, unknown-code error,
      defaults merging, mapping-key cross-check

## Phase 2 — paths and builder

- [x] paths.py: ProjectPaths mirroring coad_data/ layout
      byte-for-byte + covariate_matrices/ outputs;
      bigwig_files() ported identical (pattern list + per-pattern
      sort — column order depends on it)
- [x] builder.py: build_covariate_matrix(project, data_dir, ...)
      → CovariateMatrices(full, simple, tcga) NamedTuple; concat
      order identical to coad covariates.py:134-144; per-source
      skip-with-warning; CLR gating N<3; ATAC prefix from registry
      with fallback; simple/tcga from spec (degrade KeyError to
      warning); build_manifest.json
      (fraction_bigwigs/wavelet dispatch raises until Phase 3)
- [x] builder tests with stubbed loaders (assembly, skip paths,
      CLR gate, prefix, 3-tuple unpack, duplicate-column warning)
- [x] __init__.py lazy (PEP 562) exports of the new public API

## Phase 3 — replication timing source types

- [x] Extract _normalize_fraction_bins (behavior-preserving,
      unit-pinned)
- [x] load_repliseq_fractions_bins_from_bigwigs (N tracks
      early→late, 50 kb bins, autosomes; move _fetch_stat to
      covariates_utilities as fetch_bigwig_stat)
- [x] Wavelet type: single bigWig → per-gene rt_wavelet column
      (no mrt/clr)
- [x] source_type dispatch ("auto"|"mat"|"fraction_bigwigs"|
      "wavelet") through the load/generate functions;
      mrt_fraction_cols optional (default all)
- [x] Rename repli_seq_hct → repliseq_source (+ DeprecationWarning
      shim for one release)
- [x] Wire fraction_bigwigs + wavelet into builder._load_repliseq
- [x] Tests: N=2 closed form mrt = 0.25 + 0.5·f2; N=6 midpoints
      [1,3,5,7,9,11]/12; synthetic bigWigs in tmp_path;
      normalization refactor pin; subset renormalization

## Phase 4 — download layer

- [x] gdc.py: query_gdc_files (paginated, sorted), filters
      matching the checked-in COAD manifest (524 files; BRCA 1231)
- [x] write_gdc_manifest (gdc-client 5-col format) +
      write_gdc_sample_sheet (11 GDC columns, date-suffixed name);
      golden-row tests against real COAD sheet rows
- [x] download_file: skip-if-present, .part + atomic replace,
      HTTP Range resume, md5 verify
- [x] encode.py: resolve_encode_file → S3 URL (never stream
      bigWigs remotely); repliseq downloader (mat | fraction
      bigwigs | wavelet)
- [x] Roadmap downloader (eids × marks, tolerate 404 unless in
      required_marks)
- [x] ATAC tarball downloader (stream, extract only *.bw/*.bigWig
      flattened, guard traversal, delete tarball)
- [x] download_covariates orchestrator → DownloadReport; mocked
      tests

## Phase 4b — validation and column dictionary

- [x] cov_matrix_columns.csv writer in builder: one row per
      column — source, description, sample/detail, cell line or
      epigenome, assembly/GTF, units, fix_all transform (with
      pseudo-count); correct for partial builds
- [x] validate.py data-sanity tier: NaN fraction thresholds,
      no constant columns, value ranges (mrt/fractions in [0,1],
      sums≈1, TPM ≥ 0), expected column counts vs registry,
      unique version-stripped ENSG index
- [x] Biological-plausibility tier (Spearman, direction not
      magnitude): expression↔mrt negative; active marks
      (h3k4me3/h3k27ac/h3k9ac) + ATAC ↔ expression positive;
      repressive (h3k9me3/h3k27me3) negative; same mark across
      EIDs positive; housekeeping panel (GAPDH, ACTB, RPL13A,
      B2M) expressed + early + accessible
- [x] Auto-run at end of build (log report; strict=True raises);
      tests with planted violations

## Phase 5 — CLI, update checks, OSF fetch

- [x] cli.py console entry: download / build / fetch / validate /
      check-updates / download-gtex / projects
- [x] data/sources.json + check_updates() (GTEx GCT etag/length;
      GDC STAR file counts COAD=524 BRCA=1231; ATAC file meta via
      GET /files/<uuid> — HEAD on /data 400s; Roadmap head; ENCODE
      file status/md5; GEO mat head); update_file=True rewrites
      known blocks
- [x] fetch.py + data/osf.json: fetch_covariate_matrix(project,
      which="pca") default = PCA matrix; full/intermediates
      opt-in; CovariateArtifactsUnavailable(FileNotFoundError)
      with actionable message; plain requests
- [x] ensure_gtex_gct resolution order (explicit → packaged →
      XDG cache → download); stop writing into site-packages;
      setup.sh → thin wrapper
- [x] DEVELOPMENT.md: generalized workflow; fix stale
      run_pca_on_covariates location claim

## Phase 6 — coad_analysis migration (regression on gauss)

- [x] Regression snapshot script: current build() → parquet +
      columns JSON, BEFORE migration
- [x] covariates.py → thin build_covariate_matrix("COAD", …) call
- [x] Regression gate: identical columns (order), index, dtypes,
      np.array_equal(equal_nan=True); no cache mtime changes
- [x] Trim locations.py (keep coad_data/results/figures/
      all_maf_files; note the five unused HCT116 RT-extras for a
      future repliseq.extras registry block)
- [x] setup scripts: keep MAF step, deprecation banner for the
      covariate sections
- [x] Diff API-generated manifest vs checked-in before removing it
- [x] Update coad_analysis CLAUDE.md

## Phase 7 — BRCA bring-up (on gauss; code via scp, no GitHub auth)

- [ ] Download staged: repliseq (~52 MB) → roadmap (~8 GB, expect
      3 known 404s) → gexp (1231 files ~5 GB) → atac (~15 GB est.)
- [ ] Build BRCA; smoke: gtex_breast_mammary_tissue present;
      6 clr_rt_*; brca_* = 2×n_bigwigs; 36 roadmap columns;
      simple has 4 components
- [ ] validate_covariates("BRCA") passes both tiers
- [ ] cov_matrix_columns.csv width matches matrix; all blocks
      represented
- [ ] Registry adjustments the run reveals; document sizes

## Phase 8 — packaging (PyPI parity with sigmutsel)

- [ ] GitHub Actions test workflow (ubuntu+macos, Python 3.12,
      ruff + black + pytest)
- [ ] README rewrite (download/build/fetch workflow) + TUTORIAL.md
      (BRCA worked example) + CHANGELOG.md (start 0.2.0) +
      CONTRIBUTING.md
- [ ] Dormant publish.yml (activate after sigmutsel is on PyPI)
- [ ] Resolve run_pca_on_covariates duplication with sigmutsel

## Later (out of scope this round)

- [ ] Create OSF project; upload COAD + BRCA artifacts
      (pca + full + simple + tcga + column dictionary +
      intermediates); wire index.json
- [ ] brca_analysis repo
- [ ] More tissues via registry rows (CESC/HeLa-S3, LIHC/HepG2,
      PRAD/LNCaP wavelet, KIRC/Caki2 wavelet, LUAD/A549 wavelet,
      SKCM, …)
- [ ] repliseq.extras registry block (IZ sites, TTRs, termination
      sites, Twidth — COAD files exist, no reader wired)
