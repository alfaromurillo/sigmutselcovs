"""Gene expression covariate processing.

This module provides functions to load and process gene expression
data from GTEx (Genotype-Tissue Expression) and TCGA (The Cancer
Genome Atlas) for use as covariates in mutation rate modeling.
"""

import pandas as pd
import logging
from pathlib import Path

from covariates_utilities import annotate_with_gene_features


logger = logging.getLogger(__name__)


def import_gtex(
        loc_cov_gtex: str | Path,
        variants_df: pd.DataFrame | None = None,
        columns: str | list[str] | None = None,
        mapping_path: str | Path | None = None,
        strip_gene_id_version: bool = True) -> pd.DataFrame:
    """Import GTEx expression data.

    This function reads a GTEx **gene-level** expression matrix in GCT
    format (two header lines; gene IDs in the 'Name' column),
    restricts to `columns`, and:

    - **Gene mode**: if `variants_df` is None, returns a DataFrame indexed by
      `ensembl_gene_id` with the requested columns.
    - **Variant mode**: if `variants_df` is provided, calls
      `annotate_with_gene_features(...)` to merge the requested columns onto
      the variants table (the variants must be joinable by gene).

    Parameters
    ----------
    loc_cov_gtex : str | pathlib.Path
        Path to a GTEx GCT file, e.g. `GTEx_*_gene_tpm.gct.gz`. Must have two
        header lines and a 'Name' column containing Ensembl gene IDs.
    variants_df : pd.DataFrame or None
        If provided, a variants table compatible with
        func:`covariates_utilities.annotate_with_gene_features` (i.e.,
        it can be joined to a gene-level feature table).
    columns : str | list[str] | None
        Which GTEx tissue columns to include.
          - None (default): include all GTEx tissue columns available in
            the file (excluding 'Description').
          - list[str]: explicit GTEx column names to select.
          - str: a TCGA study code; requires `mapping_path` to be provided
            and point to a JSON mapping from codes to GTEx columns.
    mapping_path : str | pathlib.Path | None
        Optional path to a JSON mapping file that maps TCGA study
        codes (keys) to lists of GTEx column names (values). The
        mapping is inherently dependent on the specific GTEx summary
        file used (columns present), so when you update to a new GTEx
        file you should also update this mapping accordingly. If
        ``None`` and ``tcga_study_code`` is a string, a built-in
        default mapping is used for common studies.
    strip_gene_id_version : bool
        If True, strip the `.version` suffix from Ensembl IDs in the
        GTEx table (e.g., `ENSG00000123456.17 →
        ENSG00000123456`). This helps matching to variant annotations
        that typically use versionless IDs. Only relevant when
        `variants`_df is None

    Returns
    -------
    pd.DataFrame
        - **Gene mode**: DataFrame indexed by `ensembl_gene_id` with the
          selected tissue column(s).
        - **Variant mode**: A copy of `variants_df` with the selected
          GTEx column(s) attached. Columns are labeled as
          `gtex_{lowercased_column}` with non-alphanumeric chars
          replaced by underscores.

    Notes
    -----
    - This function expects **gene-level** GTEx (not sample-level per subject).
    - The helper `annotate_with_gene_features(...)` must accept a
      gene-indexed feature table (index = `ensembl_gene_id`) and a
      list of column names to attach.

    """
    gtex = pd.read_table(loc_cov_gtex,
                         skiprows=2,
                         index_col='Name')

    if isinstance(columns, list):
        cols = columns
    elif columns is None:
        # all GTEx tissue columns except 'Description'
        cols = [c for c in gtex.columns if c != 'Description']
    else:
        code = columns.upper()
        if mapping_path is None:
            raise ValueError(
                "columns is a TCGA study code but "
                "mapping_path was not provided.")
        try:
            import json
            with open(mapping_path, 'r') as fh:
                mapping = json.load(fh)
            mapping = {str(k).upper(): v for k, v in mapping.items()}
        except Exception as exc:
            raise ValueError("Failed to load mapping "
                             f"from {mapping_path}: {exc}")
        if code not in mapping:
            raise ValueError(
                f"Unsupported study code '{code}' "
                f"in mapping {mapping_path}.")
        cols = mapping[code]

    # keep only columns that exist in the GTEx file; warn if missing
    present_cols = [c for c in cols if c in gtex.columns]
    missing = [c for c in cols if c not in gtex.columns]
    if missing:
        logger.warning("GTEx missing columns: %s", ", ".join(missing))
    cols = present_cols

    def _label_for(col: str) -> str:
        return 'gtex_' + (
            col.lower()
            .replace(' ', '_')
            .replace('-', '_')
            .replace('.', '_'))

    if variants_df is None:
        out = gtex.copy()
        out = out.reset_index()
        if strip_gene_id_version:
            out['ensembl_gene_id'] = (out['Name']
                                      .astype('string')
                                      .str
                                      .replace(r'\.\d+$', '', regex=True))
        else:
            out['ensembl_gene_id'] = out['Name']
        out = out.set_index('ensembl_gene_id')
        out = out[cols]
        # rename columns to gtex_* in gene mode too
        out = out.rename(columns={c: _label_for(c) for c in cols})

    else:
        # Variant case
        labels = [_label_for(c) for c in cols]
        out = annotate_with_gene_features(
            variants_df,
            gtex,
            cols,
            labels=labels)

    return out


def import_tcga_gene_expression(
        loc_dir: str | Path,
        cols: list[str] | None = None,
        strip_gene_id_version: bool = True) -> pd.DataFrame:
    """Load TCGA STAR-count files and add Tumor_Sample_Barcode.

    Read all GDC STAR gene-count TSVs under *loc_dir*, map each file
    UUID to the 16-char TCGA barcode using the sample sheet, and return
    a long table. Rename 'gene_id' to 'ensembl_gene_id'. If
    *strip_gene_id_version* is True, remove the final '.v' from the
    Ensembl IDs.

    Parameters
    ----------
    loc_dir
        Root directory with UUID subfolders and gdc_sample_sheet*.tsv.
    cols
        Columns to keep. Defaults to:
        'ensembl_gene_id', 'gene_name', 'tpm_unstranded'
    strip_gene_id_version
        If True (default) drop the '.version' suffix from
        'ensembl_gene_id'.

    Returns
    -------
    pd.DataFrame
        Columns: 'Tumor_Sample_Barcode' plus selected gene fields.
        Rows with QC counters ('N_*') are removed.
    """
    root = Path(loc_dir)

    if cols is None:
        cols = ['ensembl_gene_id', 'gene_name', 'tpm_unstranded']
    elif cols == 'full':
        cols = ['ensembl_gene_id',
                'gene_name',
                'unstranded',
                'stranded_first',
                'stranded_second',
                'tpm_unstranded',
                'fpkm_unstranded',
                'fpkm_uq_unstranded']

    # ── load sample sheet → UUID → barcode map ─────────────────────
    ss_paths = sorted(root.glob('gdc_sample_sheet*.tsv'))
    if not ss_paths:
        raise FileNotFoundError("No gdc_sample_sheet*.tsv found.")
    ss = pd.read_csv(ss_paths[-1], sep='\t', dtype=str)

    file_id_col = 'File ID' if 'File ID' in ss.columns else 'FileID'
    if 'Sample ID' not in ss.columns:
        raise KeyError("Sample sheet missing 'Sample ID' column.")

    id_to_bar = (ss[[file_id_col, 'Sample ID']]
                 .dropna()
                 .drop_duplicates(subset=[file_id_col])
                 .set_index(file_id_col)['Sample ID']
                 .astype(str)
                 .to_dict())

    # ── find STAR count TSVs ───────────────────────────────────────
    tsvs = list(root.glob('**/*augmented_star_gene_counts.tsv'))

    frames: list[pd.DataFrame] = []
    for p in tsvs:
        file_id = p.parent.name
        barcode = id_to_bar.get(file_id)
        if barcode is None:
            logger.warning("No barcode for %s; skipping %s", file_id, p)
            continue

        df = pd.read_csv(p, sep='\t', comment='#')

        # rename 'gene_id' → 'ensembl_gene_id' before filtering
        if 'gene_id' in df.columns and \
           'ensembl_gene_id' not in df.columns:
            df = df.rename(columns={'gene_id': 'ensembl_gene_id'})

        # drop QC counters
        if 'ensembl_gene_id' in df.columns:
            mask = ~df['ensembl_gene_id'].astype(str).str.startswith('N_')
            df = df[mask]

        # optionally strip Ensembl version
        if strip_gene_id_version and 'ensembl_gene_id' in df.columns:
            df['ensembl_gene_id'] = (df['ensembl_gene_id']
                                     .astype('string')
                                     .str
                                     .replace(r'\.\d+$', '', regex=True))

        # keep only requested columns that exist
        keep = [c for c in cols if c in df.columns]
        df = df[keep]

        # coerce numerics
        for c in ['unstranded', 'tpm_unstranded', 'fpkm_unstranded',
                  'fpkm_uq_unstranded']:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors='coerce')

        df.insert(0, 'Tumor_Sample_Barcode', barcode)
        frames.append(df)

    if not frames:
        logger.warning("No STAR gene-count TSVs parsed under %s", root)
        return pd.DataFrame(columns=['Tumor_Sample_Barcode'] + cols)

    out = pd.concat(frames, ignore_index=True)
    return out


def load_or_generate_mean_tcga_gexp(
        location_csv: str | Path,
        tcga_dir: str | Path,
        *,
        cols: list[str] | None = None,
        strip_gene_id_version: bool = True,
        force_generation: bool = False,
        float_format: str = "%.6g"
        ) -> pd.Series:
    """Load or generate the mean TCGA gene expression (TPM) per Ensembl gene.

    If the CSV exists at `location_db` and `force_generation` is
    False, load it. Otherwise, parse the GDC download at `tcga_dir`
    using `import_tcga_gene_expression`, compute the mean TPM across
    all samples, save to CSV, and return the result.

    Parameters
    ----------
    location_db : str | Path
        Path to the CSV file to read/write (e.g., 'tcga_mean_tpm.csv').
    tcga_dir : str | Path
        Directory containing UUID subfolders and a `gdc_sample_sheet*.tsv`.
    cols : list[str] | None
        Columns to request from STAR TSVs during generation. If None,
        defaults used by `import_tcga_gene_expression` are applied.
    strip_gene_id_version : bool
        Strip '.version' from Ensembl IDs during generation.
    force_generation : bool
        If True, regenerate even if the CSV already exists.
    float_format : str
        Format used when writing floats to CSV (default '%.6g').

    Returns
    -------
    pd.Series
        Index: `ensembl_gene_id`; values: mean TPM (`tpm_unstranded`).
    """
    location_csv = Path(location_csv)

    if location_csv.exists() and not force_generation:
        logger.info("Loading mean TCGA expression from %s", location_csv)
        tbl = pd.read_csv(location_csv, index_col=0)
        if tbl.shape[1] == 1:
            ser = tbl.iloc[:, 0]
        else:
            # if someone saved with an unexpected schema, try to pick
            # the TPM column
            col = ("tpm_unstranded"
                   if "tpm_unstranded" in tbl.columns
                   else tbl.columns[0])
            ser = tbl[col]
        ser = ser.astype(float)
        if ser.name is None:
            ser = ser.rename("tpm_unstranded")
        logger.info("... done loading mean TCGA expression.")
        return ser

    logger.info("Generating mean TCGA expression from %s", tcga_dir)
    df = import_tcga_gene_expression(
        tcga_dir,
        cols=cols,
        strip_gene_id_version=strip_gene_id_version)

    mean_ser = (
        df.groupby("ensembl_gene_id", sort=True)["tpm_unstranded"]
          .mean()
          .astype(float)
          .rename("tpm_unstranded"))

    # write as a 2-column CSV for readability
    mean_ser.to_frame().to_csv(location_csv, float_format=float_format)
    logger.info("Saved mean TCGA expression to %s", location_csv)
    logger.info("... done generating mean TCGA expression.")
    return mean_ser
