"""import_tcga_gene_expression's warning for files it cannot name.

The builder reads expression once per tissue type, and each pass used
to warn "No barcode" for every file of the other tissue type -- every
file in the cohort, between the two passes. Only a file missing from
the sample sheet altogether deserves the warning.
"""

import logging

import pandas as pd

from sigmutselcovs.covariates_gene_expression import (
    import_tcga_gene_expression,
)


def _write_counts(root, file_id, tpm):
    d = root / file_id
    d.mkdir()
    pd.DataFrame(
        {
            "gene_id": ["ENSG00000000003.15", "N_unmapped"],
            "gene_name": ["TSPAN6", ""],
            "tpm_unstranded": [tpm, 0.0],
        }
    ).to_csv(
        d / f"{file_id}.rna_seq.augmented_star_gene_counts.tsv",
        sep="\t",
        index=False,
    )


def test_other_tissue_files_are_skipped_quietly(tmp_path, caplog):
    pd.DataFrame(
        {
            "File ID": ["f-tumor", "f-normal"],
            "Sample ID": ["TCGA-AA-0001-01A", "TCGA-AA-0001-11A"],
            "Tissue Type": ["Tumor", "Normal"],
        }
    ).to_csv(tmp_path / "gdc_sample_sheet.tsv", sep="\t", index=False)
    _write_counts(tmp_path, "f-tumor", 5.0)
    _write_counts(tmp_path, "f-normal", 7.0)
    _write_counts(tmp_path, "f-stray", 9.0)

    with caplog.at_level(logging.WARNING):
        out = import_tcga_gene_expression(
            tmp_path, tissue_type="Tumor"
        )

    assert set(out["Tumor_Sample_Barcode"]) == {"TCGA-AA-0001-01A"}
    warned = [r.getMessage() for r in caplog.records]
    assert not any("f-normal" in m for m in warned)
    assert any("f-stray" in m for m in warned)
