from pathlib import Path

_HERE = Path(__file__).resolve().parent

location_covariates_data = (_HERE / "." / "data").resolve()

location_cov_gene_expression_gtex = (
    location_covariates_data
    / "gene_expression"
    / "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
)

location_gtex_tcga_mapping = (
    location_covariates_data
    / "gene_expression"
    / "gtex_tcga_mapping.json"
)

location_gencode19_gtf = (
    location_covariates_data
    / "gencode"
    / "gencode.v19.annotation.gtf.gz"
)

location_gencode38_gtf = (
    location_covariates_data
    / "gencode"
    / "gencode.v38.annotation.gtf.gz"
)
