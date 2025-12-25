#!/bin/bash

echo "Gene expression from GTEx..."
cd data/
mkdir -p gene_expression
cd gene_expression
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz
gunzip -f GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz
