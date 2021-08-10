#!/bin/sh

wget https://data.broadinstitute.org/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz
wget https://www.bcgsc.ca/downloads/POG570/POG570_TPM_expression.txt.gz
wget https://www.bcgsc.ca/downloads/POG570/Table_S1_Demographics.xlsx
wget https://www.bcgsc.ca/downloads/POG570/Table_S2_Treatment.xlsx
wget https://www.bcgsc.ca/downloads/POG570/POG570_small_mutations.txt.gz
wget -O CCLE_mutations.csv https://ndownloader.figshare.com/files/21521967
wget -O CCLE_gene_cn.csv https://ndownloader.figshare.com/files/27902124
wget -O primary-screen-replicate-collapsed-treatment-info.csv https://ndownloader.figshare.com/files/20237715
wget -O primary-screen-cell-line-info.csv https://ndownloader.figshare.com/files/20237718
wget -O primary-screen-replicate-collapsed-logfold-change.csv https://ndownloader.figshare.com/files/20237709
wget -O media-2.tsv https://www.medrxiv.org/content/medrxiv/early/2021/05/08/2021.04.30.21251941/DC2/embed/media-2.tsv?download=true
wget ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_25Feb20.xlsx
wget -O CRISPR_gene_effect.csv  https://ndownloader.figshare.com/files/27902226