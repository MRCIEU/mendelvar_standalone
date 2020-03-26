#!/bin/bash
analysis=$HOME/MendelVar/ontology/freund
scripts=$HOME/bin/MendelVar/ontology/freund
ont_scripts=$HOME/bin/MendelVar/ontology
mendelvar=$HOME/MendelVar

mkdir -p $analysis

cd $analysis

git clone https://github.com/bogdanlab/gene_sets
#Freund et al. gene sets - keep as it is, map symbol to hgnc ids (can be through synonyms). Check that all symbols have a hgnc id.
#Prepare a file ready to be used by INRICH:  gene ID (using hgnc ID), gene-set ID, gene-set name
python $scripts/parse_freund.py gene_sets/mendelian_gene_sets ../hgnc_complete_set.txt freund_inrich.txt 

#Make sure our input file has only unique lines.
cat freund_inrich.txt | sort | uniq >temp
mv temp freund_inrich.txt

#Subset the INRICH annotation only to those genes which are present in our final gene-disease association table.
Rscript $ont_scripts/subset_inrich.R $mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt \
freund_inrich.txt freund_inrich_subset.txt
