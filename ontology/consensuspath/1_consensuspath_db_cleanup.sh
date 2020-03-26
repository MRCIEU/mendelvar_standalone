#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/ontology/consensuspath
scripts=$HOME/bin/mendelvar_standalone/ontology/consensuspath
ont_scripts=$HOME/bin/mendelvar_standalone/ontology
mendelvar=$HOME/MendelVar

mkdir -p $analysis

cd $analysis

#Have to download manually from the website, using HGNC symbols.
ls CPDB_pathways_genes.tab

#ConsensusPathDB - keep as it is, map symbol to hgnc ids (can be through synonyms). Check that all symbols have a hgnc id.
#Prepare a file ready to be used by INRICH:  gene ID (using hgnc ID), gene-set ID, gene-set name
python $scripts/parse_consensuspath.py CPDB_pathways_genes.tab ../hgnc_complete_set.txt \
consensuspath_inrich.txt

#Make sure our input file has only unique lines and no spaces within columns.
cat consensuspath_inrich.txt | sort | uniq | sed 's/ /_/g' >temp
mv temp consensuspath_inrich.txt

#Subset the INRICH annotation only to those genes which are present in our final gene-disease association table.
Rscript $ont_scripts/subset_inrich.R $mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt \
consensuspath_inrich.txt consensuspath_inrich_subset.txt

