#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/ontology/pathway_commons
scripts=$HOME/bin/mendelvar_standalone/ontology/pathway_commons
ont_scripts=$HOME/bin/mendelvar_standalone/ontology
mendelvar=$HOME/MendelVar

mkdir -p $analysis

cd $analysis

curl -o PathwayCommons12.All.hgnc.gmt.gz "https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.gmt.gz"
gunzip -f PathwayCommons12.All.hgnc.gmt.gz

#Pathway Commons - use PathwayCommons11.All.hgnc.gmt, map symbol to hgnc ids (can be through synonyms). Check that all symbols have a hgnc id.
#Prepare a file ready to be used by INRICH:  gene ID (using hgnc ID), gene-set ID, gene-set name
python $scripts/parse_pathway_commons.py PathwayCommons12.All.hgnc.gmt \
../hgnc_complete_set.txt pathway_commons_inrich.txt

#Make sure our input file has only unique lines.
cat pathway_commons_inrich.txt | sort | uniq >temp
mv temp pathway_commons_inrich.txt

#Subset the INRICH annotation only to those genes which are present in our final gene-disease association table.
Rscript $ont_scripts/subset_inrich.R $mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt \
pathway_commons_inrich.txt pathway_commons_inrich_subset.txt

