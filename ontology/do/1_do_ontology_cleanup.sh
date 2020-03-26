#!/bin/bash
set -e
set -o pipefail
set -u
analysis=$HOME/MendelVar/ontology/do
scripts=$HOME/bin/mendelvar_standalone/ontology/do
mendelvar=$HOME/MendelVar

mkdir -p $analysis

cd $analysis

curl -o DO_AGR_slim.obo "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/subsets/DO_AGR_slim.obo"
curl -o doid-merged.obo "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/doid-merged.obo"

python $scripts/update_disease_table_do.py doid-merged.obo $mendelvar/disease/g_d_integrated_basic_desc.txt \
$mendelvar/disease/g_d_integrated_basic_desc_added1.txt

#Prepare a GAF file for GOATOOLS
python $scripts/parse_do_gaf.py $mendelvar/disease/g_d_integrated_basic_desc_added2.txt \
do_annotation.gaf
#Prepare a file ready to be used by INRICH:  gene ID (using hgnc ID), gene-set ID, gene-set name
Rscript $scripts/generate_inrich_input_do.R doid-merged.obo do_annotation.gaf do_inrich.txt


