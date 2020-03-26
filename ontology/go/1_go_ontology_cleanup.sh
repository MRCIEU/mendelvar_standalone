#!/bin/bash
set -e
set -o pipefail
set -u
analysis=$HOME/MendelVar/ontology/go
scripts=$HOME/bin/mendelvar_standalone/ontology/go
ont_scripts=$HOME/bin/mendelvar_standalone/ontology
mendelvar=$HOME/MendelVar

mkdir -p $analysis

cd $analysis

#Human annotation
curl -o goa_human.gaf.gz "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
curl -o goa_human.gpi.gz "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gpi.gz"
gunzip -f goa_human.gaf.gz
gunzip -f goa_human.gpi.gz

#GO ontology
curl -o go-basic.obo "http://current.geneontology.org/ontology/go-basic.obo"
curl -o goslim_generic.obo "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"

#For annotation use goa_human.gaf, as some GO entries missing from the isoform file. 
#Parse the GO annotation file to bring it in line with that of DO and HPO GAFs.
#Remove IEA annotations. Switch idenifier to HGNC. Filter on qualifier (remove NOT)
python $scripts/parse_go_gaf.py goa_human.gaf goa_human.gpi go.gaf

#Prepare a file ready to be used by INRICH: gene ID (using hgnc ID), gene-set ID, gene-set name
Rscript $scripts/generate_inrich_input_go.R go-basic.obo go.gaf go_inrich.txt

#Subset the INRICH annotation only to those genes which are present in our final gene-disease association table.
Rscript $ont_scripts/subset_inrich.R $mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt \
go_inrich.txt go_inrich_subset.txt

Rscript $ont_scripts/subset_inrich.R $mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt \
go_slim_inrich.txt go_slim_inrich_subset.txt

