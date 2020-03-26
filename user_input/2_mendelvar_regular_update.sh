#!/bin/bash
set -e
set -o pipefail
set -u

scripts=$HOME/bin/mendelvar_standalone
data=$HOME/MendelVar

##OMIM
#User needs to be a registered OMIM user with a personal API Key.
#MANUAL CHANGES
#sapiKey needs to be changed to personal API key in the script 01_omim_build.sh (lines 17,18,19) 
#and in the script `omim_api.py`, (line 21).
#Download latest Monthly update list page and change the date in the file on line 22
#to update the database:
$scripts/omim/01_omim_cleanup.sh
##Genomics England
$scripts/genomics_england/01_genomics_england_cleanup.sh
##DECIPHER
$scripts/decipher/01_decipher_cleanup.sh
##ORPHANET
$scripts/orphanet/01_orphanet_parse.sh
##Disease integration
$scripts/disease/01_disease_integration.sh
##ClinVar
$scripts/clinvar/01_clinvar_cleanup.sh
##Data_cleanup - ontology
$scripts/ontology/01_data_cleanup.sh
##ConsensusPathDb
#MANUAL CHANGES - need to download the CPDB_pathways_genes.tab file from its website, using HGNC symbols. Db not often updated!
#Download the following file: "Biological pathways (as defined by source databases) with their genes identified with HGNC ID accession numbers" and save as CPDB_pathways_genes.tab inside the $HOME/MendelVar/ontology/consensuspath folder
$scripts/ontology/consensuspath/1_consensuspath_db_cleanup.sh
##PathwayCommons
#MANUAL CHANGES - need to identify, download and subsitute the latest file from website. Db not often updated!
$scripts/ontology/pathway_commons/1_pathway_commons_cleanup.sh
#Disease Ontology
$scripts/ontology/do/1_do_ontology_cleanup.sh
#Gene Ontology
$scripts/ontology/go/1_go_ontology_cleanup.sh
#Human Phenotype Ontology
$scripts/ontology/hpo/1_hpo_ontology_cleanup.sh
#Reactome
$scripts/ontology/reactome/1_reactome_cleanup.sh
#INRICH
$scripts/ontology/inrich_files/create_inrich_ontology_input.sh

