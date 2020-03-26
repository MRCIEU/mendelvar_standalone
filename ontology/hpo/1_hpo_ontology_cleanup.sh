#!/bin/bash
set -e
set -o pipefail
set -u
analysis=$HOME/MendelVar/ontology/hpo
scripts=$HOME/bin/mendelvar_standalone/ontology/hpo
mendelvar=$HOME/MendelVar

mkdir -p $analysis

cd $analysis

#HPO ontology
curl -o hp.obo "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo"
#HPO annotation
curl -o phenotype.hpoa "http://compbio.charite.de/jenkins/job/hpo.annotations.current/lastSuccessfulBuild/artifact/current/phenotype.hpoa"
curl -o ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt "http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/phenotype_to_genes.txt"

#Need to propagate parent terms for enrichment analysis with a tool, annotation table printed to user - just leaf HPO terms.
python $scripts/update_disease_table_hpo.py phenotype.hpoa $mendelvar/disease/g_d_integrated_basic_desc_added1.txt \
$mendelvar/disease/g_d_integrated_basic_desc_added2.txt

#Remove certain ontologies (like mode of inheritance, mortality/aging, frequency, clinical modifier and onset and clinical course trees from HPO – for sure! 
#To remove children of the following terms:
#0000005 Mode of inheritance,
#0031797, Clinical course
#0040279, Frequency
#0012823, Clinical modifier

#I.e. keep only phenotypic abnormality. Will need to remove those terms from our annotation too.
Rscript $scripts/hpo_terms_to_prune.R hp.obo $mendelvar/disease/g_d_integrated_basic_desc_added2.txt \
$mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt

#Prepare a GAF file for GOATOOLS
python $scripts/parse_hpo_gaf.py $mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt \
hpo_annotation.gaf
#Prepare a file ready to be used by INRICH:  gene ID (using hgnc ID), gene-set ID, gene-set name
Rscript $scripts/generate_inrich_input_hpo.R hp.obo hpo_annotation.gaf hpo_inrich.txt

#Prepare a script to find direct descendants of a given list of terms in a file. 
#Here using it to find all the 25 direct descendants of the Phenotypic abnormality term 
#to be used for HPO slim.
python $scripts/find_direct_descendants.py hpo_inrich.txt hp.obo hpo_slim.ids hpo_slim_inrich.txt

