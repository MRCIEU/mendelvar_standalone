#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/ontology/reactome
scripts=$HOME/bin/mendelvar_standalone/ontology/reactome
ont_scripts=$HOME/bin/mendelvar_standalone/ontology
mendelvar=$HOME/MendelVar

mkdir -p $analysis

cd $analysis

curl -o Ensembl2Reactome_All_Levels.txt "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt"

sed -i 's/ENSG00000183214/ENSG00000204520/g' Ensembl2Reactome_All_Levels.txt
sed -i 's/ENSG00000206312/ENSG00000204301/g' Ensembl2Reactome_All_Levels.txt
sed -i 's/ENSG00000226704/ENSG00000204390/g' Ensembl2Reactome_All_Levels.txt
sed -i 's/ENSG00000235233/ENSG00000204390/g' Ensembl2Reactome_All_Levels.txt
sed -i 's/ENSG00000278620/ENSG00000241598/g' Ensembl2Reactome_All_Levels.txt

curl -o HUMAN_9606_idmapping_selected.tab.gz "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz" 
gunzip -f HUMAN_9606_idmapping_selected.tab.gz
#Reactome - use Ensembl2Reactome_All_Levels.txt to get ENSGs, and then convert them to HGNC_IDs. Filter away IEA annotations.
#Prepare a file ready to be used by INRICH:  gene ID (using hgnc ID), gene-set ID, gene-set name
#Prepare a list of all the Ensembl IDs present in Ensembl
python $scripts/parse_reactome.py Ensembl2Reactome_All_Levels.txt \
../hgnc_complete_set.txt HUMAN_9606_idmapping_selected.tab reactome_inrich.txt | sort | uniq >missing_ensg.txt

#Make sure our input file has only unique lines.
cat reactome_inrich.txt | sort | uniq >temp
mv temp reactome_inrich.txt

#Subset the INRICH annotation only to those genes which are present in our final gene-disease association table.
Rscript $ont_scripts/subset_inrich.R $mendelvar/disease/g_d_integrated_basic_desc_added2_hpo_removed.txt \
reactome_inrich.txt reactome_inrich_subset.txt
