#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/omim
scripts=$HOME/bin/mendelvar_standalone/omim

cd $analysis

#TO RUN FOR EVERY UPDATE
#User needs to be a registered OMIM user with a personal API Key.

#Run for an update to highlight potential targets to be updated in our database.
#To process that we need both Update List from Omim website and latest edition of GeneMap to 
#filter out only disorder ids.
curl https://data.omim.org/downloads/apiKey/genemap2.txt -o genemap2.txt
curl https://data.omim.org/downloads/apiKey/mimTitles.txt -o mimTitles.txt
curl https://data.omim.org/downloads/apiKey/morbidmap.txt -o morbidmap.txt

#Download the Update List for a given month - change. Save as the "latest.html" file
wget https://www.omim.org/statistics/updates/2020/03 -O latest.html
#Apply the monthly update:
python $scripts/omim_check_for_updates.py genemap2.txt latest.html >latest_update.txt
###API key needs to be changed to personal API key in the script `omim_api.py`, line 21.
#Redownload modified JSON files
python $scripts/omim_api.py latest_update.txt 1 current_update_r 1
#Recalculate the output table
python $scripts/omim_api.py genemap_parsed_master.txt 0 current_update_r 1

#Create file with all the relationships between MIM diseases and genes (not just confirmed).
python $scripts/parse_genemap_relationships.py genemap2.txt >all_gene_disease_relationships_omim.txt

#Creates two outputs based on mimTitles.txt: 
#-File with all disease IDs currently in OMIM (as marked by Plus, Percent, Null, Number Sign)
##-File with all gene IDs currently in OMIM (as marked by Plus, Asterisk)
#-file which maps deprecated mim ids to the current ones (as marked by Caret)
python $scripts/all_disease_moved_entries_ref.py mimTitles.txt

#Create a file to be supplied to the omim_api.py file which includes all the OMIM diseases with
#no known molecular basis (i.e. Plus, Percent, Null). This is so that we can extract HPO, DO
#for these diseases which could be then annotated by genes from another database.
python $scripts/parse_genemap_unmapped_disease.py mimTitles.txt
python $scripts/omim_api.py mimtitles_parsed_master.txt 0 mimtitles_parsed_r 1

#Concatenate the files with OMIM output: diseases with mapped genes and without.
(cat current_update_r_gwas_to_include_omim.txt; tail -n +2 mimtitles_parsed_r_gwas_to_include_omim.txt) >all_omim_database.txt

#Create a dictionary with disease MIM id - phenotypicSeriesNumber and a dict with phenotypicSeriesNumber - all the associated disease IDs.
python $scripts/get_phenotypic_series.py 

#Extract OMIM disease descriptions from the Text Description field. Around 3904 descriptions.
python $scripts/get_disease_descriptions.py mimTitles_processed.txt >omim_disease_descriptions.txt


####################
# STATS USAGE ONLY
####################

#Find all disease OMIM IDs with annotated HPO term
cat all_omim_database.txt | awk -F "\t" '($7 != "NA") {print}' | cut -f1 | sort | uniq >omim_disease_ids_hpo.txt
#Find all disease OMIM IDs with annotated DO term
cat all_omim_database.txt | awk -F "\t" '($8 != "NA") {print}' | cut -f1 | sort | uniq >omim_disease_ids_do.txt
