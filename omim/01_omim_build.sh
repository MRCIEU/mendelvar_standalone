#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/omim
scripts=$HOME/bin/mendelvar_standalone/omim

mkdir -p $analysis && cd $analysis

#Download initial data from OMIM. apiKey in the 3 links below needs to be substituted to personal API key.
curl https://data.omim.org/downloads/apiKey/genemap2.txt -o genemap2.txt
curl https://data.omim.org/downloads/apiKey/mimTitles.txt -o mimTitles.txt
curl https://data.omim.org/downloads/apiKey/morbidmap.txt -o morbidmap.txt

#This produces a file where each described phenotype, even if it is mapped to the same phenotype,
#is given a separate line.
python $scripts/parse_genemap.py genemap2.txt >genemap_parsed_master.txt

#Take a list of omim phenotype identifiers, download relevant data and check them for final inclusion in our database.
#Flags: update, prefix, no_regex_filtering
python $scripts/omim_api.py genemap_parsed_master.txt 0 initial_update_r 1

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
(cat initial_update_r_gwas_to_include_omim.txt; tail -n +2 mimtitles_parsed_r_gwas_to_include_omim.txt) >all_omim_database.txt


#Need to download JSON files for remaining disorders which have gene-disease relationship, 
#but are not of mapping key 3.
#To be able to generate comprehensive phenotypic series mappings. Will result in download 
#of additional 145 files.
tail -n +2 all_gene_disease_relationships_omim.txt | cut -f4 | sort | uniq >all_disease_mim_with_gene.txt
ls -1 json/*json | cut -d"/" -f2 | cut -d"." -f1  | sort | uniq >present_mim_files.txt
comm -23 all_disease_mim_with_gene.txt present_mim_files.txt >to_download_json.txt
python $scripts/parse_genemap_no_mk.py genemap2.txt >all_genemap_parsed.txt

(head -1 all_genemap_parsed.txt; grep -w -F -f to_download_json.txt all_genemap_parsed.txt) >additional_json_to_download.txt
python $scripts/omim_api.py additional_json_to_download.txt 0 additional_json 1
#Create a dictionary with disease MIM id - phenotypicSeriesNumber and a dict with phenotypicSeriesNumber - all the associated disease IDs.
python $scripts/get_phenotypic_series.py 

#Extract OMIM disease descriptions from the Text Description field. Around 3904 descriptions.
python $scripts/get_disease_descriptions.py mimTitles_processed.txt >omim_disease_descriptions.txt


####################
# STATS USAGE ONLY
####################

#Find all disease OMIM IDs with annotated HPO term
cat all_omim_database.txt | awk -F "\t" '($7 != "NA") {print}' | cut -f1 | sort | uniq >omim_disease_ids_hpo.txt
#Find all disease OMIM IDs with annotated DP term
cat all_omim_database.txt | awk -F "\t" '($8 != "NA") {print}' | cut -f1 | sort | uniq >omim_disease_ids_do.txt
