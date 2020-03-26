#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/disease
scripts=$HOME/bin/mendelvar_standalone/disease
omim=$HOME/MendelVar/omim
omim_scripts=$HOME/bin/mendelvar_standalone/omim

mkdir -p $analysis

cd $analysis
#Download gene synonyms
curl -o hgnc_complete_set.txt "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
#Download Disease descriptions
#from UniProt
curl -o UniProt_disease.txt "https://www.uniprot.org/diseases/?query=*&format=tab&limit=10000000000&columns=id"
#Basic processing to match OMIM disease ids to descriptions.
#from DO
curl -o doid.obo https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/doid.obo
python $scripts/parse_do_descriptions.py doid.obo >do_disease_descriptions.txt

#Match disease names to OMIM Identifiers.
python $scripts/parse_uniprot_descriptions.py UniProt_disease.txt \
../omim/mimTitles.txt >UniProt_disease_clean.txt

#Basic processing to extract DO ids matched to descriptions.
#OMIM
ls ../omim/omim_disease_descriptions.txt
#Orphanet (in the orphanet file)
ls ../orphanet/orphanet_all_parsed_corrected.txt

#Integrate OMIM, Orphanet, Genomics England and Decipher tables.
python $scripts/integrate_g_d_sources.py \
../omim/all_omim_database.txt \
../orphanet/orphanet_all_parsed_corrected.txt \
../decipher/DDG2P_corrected.tsv \
../genomics_england/genomics_england_latest.txt \
hgnc_complete_set.txt \
g_d_integrated.txt

#Basic filter to filter for those entries with gene-disease associations.
cat g_d_integrated.txt | awk -F"\t" -v OFS="\t" '($1 != "NA" && $7 != "NA") {print}' >g_d_integrated_basic.txt

#Add disease descriptions to master table.
python $scripts/add_disease_descriptions.py \
../omim/omim_disease_descriptions.txt \
../orphanet/orphanet_all_parsed_corrected.txt \
UniProt_disease_clean.txt \
do_disease_descriptions.txt \
g_d_integrated_basic.txt \
g_d_integrated_basic_desc.txt


##############################Stats only
#Filter the table for presence of both gene and disease (no missing gene or disease mim id)
cat g_d_integrated_basic.txt | awk -F"\t" -v OFS="\t" '($1 != "NA" && $5 != "NA") {print}' >g_d_integrated_filt.txt
python $omim_scripts/remove_present_in_omim.py \
--omim $omim/all_gene_disease_relationships_omim.txt \
--other g_d_integrated_filt.txt \
--omim_gene gene_mim \
--omim_disease disease_mim \
--other_gene omim_gene_id  \
--other_disease omim_disease_id \
--pheno_series $omim/phenotype_series2.json >g_d_integrated_filt_unique.txt

#See how many new gene-disease relationships versus OMIM alone.
#See how many new genes versus OMIM alone.
cat g_d_integrated_basic.txt | cut -f4 | sort | uniq >g_d_genes_all
cat g_d_integrated_basic.txt | cut -f2 | sort | uniq >g_d_genes_symbol
cat $omim/all_omim_database.txt | cut -f5 | sort | uniq >omim_genes
comm -23 g_d_genes_all omim_genes

#Compare with Orphanet list.
comm -23 ../orphanet/orphanet_unique g_d_genes_all
#Compare with Decipher list.
comm -23 ../decipher/decipher_genes g_d_genes_symbol
#Compare with GE list.
comm -23 ../genomics_england/ge_genes g_d_genes_all

#All the entries with OMIM disease number which have a HPO term assigned
cat g_d_integrated_basic_desc_added2.txt | awk -F "\t" '($8 != "NA") {print}' | cut -f5 | sort | uniq >g_d_hpo_omim
#All the entries with OMIM disease number which have a DO term assigned
cat g_d_integrated_basic_desc_added2.txt | awk -F "\t" '($9 != "NA") {print}' | cut -f5 | sort | uniq >g_d_do_omim

#How many HPO annotations from other sources compared to OMIM?
comm -23 g_d_hpo_omim ../omim/omim_disease_ids_hpo.txt
#How many DO annotations from other sources compared to OMIM?
comm -23 g_d_do_omim ../omim/omim_disease_ids_do.txt

#How many gene-disease associations have a description?
cat g_d_integrated_basic_desc.txt | cut -f12 | awk '($1 != "NA") {print}' | wc
cat g_d_integrated_basic_desc.txt | awk '($12 != "NA") {print}' | wc
######################################################



