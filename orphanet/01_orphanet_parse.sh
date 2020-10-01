#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/orphanet
scripts=$HOME/bin/MendelVar_production/orphanet
omim_scripts=$HOME/bin/MendelVar_production/omim
omim=$HOME/MendelVar/omim

cd $analysis

#Rare diseases and cross referencing
curl -O http://www.orphadata.org/data/xml/en_product1.xml

#LINEARISATION OF DISORDERS
curl -O http://www.orphadata.org/data/xml/en_product7.xml

#RARE DISEASES WITH THEIR ASSOCIATED GENES
curl -O http://www.orphadata.org/data/xml/en_product6.xml

curl -O http://www.orphadata.org/data/xml/en_product4.xml


#Disease info
python $scripts/orphanet_xml1_bs.py >orphanet_xml1_parsed
#Gene-disease relationships
python $scripts/orphanet_xml6_bs.py >orphanet_xml6_parsed
#Disease-HPO associations
python $scripts/orphanet_xml4_bs.py >orphanet_xml4_parsed

python $scripts/join_xml_sources.py

python $scripts/check_correctness_orphanet.py \
--genemap ../omim/genemap2.txt \
--title ../omim/mimTitles_processed.txt \
--genes_title ../omim/mimTitles_genes.txt \
--caret ../omim/omim_moved.txt \
--orphanet orphanet_all_parsed.txt \
--hgnc ../decipher/hgnc_complete_set.txt \
--output orphanet_all_parsed_corrected.txt

###STATS ONLY USAGE
#Filter out Orphanet gene-disease associations to only those
#which contain a gene, and see which ones are already present in OMIM.

#both gene and disease mim present
awk -F "\t" -v OFS="\t" '($5 != "NA" && $2 != "NA") {print $0}' orphanet_all_parsed_corrected.txt >orphanet_all_parsed_present_omim.txt

python $omim_scripts/remove_present_in_omim.py \
--omim $omim/all_gene_disease_relationships_omim.txt \
--other orphanet_all_parsed_present_omim.txt \
--omim_gene gene_mim \
--omim_disease disease_mim \
--other_gene omim_gene_id \
--other_disease omim_disease_id \
--pheno_series $omim/phenotype_series2.json >orphanet_all_parsed_present_omim_unique.txt

#Around 1679 OMIM-independent gene-disease associations (counting those that have both gene and disease MIM), 
#202 of which have "Not yet assessed" status.
#Around 489 of unique genes (not present in OMIM) among all associations 
cat orphanet_all_parsed2_all.txt | cut -f6 | sort | uniq >orphanet_genes
cat $omim/all_omim_database.txt | cut -f5 | sort | uniq >omim_genes
comm -23 orphanet_genes omim_genes >orphanet_unique

###STATS ONLY USAGE

#disease mim present
awk -F "\t" -v OFS="\t" '($2 != "NA") {print $0}' orphanet_all_parsed_corrected.txt >orphanet_disease_parsed_present_omim.txt

python $omim_scripts/additional_hpo_terms.py \
--omim $omim/all_omim_database.txt \
--other orphanet_disease_parsed_present_omim.txt \
--omim_hpo hpo \
--omim_disease '# omim_disease_id' \
--other_hpo hpo \
--other_disease omim_disease_id \
--pheno_series $omim/phenotype_series2.json \
--output orphanet_hpo_parsed_present_omim_unique.txt

#94838 additional HPO terms spanning 4855 diseases with MIM identifiers found in Orphanet compared to OMIM.
