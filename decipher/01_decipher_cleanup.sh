#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/decipher
scripts=$HOME/bin/mendelvar_standalone/decipher
omim=$HOME/MendelVar/omim
omim_scripts=$HOME/bin/mendelvar_standalone/omim

mkdir -p $analysis

cd $analysis

#Get decipher files.
wget -N "http://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz"
rm -rf DDG2P.csv
gunzip DDG2P.csv.gz
curl "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt" -o hgnc_complete_set.txt 

#Filter the results.
#Print only select columns.
python $scripts/parse_decipher.py --delimit csv --input DDG2P.csv >DDG2P.tsv

#The MIM disease numbers are often wrong in Decipher (do not exist).
#The genes often say not found in OMIM, but they are actually using a different synonym.
python $scripts/check_correctness_decipher.py \
--genemap ../omim/genemap2.txt \
--mimtitle ../omim/mimTitles_processed.txt \
--title ../omim/mimTitles.txt \
--genes_title ../omim/mimTitles_genes.txt \
--caret ../omim/omim_moved.txt \
--decipher DDG2P.csv \
--hgnc hgnc_complete_set.txt \
--output  DDG2P_corrected.csv

#Rerun parsing filtering on corrected file.
python $scripts/parse_decipher.py --delimit tsv --input DDG2P_corrected.csv >DDG2P_corrected.tsv

###STATS ONLY USAGE
#Filter out Decipher gene-disease associations to only those
#which contain a gene, and see which ones are already present in OMIM.

#both gene and disease mim present
awk -F "\t" -v OFS="\t" '($2 != "NA" && $4 != "NA") {print $0}' DDG2P_corrected.tsv >DDG2P_all_parsed_present_omim.txt

python $omim_scripts/remove_present_in_omim.py \
--omim $omim/all_gene_disease_relationships_omim.txt \
--other DDG2P_all_parsed_present_omim.txt \
--omim_gene gene_mim \
--omim_disease disease_mim \
--other_gene omim_gene_id  \
--other_disease omim_disease_id \
--pheno_series $omim/phenotype_series2.json >DDG2P_all_parsed_present_omim_unique.txt

#Number of genes unique to Decipher
cat DDG2P_corrected.tsv | cut -f1 | sort | uniq >decipher_genes
cat $omim/all_omim_database.txt | cut -f4 | sort | uniq >omim_genes
comm -23 decipher_genes omim_genes

#Compare HPO term sets present in DECIPHER to the ones in OMIM.
#See if any additional ones present.
python $omim_scripts/additional_hpo_terms.py \
--omim $omim/all_omim_database.txt \
--other DDG2P_disease_parsed_present_omim.txt \
--omim_hpo hpo \
--omim_disease '# omim_disease_id' \
--other_hpo hpo \
--other_disease omim_disease_id \
--pheno_series $omim/phenotype_series2.json \
--output DDG2P_hpo_parsed_present_omim_unique.txt


