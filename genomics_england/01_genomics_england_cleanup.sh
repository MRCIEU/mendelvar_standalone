#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/genomics_england
scripts=$HOME/bin/mendelvar_standalone/genomics_england
omim=$HOME/MendelVar/omim
omim_scripts=$HOME/bin/mendelvar_standalone/omim

mkdir -p $analysis

cd $analysis

#Parse Genomics England API.
#If no OMIM ID present, see if we can match it to OMIM by brute force.
python $scripts/genomics_england_api.py \
--genemap ../omim/genemap2.txt \
--title $omim/mimTitles.txt \
--caret $omim/omim_moved.txt \
--mimtitle $omim/mimTitles_processed.txt \
--genes_title $omim/mimTitles_genes.txt \
--hgnc ../decipher/hgnc_complete_set.txt \
--output genomics_england_latest.txt >log

#Remove the phenotypic series 268000 retinitis pigmentosa entries 
#which duplicate the more detailed ones from OMIM.
cat genomics_england_latest.txt | awk -F"\t" -v OFS="\t" '($1 != "268000") {print}' >temp
mv temp genomics_england_latest.txt 


###STATS ONLY USAGE
#Filter out Genomics England gene-disease associations to only those
#which contain a gene, and see which ones are already present in OMIM.
python $omim_scripts/remove_present_in_omim.py \
--omim $omim/all_gene_disease_relationships_omim.txt \
--other genomics_england_latest.txt \
--omim_gene gene_mim \
--omim_disease disease_mim \
--other_gene omim_gene_id  \
--other_disease '# omim_disease_id' \
--pheno_series $omim/phenotype_series2.json >genomics_england_unique.txt

#Number of genes unique to Genomics England
cat genomics_england_latest.txt | cut -f7 | sort | uniq >ge_genes
cat $omim/all_omim_database.txt | cut -f5 | sort  | uniq >omim_genes
comm -23 ge_genes omim_genes
