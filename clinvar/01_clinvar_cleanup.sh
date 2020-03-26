#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/clinvar
scripts=$HOME/bin/mendelvar_standalone/clinvar

mkdir -p $analysis

cd $analysis

curl ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.md5 -o variant_summary.txt.gz.md5
curl ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz -o variant_summary.txt.gz

#Keep only variants that contain "pathogenic", "likely pathogenic" or "risk factor" among its effects.
#Targetting only single genes and a small number of genes, no subsets of many genes.
python $scripts/parse_clinvar_tabbed.py --clinvar variant_summary.txt.gz --md5 variant_summary.txt.gz.md5 \
 --uncertain 1 --gene 1 --omim ../omim/morbidmap.txt 

#Compare allele IDS between the 3 files (hg18_clinvar.txt hg19_clinvar.txt hg38_clinvar.txt)
tail -n +2 hg19_clinvar.txt | cut -f18  | sort | uniq | sort >allele_id_hg19.txt
tail -n +2 hg38_clinvar.txt | cut -f18  | sort | uniq | sort >allele_id_hg38.txt
tail -n +2 hg18_clinvar.txt | cut -f18  | sort | uniq | sort >allele_id_hg18.txt
comm -12 allele_id_hg19.txt  allele_id_hg38.txt | wc

comm -3 allele_id_hg19.txt  allele_id_hg38.txt | wc

comm -12 allele_id_hg19.txt allele_id_hg18.txt | wc

comm -12 allele_id_hg38.txt allele_id_hg18.txt | wc


#To be indexed by Giggle
function giggle_index {
a=$1
cat $a | sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > ${a}.gz
giggle index -f -i ${a}.gz -o ${a}_index
}

a="hg38_clinvar.bed"
rm -rf ${a}_index 
giggle_index $a

a="hg19_clinvar.bed"
rm -rf ${a}_index
giggle_index $a
