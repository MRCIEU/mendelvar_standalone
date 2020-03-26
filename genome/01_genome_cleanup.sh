#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/genome
scripts=$HOME/bin/MendelVar/genome
giggle=$HOME/bin/eczema_gwas_fu/annotation/giggle

mkdir -p $analysis

cd $analysis

#Genome annotation for version hg38  - matching APPRIS release
curl -O ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
curl -O ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz
#Genome annotation for version hg37 - matching APPRIS release
curl -O ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
curl -O ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

#HUGO gene set
curl -o hgnc_complete_set.txt "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"

#APPRIS isoforms for version hg19 (Matched with Gencode version 19)
curl -O http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh37/appris_data.principal.txt
mv appris_data.principal.txt hg19_appris_data.principal.txt
#APPRIS isoforms for version hg38 (Matched with Gencode version 31)
curl -O http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt
mv appris_data.principal.txt hg38_appris_data.principal.txt

for a in *gz; do gunzip $a; done
#Downloaded USCS cytoband annotation in bed format from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
#Annotate with cytobands
sed 's/chr//g' ucsc_hg19_cytoBand >ucsc_hg19_cytoBand.bed

#Downloaded USCS cytoband annotation in bed format from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
#Annotate with cytobands
sed 's/chr//g' ucsc_hg38_cytoBand >ucsc_hg38_cytoBand.bed

#Convert to 1-indexed notaton to use in giggle, inclusive of last position! - ie. GFF3-style.
#Needs to correct for inclusion of scaffolds with less than 5 fields in the row.
#Add an additional column with full chromosome location

cat ucsc_hg38_cytoBand.bed | awk -v OFS="\t" '(NF == 5) {print $1, $2+1, $3, $4, $5, $1$4}' >ucsc_hg38_cytoBand_final.bed
cat ucsc_hg19_cytoBand.bed | awk -v OFS="\t" '(NF == 5) {print $1, $2+1, $3, $4, $5, $1$4}' >ucsc_hg19_cytoBand_final.bed

#Test in giggle
a="ucsc_hg38_cytoBand_final.bed"
qsub -v input=$a,output=${a}_index $giggle/sub_giggle_index.sh 
a="ucsc_hg19_cytoBand_final.bed"
qsub -v input=$a,output=${a}_index $giggle/sub_giggle_index.sh 


bgzip -c test_hg38.bed > test_hg38.bed.gz
a="ucsc_hg38_cytoBand_final.bed"
qsub -v input=${a}_index,query=test_hg38.bed.gz $giggle/sub_giggle.sh

#Parse GFF3 file. Extract: chrom, start, end, strand, transcript ID, gene ID, gene_name, hgnc_id.
#Only CDS and transcript
#Warning, the field hgnc_id not present in version 19!
#Multiple CDS defined for each gene in the file, so using start codon and stop codon to define it ourselves.

#Not all the genes have a defined CDS, for instance pseudogenes often do not. 
#Possibility of using gene annotation but decided to use APPRIS transcript annotations.

#Warning, gencode v19 does not containc hgnc_id field.
python $scripts/parse_gff3.py gencode.v19.annotation.gff3 \
gencode.v19.annotation_filtered_cds.txt gencode.v19.annotation_filtered_transcript.txt
python $scripts/parse_gff3.py gencode.v31.annotation.gff3 \
gencode.v31.annotation_filtered_cds.txt gencode.v31.annotation_filtered_transcript.txt
#Check if we get unique transcript IDS.
cat gencode.v19.annotation_filtered_transcript.txt | cut -f6  >all_trans_ids_hg19.txt
cat gencode.v19.annotation_filtered_cds.txt | cut -f6 >all_cds_ids_hg19.txt

cat gencode.v31.annotation_filtered_transcript.txt | cut -f6  >all_trans_ids_hg38.txt
cat gencode.v31.annotation_filtered_cds.txt | cut -f6 >all_cds_ids_hg38.txt

#Less than 200 duplicates present.
cat all_trans_ids_hg19.txt | sort | uniq -d >unique_trans_ids_hg19.txt
grep -F -f unique_trans_ids_hg19.txt gencode.v19.annotation_filtered_transcript.txt

cat all_trans_ids_hg38.txt | sort | uniq -d >unique_trans_ids_hg38.txt
grep -F -f unique_trans_ids_hg38.txt gencode.v31.annotation_filtered_transcript.txt
#Grep the original file for these IDs to see what they represent.
#These are genes with same names on chromosome X and Y. so that is fine.

#Only 20738 genes present in APPRIS. Around 37k more "genes" annotated in 
#Gencode - RNA genes, hypothetical genes etc, most of which, other than 40,
#do not have a defined CDS (start and stop codon).
#Including them just in case - going purely by the longest transcript.
#hg19
python $scripts/filter_gff3_apris.py gencode.v19.annotation_filtered_transcript.txt \
hg19_appris_data.principal.txt gencode.v19.annotation_filtered_transcript_appris.txt
#For CDS, lots of genes present in APPRIS missing. 
#Note, that for the same gene, different transcript id can be selected in the CDS file
#versus the transcript file in case of ties (equal APPRIS isoform score or missing in APPRIS), 
#as measuring length of transcripts versus CDS in respective files to select top hit.
python $scripts/filter_gff3_apris.py gencode.v19.annotation_filtered_cds.txt \
hg19_appris_data.principal.txt gencode.v19.annotation_filtered_cds_appris.txt

#hg38
python $scripts/filter_gff3_apris.py gencode.v31.annotation_filtered_transcript.txt \
hg38_appris_data.principal.txt gencode.v31.annotation_filtered_transcript_appris.txt
python $scripts/filter_gff3_apris.py gencode.v31.annotation_filtered_cds.txt \
hg38_appris_data.principal.txt gencode.v31.annotation_filtered_cds_appris.txt

##hg19
#All apris gene ids
cat hg19_appris_data.principal.txt | cut -f2 | sort | uniq >appris_hg19_geneids
#All filtered gene ids
cat gencode.v19.annotation_filtered_cds_appris.txt | cut -f7 | sort | uniq >gencode.v19.cds_geneids
cat gencode.v19.annotation_filtered_transcript_appris.txt | cut -f7 | sort | uniq >gencode.v19.transcript_geneids
#GeneIDS unique to Appris
comm -23 appris_hg19_geneids gencode.v19.cds_geneids
comm -23 appris_hg19_geneids gencode.v19.transcript_geneids
#GeneIDS unique to Gencode
comm -13 appris_hg19_geneids gencode.v19.cds_geneids
comm -13 appris_hg19_geneids gencode.v19.transcript_geneids
##hg38
#All apris gene ids
cat hg38_appris_data.principal.txt | cut -f2 | sort | uniq >appris_hg38_geneids
#All filtered gene ids
cat gencode.v31.annotation_filtered_cds_appris.txt | cut -f7 | sort | uniq >gencode.v31.cds_geneids
cat gencode.v31.annotation_filtered_transcript_appris.txt | cut -f7 | sort | uniq >gencode.v31.transcript_geneids
#GeneIDS unique to Appris
comm -23 appris_hg38_geneids gencode.v31.cds_geneids
comm -23 appris_hg38_geneids gencode.v31.transcript_geneids
#GeneIDS unique to Gencode
comm -13 appris_hg38_geneids gencode.v31.cds_geneids
comm -13 appris_hg38_geneids gencode.v31.transcript_geneids


#Add hgnc ids to the old genome, v19. Make sure genomic range data displays current HGNC symbols. 
python $scripts/clean_v31_names.py gencode.v31.annotation_filtered_transcript_appris.txt \
hgnc_complete_set.txt >gencode.v31.annotation_filtered_transcript_appris_corrected.txt
python $scripts/clean_v19_names.py gencode.v19.annotation_filtered_transcript_appris.txt \
gencode.v31.annotation_filtered_transcript_appris.txt hgnc_complete_set.txt \
 >gencode.v19.annotation_filtered_transcript_appris_corrected.txt


#Add individualy looked up gene entries. 3 for hg38 and 6 for hg37.
cat gencode.v31.annotation_filtered_transcript_appris_corrected.txt >gencode.v31.annotation_filtered_transcript_appris_final.txt
cat gchr38_toadd.txt >>gencode.v31.annotation_filtered_transcript_appris_final.txt
cat gencode.v19.annotation_filtered_transcript_appris_corrected.txt >gencode.v19.annotation_filtered_transcript_appris_final.txt
cat gchr37_toadd.txt >>gencode.v19.annotation_filtered_transcript_appris_final.txt

#Index the files with coordinates for genes (transcripts)
function giggle_index {
a=$1
cat $a | sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > ${a}.gz
giggle index -f -i ${a}.gz -o ${a}_index
}

a="gencode.v31.annotation_filtered_transcript_appris_final.txt"
giggle_index $a
a="gencode.v19.annotation_filtered_transcript_appris_final.txt"
giggle_index $a

