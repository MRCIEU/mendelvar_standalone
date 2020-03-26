#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/ref_panel
scripts=$HOME/bin/MendelVar/ref_panel
giggle=$HOME/bin/eczema_gwas_fu/annotation/giggle
liftover=$HOME/bin/liftover
giggle=$HOME/bin/eczema_gwas_fu/annotation/giggle

mkdir -p $analysis

cd $analysis

########GRCh37 data
#Download the reference genome for 1000 Genomes
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai

#Download the missing X and Y chromosomes.
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz

gunzip hs37d5.fa.gz
mv hs37d5.fa.gz.fai hs37d5.fa.fai

#Remove monorphic SNPs in CEU, AFR, AMR, EAS, SAS
for chrom in {1..22}
do
	my_input=ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 
	bash $scripts/sub_vcf_no_mono_1k.sh $my_input $chrom
done
#X and Y chromosomes
bash $scripts/sub_vcf_no_mono_1k_x_y.sh


#(3) Exclusion when a variant is more than biallelic (e.g., three alleles)
for chrom in {1..22}
do
	bash $scripts/sub_keep_only_biallelic_unique_position.sh $chrom 
done

#X and Y chromosomes
bash $scripts/sub_keep_only_biallelic_unique_position_x_y.sh

#Prepare a table of positions available for each population genome-wide. Report in BED format to be used with BEDtools. GFF-style.
populations=(AFR EUR CEU SAS EAS AMR)
for p in "${populations[@]}"
do
my_file=${p}_sites_1k.bed
rm -rf $my_file
touch $my_file
done

for p in "${populations[@]}"
do
	for a in ALL.chr*${p}*_no_mono_uniq.vcf.gz
	do
	out=${p}_sites_1k.bed 
	zcat $a | awk -v OFS="\t" '{print $1, $2, $2, $3}' | grep -v "#" | sed 's/^/chr/g' | sort -k 1,1 -k2,2n >>$out
	done
done

#Create giggle index. 
function giggle_index {
a=$1
cat $a | sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > ${a}.gz
giggle index -f -i ${a}.gz -o ${a}_index
}

for a in EUR_sites_1k.bed CEU_sites_1k.bed SAS_sites_1k.bed EAS_sites_1k.bed AMR_sites_1k.bed AFR_sites_1k.bed
do
giggle_index $a
done

#Download hotspots for recombination (hg37).
#Of course, no Y chromosome present.
mkdir hotspots and cd hotspots
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
gunzip genetic_map_HapMapII_GRCh37.tar.gz
tar -xvf genetic_map_HapMapII_GRCh37.tar
rm -rf genetic_map_GRCh37_chrX_par1.txt genetic_map_GRCh37_chrX_par2.txt

rm -rf all_recombination_hotspots_grch37.bed
touch all_recombination_hotspots_grch37.bed
for a in genetic*txt
do
tail -n +2 $a | awk -v OFS="\t" '($3 > 3){print $1, $2, $2}' >>all_recombination_hotspots_grch37.bed
done

#sortBed -i all_recombination_hotspots_grch37.bed >all_recombination_hotspots_grch37_sorted.bed
sort -k 1,1 -k2,2n all_recombination_hotspots_grch37.bed | uniq > all_recombination_hotspots_grch37_sorted.bed

#Convert to hg38. First need to convert to BED-style 0-indexed, n+1 style
cat all_recombination_hotspots_grch37_sorted.bed | awk -v OFS="\t" '{print $1, $2-1, $3}' >all_recombination_hotspots_grch37_liftover.bed
liftOver all_recombination_hotspots_grch37_liftover.bed $liftover/hg19ToHg38.over.chain all_recombination_hotspots_grch38_liftover.bed all_recombination_hotspots_grch38_unmapped.bed 
#Back to GTF-style BED
cat all_recombination_hotspots_grch38_liftover.bed | awk -v OFS="\t" '{print $1, $2+1, $3}' >all_recombination_hotspots_grch38_sorted.bed
ls all_recombination_hotspots_grch37_sorted.bed
ls all_recombination_hotspots_grch38_sorted.bed

########GRCh38 data
#All the chroms
ftp=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions

for a in {1..22}
do
wget $ftp/ALL.chr${a}_GRCh38_sites.20170504.vcf.gz
done

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrX_GRCh38_sites.20170504.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chrY_GRCh38_sites.20170504.vcf.gz

#Create a BED file for giggle GRCh38 searches. Use the already filtered GRch37 sites (for monorphic, multiples at the same position etc.)
rm -rf all_sites_GRCh38.bed
touch all_sites_GRCh38.bed
for a in {1..22}
do
my_file=ALL.chr${a}_GRCh38_sites.20170504.vcf.gz
zcat $my_file | awk -v OFS="\t" '{print $1, $2, $2, $3}' | grep -v "#" | sed 's/^/chr/g' | sort -k 1,1 -k2,2n >>all_sites_GRCh38.bed
done

for b in X Y
do
my_file=ALL.chr${b}_GRCh38_sites.20170504.vcf.gz
zcat $my_file | awk -v OFS="\t" '{print $1, $2, $2, $3}' | grep -v "#" | sed 's/^/chr/g' | sort -k 1,1 -k2,2n >>all_sites_GRCh38.bed
done


#Left-join GRCh38 sites to GRCh37 sites.
populations=(AFR EUR CEU SAS EAS AMR)
for p in "${populations[@]}"
do
python $scripts/join_grch38_to_grch37.py ${p}_sites_1k.bed ${p}_sites_GRCh38.bed
done


populations=(AFR EUR CEU SAS EAS AMR)
for p in "${populations[@]}"
do
#Get rid of header
cat ${p}_sites_GRCh38.bed | tail -n +2 >temp
mv temp ${p}_sites_GRCh38.bed
giggle_index ${p}_sites_GRCh38.bed
done

#Generate SNP map file for INRICH.
for a in CEU EUR AFR AMR SAS EAS
do
cat ${a}_sites_1k.bed | cut -f1,2 | sed 's/chr//' >${a}_sites_1k_GRCh37.inrich
tail -n +2 ${a}_sites_GRCh38.bed | cut -f1,2 | sed 's/chr//' >${a}_sites_1k_GRCh38.inrich 
done

#NOTE: For indels, 1k positions often vary off be a couple of bases from those on dbSNP. 
#For that reason, they are often excluded from LDlink. For SNV, positions match.
