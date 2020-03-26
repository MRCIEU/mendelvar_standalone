#!/bin/bash

my_chrom=$1

input_CEU=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_CEU_no_mono.vcf.gz
input_AFR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AFR_no_mono.vcf.gz
input_AMR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AMR_no_mono.vcf.gz
input_EAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EAS_no_mono.vcf.gz
input_SAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_SAS_no_mono.vcf.gz
input_EUR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.vcf.gz

inter_CEU=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_CEU_no_mono_inter.vcf.gz
inter_AFR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AFR_no_mono_inter.vcf.gz
inter_AMR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AMR_no_mono_inter.vcf.gz
inter_EAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EAS_no_mono_inter.vcf.gz
inter_SAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_SAS_no_mono_inter.vcf.gz
inter_EUR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_inter.vcf.gz

output_CEU=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_CEU_no_mono_uniq.vcf.gz
output_AFR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AFR_no_mono_uniq.vcf.gz
output_AMR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AMR_no_mono_uniq.vcf.gz
output_EAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EAS_no_mono_uniq.vcf.gz
output_SAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_SAS_no_mono_uniq.vcf.gz
output_EUR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_uniq.vcf.gz


bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_CEU -o $inter_CEU
vcftools --gzvcf $inter_CEU --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_CEU

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_EUR -o $inter_EUR
vcftools --gzvcf $inter_EUR --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_EUR

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_AFR -o $inter_AFR
vcftools --gzvcf $inter_AFR --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_AFR

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_AMR -o $inter_AMR
vcftools --gzvcf $inter_AMR --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_AMR

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_EAS -o $inter_EAS
vcftools --gzvcf $inter_EAS --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_EAS

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_SAS -o $inter_SAS
vcftools --gzvcf $inter_SAS --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_SAS

rm -rf $inter_CEU $inter_EUR $inter_AFR $inter_AMR $inter_EAS $inter_SAS