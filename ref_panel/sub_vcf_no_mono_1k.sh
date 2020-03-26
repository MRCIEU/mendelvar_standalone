#!/bin/bash

#Remove all monomorphic SNPs
my_file=$1
my_chrom=$2
table=integrated_call_samples_v3.20130502.ALL.panel
to_remove_CEU=$(awk '($2 != "CEU") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_AFR=$(awk '($3 != "AFR") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_AMR=$(awk '($3 != "AMR") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_EAS=$(awk '($3 != "EAS") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_SAS=$(awk '($3 != "SAS") {print "--remove-indv " $1}' $table | tr '\n' ' ')

output_CEU=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_CEU_no_mono.vcf.gz
output_AFR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AFR_no_mono.vcf.gz
output_AMR=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AMR_no_mono.vcf.gz
output_EAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EAS_no_mono.vcf.gz
output_SAS=ALL.chr${my_chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_SAS_no_mono.vcf.gz

vcftools --gzvcf $my_file $to_remove_CEU --mac 1 --recode  --stdout | gzip -c >$output_CEU
vcftools --gzvcf $my_file $to_remove_AFR --mac 1 --recode  --stdout | gzip -c >$output_AFR
vcftools --gzvcf $my_file $to_remove_AMR --mac 1 --recode  --stdout | gzip -c >$output_AMR
cftools --gzvcf $my_file $to_remove_EAS --mac 1 --recode  --stdout | gzip -c >$output_EAS
vcftools --gzvcf $my_file $to_remove_SAS --mac 1 --recode  --stdout | gzip -c >$output_SAS

