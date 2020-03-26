#!/bin/bash

#Remove all monomorphic SNPs
table=integrated_call_samples_v3.20130502.ALL.panel
to_remove_CEU=$(awk '($2 != "CEU") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_AFR=$(awk '($3 != "AFR") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_EUR=$(awk '($3 != "EUR") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_AMR=$(awk '($3 != "AMR") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_EAS=$(awk '($3 != "EAS") {print "--remove-indv " $1}' $table | tr '\n' ' ')
to_remove_SAS=$(awk '($3 != "SAS") {print "--remove-indv " $1}' $table | tr '\n' ' ')

x=ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
y=ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz


for my_file in $x $y
do
output_CEU=${my_file%.vcf.gz}_CEU_no_mono.vcf.gz
output_EUR=${my_file%.vcf.gz}_EUR_no_mono.vcf.gz
output_AFR=${my_file%.vcf.gz}_AFR_no_mono.vcf.gz
output_AMR=${my_file%.vcf.gz}_AMR_no_mono.vcf.gz
output_EAS=${my_file%.vcf.gz}_EAS_no_mono.vcf.gz
output_SAS=${my_file%.vcf.gz}_SAS_no_mono.vcf.gz
vcftools --gzvcf $my_file $to_remove_CEU --mac 1 --recode  --stdout | gzip -c >$output_CEU
vcftools --gzvcf $my_file $to_remove_EUR --mac 1 --recode  --stdout | gzip -c >$output_EUR
vcftools --gzvcf $my_file $to_remove_AFR --mac 1 --recode  --stdout | gzip -c >$output_AFR
vcftools --gzvcf $my_file $to_remove_AMR --mac 1 --recode  --stdout | gzip -c >$output_AMR
vcftools --gzvcf $my_file $to_remove_EAS --mac 1 --recode  --stdout | gzip -c >$output_EAS
vcftools --gzvcf $my_file $to_remove_SAS --mac 1 --recode  --stdout | gzip -c >$output_SAS
done