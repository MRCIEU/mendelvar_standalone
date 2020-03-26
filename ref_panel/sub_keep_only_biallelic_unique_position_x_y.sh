#!/bin/bash
x=ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
y=ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz

for my_file in $x $y
do
input_CEU=${my_file%.vcf.gz}_CEU_no_mono.vcf.gz
input_EUR=${my_file%.vcf.gz}_EUR_no_mono.vcf.gz
input_AFR=${my_file%.vcf.gz}_AFR_no_mono.vcf.gz
input_AMR=${my_file%.vcf.gz}_AMR_no_mono.vcf.gz
input_EAS=${my_file%.vcf.gz}_EAS_no_mono.vcf.gz
input_SAS=${my_file%.vcf.gz}_SAS_no_mono.vcf.gz
output_CEU=${my_file%_uniq.vcf.gz}_CEU_no_mono_uniq.vcf.gz
output_EUR=${my_file%_uniq.vcf.gz}_EUR_no_mono_uniq.vcf.gz
output_AFR=${my_file%_uniq.vcf.gz}_AFR_no_mono_uniq.vcf.gz
output_AMR=${my_file%_uniq.vcf.gz}_AMR_no_mono_uniq.vcf.gz
output_EAS=${my_file%_uniq.vcf.gz}_EAS_no_mono_uniq.vcf.gz
output_SAS=${my_file%_uniq.vcf.gz}_SAS_no_mono_uniq.vcf.gz
bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_CEU -o test.gz
vcftools --gzvcf test.gz --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_CEU

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_EUR -o test.gz
vcftools --gzvcf test.gz --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_EUR

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_AFR -o test.gz
vcftools --gzvcf test.gz --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_AFR

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_AMR -o test.gz
vcftools --gzvcf test.gz --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_AMR

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_EAS -o test.gz
vcftools --gzvcf test.gz --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_EAS

bcftools norm -c wx -d all -f hs37d5.fa -m +any -O z $input_SAS -o test.gz
vcftools --gzvcf test.gz --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c >$output_SAS
done




