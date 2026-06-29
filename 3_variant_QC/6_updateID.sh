#!/bin/bash

# Define input file
input_file=/path/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95.vcf.gz
output_file=/path/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz

bsub -o updateID.o -e updateID.e -M 100000 -n 1 -R 'select[mem>100000] rusage[mem=100000]' "bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -O z -o $output_file $input_file"
