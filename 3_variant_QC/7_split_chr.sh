#!/bin/bash

# Define input file
input_vcf=/path/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz

# List of chromosomes to extract
chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

mkdir /path/updateID_byCHR

# Loop through chromosomes and extract each into its own file
for chr in $chroms; do
    bsub -e log/${chr}.e -o log/${chr}.o -R"select[mem>10000] rusage[mem=10000]" -M10000 bcftools view -r chr${chr} ${input_vcf} -O z -o /path/updateID_byCHR/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_chr${chr}.vcf.gz
done
