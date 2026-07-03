#!/bin/bash 

# $1, chr
# $2, start
# $3, end
# $4, snp


VCF_IN=/path/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_${1}.vcf.gz
OUTPUT_PREFIX=/path/SNPs/${4}

# ## 1. subset variants from a vcf
# get variants for eQTL plot

vcftools --gzvcf ${VCF_IN} \
	--chr ${1} \
	--from-bp ${2} \
	--to-bp ${3} \
	--recode --out ${OUTPUT_PREFIX}
gzip -c ${OUTPUT_PREFIX}".recode.vcf" > ${OUTPUT_PREFIX}".recode.vcf.gz"
rm ${OUTPUT_PREFIX}".recode.vcf"

# ## 2. convert genotypes to 012
# Genotypes are represented as 0, 1 and 2, 
# where the number represent that number of non-reference alleles. 
# Missing genotypes are represented by -1. 

vcftools --gzvcf ${OUTPUT_PREFIX}".recode.vcf.gz" \
	--012 \
	--out ${OUTPUT_PREFIX}
gzip -c ${OUTPUT_PREFIX}".012" > ${OUTPUT_PREFIX}".012.gz"
rm ${OUTPUT_PREFIX}".012"
gzip -c ${OUTPUT_PREFIX}".012.pos" > ${OUTPUT_PREFIX}".012.pos.gz"
rm ${OUTPUT_PREFIX}".012.pos"

