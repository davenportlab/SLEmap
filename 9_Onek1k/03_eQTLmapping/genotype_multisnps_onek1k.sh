#!/bin/bash 

# $1, chr
# $2, celltype
# $3, gene

VCF_IN=/path/onek1k/Genotypes/byCHR/onek1k_imputed_allchr_afterQC2_updatedID_chr${1}.vcf.gz
OUTPUT_PREFIX=/path/onek1k/Colocalisation/indep_coloc/${2}/genotypes/${3}_genotypes

# ## 1. subset variants from a vcf
# get variants for eQTL plot

vcftools --gzvcf ${VCF_IN} \
	--snps /path/onek1k/Colocalisation/indep_coloc/${2}/snps/${3}_snps.txt \
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


# datamash -W transpose <  $2".012.pos" > $2".012.pos.transpose.txt"
