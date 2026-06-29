#!/bin/bash 

MAF=0.05
MISS=0.95
HWE=0.000001
max_maf=0.95

DIR_MAIN="/path/wgs/01_5.variants_for_eQTLs"
DIR_OUT=${DIR_MAIN}"/2.eQTL_hwe_0.000001_maf_0.05_miss_0.95"
VCF_OUT_allChrs=${DIR_OUT}"/281_samples.allChrs.hwe_"${HWE}"_maf_"${MAF}"_miss_"${MISS}
echo ${VCF_OUT_allChrs}

## Calculate allele frequencies (using gzipped VCF)
vcftools --gzvcf ${VCF_OUT_allChrs}".vcf.gz" \
	--freq \
	--out ${VCF_OUT_allChrs}".vcf" 
gzip ${VCF_OUT_allChrs}".vcf.frq"
	
vcftools --gzvcf ${VCF_OUT_allChrs}".vcf.gz" \
	--freq2 \
	--out ${VCF_OUT_allChrs}".vcf"
gzip ${VCF_OUT_allChrs}".vcf.frq2"	





