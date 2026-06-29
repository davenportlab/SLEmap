#!/bin/bash 

DIR_MAIN="/path/wgs/01_5.variants_for_eQTLs"
DIR_OUT=${DIR_MAIN}"/2.eQTL_hwe_0.000001_maf_0.05_miss_0.95"
mkdir -p ${DIR_OUT}

MAF=0.05
MISS=0.95
HWE=0.000001
max_maf=0.95
echo ${max_maf}

CHR="chr"${LSB_JOBINDEX}
echo ${CHR}

if [ ${LSB_JOBINDEX} == 23 ]; then
    CHR="chrX"
    echo ${CHR}
elif [ ${LSB_JOBINDEX} == 24 ]; then
    CHR="chrY"
    echo ${CHR}
fi

VCF_IN=${DIR_MAIN}"/1.281_samples.vcf.gz"
VCF_OUT=${DIR_OUT}"/281_samples."${CHR}".hwe_"${HWE}"_maf_"${MAF}"_miss_"${MISS}".vcf.gz"


vcftools --gzvcf $VCF_IN \
	--chr ${CHR} \
	--hwe $HWE \
	--maf $MAF \
	--max-maf $max_maf \
	--max-missing $MISS \
	--recode --stdout | gzip -c > \
	$VCF_OUT

bcftools stats ${VCF_OUT} > ${DIR_OUT}"/281_samples."${CHR}".hwe_"${HWE}"_maf_"${MAF}"_miss_"${MISS}".vcf.stats.txt"


