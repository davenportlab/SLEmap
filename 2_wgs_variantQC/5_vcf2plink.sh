#!/bin/bash 

MAF=0.05
MISS=0.95
HWE=0.000001
max_maf=0.95

DIR_MAIN="/path/wgs/01_5.variants_for_eQTLs"
DIR_OUT_1=${DIR_MAIN}"/2.eQTL_hwe_0.000001_maf_0.05_miss_0.95"
DIR_OUT_2=${DIR_MAIN}"/3.eQTL_hwe_0.000001_maf_0.05_miss_0.95.plink"

VCF_OUT_allChrs="/281_samples.allChrs.hwe_"${HWE}"_maf_"${MAF}"_miss_"${MISS}

VCF_INPUT=${DIR_OUT_1}"/"${VCF_OUT_allChrs}".vcf.gz"
PLINK_OUTPUT=${DIR_OUT_2}"/"${VCF_OUT_allChrs}

bcftools query -l ${VCF_INPUT} > ${DIR_OUT_2}"/sample_names.txt"
awk '{print $1, $1, 2}' ${DIR_OUT_2}"/sample_names.txt" > ${DIR_OUT_2}"/sample_sex_update.txt"


## make bed/bim/fam
plink2 --vcf ${VCF_INPUT} \
	--update-sex ${DIR_OUT_2}"/sample_sex_update.txt" \
	--split-par hg38 \
	--make-bed \
	--set-missing-var-ids @:# \
	--max-alleles 2 \
	--out ${PLINK_OUTPUT}

