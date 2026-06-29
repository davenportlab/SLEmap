#!/bin/bash

DIR_wgs="/path/"
DIR_out=${DIR_wgs}"/01_5.variants_for_eQTLs"
only_for_eQTLs=${DIR_out}"/281_samples_for_eQTLs.txt"
INPUT_VCF="/path/SLEmap.286_samples.vcf.gz"

bcftools view \
	-S ${only_for_eQTLs} \
	-o ${DIR_out}"/1.281_samples.vcf.gz" \
	--threads 5 \
	${INPUT_VCF}

tabix -p vcf ${DIR_out}"/1.281_samples.vcf.gz"
bcftools stats ${DIR_out}"/1.281_samples.vcf.gz" > ${DIR_out}"/1.281_samples.vcf.stats.txt" 



