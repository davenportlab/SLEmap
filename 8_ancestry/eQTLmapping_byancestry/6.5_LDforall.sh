#!/bin/bash

#1=snp
#2=chr
#3=frombp
#4=tobp

plink2 \
    --vcf /path/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --ld-snp $1 \
    --r2-unphased \
    --ld-window-r2 0 \
    --ld-window 99999 \
    --split-par 'hg38' \
    --update-sex /path/update_sex.txt \
    --out /path/coloc/0_calcLD/all/$1
