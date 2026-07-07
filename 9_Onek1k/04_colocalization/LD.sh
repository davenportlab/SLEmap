#!/bin/bash

#1=snp
#2=chr
#3=frombp
#4=tobp


plink2 \
    --vcf /path/onek1k_imputed_allchr_afterQC2_new_filter.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --ld-snp $1 \
    --r2-unphased \
    --ld-window-r2 0 \
    --ld-window 99999 \
    --split-par 'hg38' \
    --out /path/Onek1k/Genotypes/calcLD/$1

plink2 \
    --vcf /path/onek1k_imputed_allchr_afterQC2_new_filter.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --split-par 'hg38' \
    --write-snplist \
    --out /path/Onek1k/Genotypes/calcLD/$1
