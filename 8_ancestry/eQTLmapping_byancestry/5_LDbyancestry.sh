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
    --silent \
    --out /path/ancestry_coloc/0_calcLD/all/$1

plink2 \
    --vcf /path/ancestry_coloc/AFR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_AFR_n102_ancestrymaf_0.05.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --ld-snp $1 \
    --r2-unphased \
    --ld-window-r2 0 \
    --ld-window 99999 \
    --silent \
    --out /path/ancestry_coloc/0_calcLD/AFR/$1

plink2 \
    --vcf /path/ancestry_coloc/EUR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_EUR_n81_ancestrymaf_0.05.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --ld-snp $1 \
    --r2-unphased \
    --ld-window-r2 0 \
    --ld-window 99999 \
    --silent \
    --out /path/ancestry_coloc/0_calcLD/EUR/$1

plink2 \
    --vcf /path/ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62_ancestrymaf_0.05.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --ld-snp $1 \
    --r2-unphased \
    --ld-window-r2 0 \
    --ld-window 99999 \
    --silent \
    --out /path/ancestry_coloc/0_calcLD/SAS/$1

    