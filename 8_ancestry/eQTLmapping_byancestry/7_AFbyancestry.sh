#!/bin/bash

plink2 --vcf /path/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz --freq --out /path/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID
plink2 --vcf /path/ancestry_coloc/AFR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_AFR_n102_ancestrymaf_0.05.vcf.gz --freq --out /path/ancestry_coloc/AFR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_AFR_n102_ancestrymaf_0.05
plink2 --vcf /path/ancestry_coloc/EUR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_EUR_n81_ancestrymaf_0.05.vcf.gz --freq --out /path/ancestry_coloc/EUR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_EUR_n81_ancestrymaf_0.05
plink2 --vcf /path/ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62_ancestrymaf_0.05.vcf.gz --freq --out /path/ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62_ancestrymaf_0.05