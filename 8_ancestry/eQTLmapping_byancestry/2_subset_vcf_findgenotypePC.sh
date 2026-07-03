#!/bin/bash

module load HGI/softpack/groups/team282/data_QC_wgs/1

bcftools view -S /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/inputs/AFR_samples.txt -o /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_AFR_n102.vcf -O v /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/wl2/wgs/01_5.variants_for_eQTLs/2.eQTL_hwe_0.000001_maf_0.05_miss_0.95/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz
bcftools view -S /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/inputs/EUR_samples.txt -o /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_EUR_n81.vcf -O v /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/wl2/wgs/01_5.variants_for_eQTLs/2.eQTL_hwe_0.000001_maf_0.05_miss_0.95/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz
bcftools view -S /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/inputs/SAS_samples.txt -o /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62.vcf -O v /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/wl2/wgs/01_5.variants_for_eQTLs/2.eQTL_hwe_0.000001_maf_0.05_miss_0.95/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz


##SAS
plink2 \
    --vcf /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62.vcf \
    --make-bed \
    --maf 0.05 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/FindGenotypePCs/genotypes

plink2 \
    --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/FindGenotypePCs/genotypes \
    --indep-pairwise 50 5 0.2 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/FindGenotypePCs/prune

plink2 \
    --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/FindGenotypePCs/genotypes \
    --extract /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/FindGenotypePCs/prune.prune.in \
    --pca 20 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/FindGenotypePCs/pca

plink2 \
  --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/FindGenotypePCs/genotypes \
  --recode vcf bgz \
  --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62_ancestrymaf_0.05

##AFR
plink2 \
    --vcf /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_AFR_n102.vcf \
    --make-bed \
    --maf 0.05 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/FindGenotypePCs/genotypes

plink2 \
    --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/FindGenotypePCs/genotypes \
    --indep-pairwise 50 5 0.2 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/FindGenotypePCs/prune

plink2 \
    --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/FindGenotypePCs/genotypes \
    --extract /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/FindGenotypePCs/prune.prune.in \
    --pca 20 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/FindGenotypePCs/pca

plink2 \
  --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/FindGenotypePCs/genotypes \
  --recode vcf bgz \
  --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_AFR_n102_ancestrymaf_0.05

##EUR
plink2 \
    --vcf /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_EUR_n81.vcf \
    --make-bed \
    --maf 0.05 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/FindGenotypePCs/genotypes

plink2 \
    --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/FindGenotypePCs/genotypes \
    --indep-pairwise 50 5 0.2 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/FindGenotypePCs/prune

plink2 \
    --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/FindGenotypePCs/genotypes \
    --extract /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/FindGenotypePCs/prune.prune.in \
    --pca 20 \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/FindGenotypePCs/pca

plink2 \
  --bfile /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/FindGenotypePCs/genotypes \
  --recode vcf bgz \
  --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_EUR_n81_ancestrymaf_0.05
