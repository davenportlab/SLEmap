#!/bin/bash
#module load HGI/softpack/groups/team282/data_QC_wgs/1

#1=snp
#2=chr
#3=frombp
#4=tobp

##example: bash /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/20_ancestry_coloc/6_LDbyancestry.sh chr11_64360951 11 63860951 64860951


plink2 \
    --vcf /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/eb30/results/oneK1K_imputation/Imputation_TOPMed_run3/plink_conversion_QC2/onek1k_imputed_allchr_afterQC2_new_filter.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --ld-snp $1 \
    --r2-unphased \
    --ld-window-r2 0 \
    --ld-window 99999 \
    --split-par 'hg38' \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Genotypes/calcLD/$1 \
|| plink2 \
    --vcf /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/eb30/results/oneK1K_imputation/Imputation_TOPMed_run3/plink_conversion_QC2/onek1k_imputed_allchr_afterQC2_new_filter.vcf.gz \
    --chr $2 \
    --from-bp $3 \
    --to-bp $4 \
    --split-par 'hg38' \
    --write-snplist \
    --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Genotypes/calcLD/$1
