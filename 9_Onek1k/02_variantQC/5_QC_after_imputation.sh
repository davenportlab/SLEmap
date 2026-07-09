#4. QC2
#step1: Imputation quality R2 > 0.8
plink2 \
  --pfile onek1k_imputed_allchr_SNPsonly_renamed \
  --extract-if-info 'R2>0.8' \
  --make-pgen \
  --out step1_rerun
#17190177 variants remaining after main filters

#step 2: SNP call rate > 0.95 (missing <= 0.05)
plink2 --pfile step1_rerun --geno 0.05 --make-pgen --out step2_rerun
#17190177 variants remaining after main filters

#step 3: Minor Allele Frequency (MAF) > 0.05
plink2 --pfile step2_rerun --maf 0.05 --allow-extra-chr --make-pgen --out step3_rerun
#5635462 variants remaining after main filters

#step 5: Hardy–Weinberg Equilibrium (HWE) p-value > 10^-6
plink2 --pfile step3_rerun --hwe 1e-6 --make-pgen --out step4_rerun
#5635328 variants remaining after main filters

#convert to vcf
plink2 --pfile step4_rerun --set-all-var-ids 'chr@_#_$r_$a' --make-pgen --export vcf-4.2 bgz id-paste=iid --out onek1k_imputed_allchr_afterQC2_new_filter

#create index
tabix -p vcf onek1k_imputed_allchr_afterQC2_new_filter.vcf.gz

#2. LD pruning  
 plink2 \
    --pfile step4_rerun \
    --indep-pairwise 50 5 0.2 \
    --out ../find_genotype_pc/prune_rerun
# --indep-pairwise (1 compute thread): 5297231/5635328 variants removed.

#3. creating PCA
plink2 \
    --pfile step4_rerun \
    --extract ../find_genotype_pc/prune_rerun.prune.in \
    --pca 20 \
    --out ../find_genotype_pc/onek1k_allchr_pca_rerun
# 338097 variants remaining after main filters.
# Excluding 9486 variants on non-autosomes from GRM construction.