## code (command line) and information about the results

### PREPARATION ####
## copy data and check the number of individuals and SNPs
cd /path/onek1k/powel_provided_array_genotypes

# step 0: copy the data to the desired directory + check the number of variants and individuals
plink2 --bfile OneK1K_AllChr --make-pgen --out /path/genotype_oneK1K_plink2_v6_correctdata/step0

## code that can be used to do quick check on numbers
wc -l <name>.psam #need to -1 
wc -l <name>.pvar #need to -1 

## change the directory our directory
cd /path/genotype_oneK1K_plink2_v6_correctdata

### QC1 steps ####
#step 1: Keep only biallelic SNPs
## take the ID only
awk 'NR>1 && length($4)==1 && length($5)==1 && $4~/^[ACGT]$/ && $5~/^[ACGT]$/ {print $3}' step0.pvar > keep_snps.txt
## incase if you want to check the filtered SNPs
awk 'length($4)==1 && length($5)==1 && $4~/^[ACGT]$/ && $5~/^[ACGT]$/' step0.pvar > step0_biallelic_snps.pvar
##filter out the non biallelic SNPs
plink2 --pfile step0 --extract keep_snps.txt --make-pgen --out step1
##run this to make sure no more variant will be removed -- yes, no more
plink2 --pfile step1 --snps-only just-acgt --max-alleles 2 --make-pgen --out step1_1 

#step 2: SNP call rate > 0.97 (missing < 0.03)
plink2 --pfile step1 --geno 0.03 --make-pgen --out step2

#step 3: individual call rate > 0.97 (missing < 0.03)
plink2 --pfile step2 --mind 0.03 --make-pgen --out step3

#step 4: Minor Allele Frequency (MAF) > 0.01
plink2 --pfile step3 --maf 0.01 --allow-extra-chr --make-pgen --out step4

#step 5: Hardy–Weinberg Equilibrium (HWE) p-value > 10^-6
plink2 --pfile step4 --hwe 1e-3 --make-pgen --out step5

#step 6: Heterozygosity rate outliers (±3 SD)
plink2 --pfile step5 --het --out step6 #calculate per-individual heterozygosity
##(filtered the result in python)
plink2 --pfile step5 --remove het_outliers.txt --make-pgen --out step6_filtered

#step 7: GRM <0,125
plink2 --pfile step6_filtered --king-cutoff 0.125 --out step7
plink2 --pfile step6_filtered --remove step7.king.cutoff.out.id --make-pgen --out step7_filtered

#step 8: Sex mismatch check -- this is additional, not to filter anything just to check since SLEmap did it
plink2 --pfile step7_filtered --check-sex --out sexcheck

awk 'NR>1{
    f=$6; y=$7;
    if(f<0.2 && y<0.36) sex="F"
    else if(f>0.8 && y>0.7) sex="M"
    else sex="U"
    print $1,$2,sex
}' sexcheck.sexcheck > sex_inferred_Y.txt

#step 9: subset
##the selected donors is listed in this file: selected_donor_for_imputation.txt
plink2 --pfile step7_filtered --keep ../selected_donor_for_imputation.txt --make-pgen --out oneK1K_genotype_QC1_subset