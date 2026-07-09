#step 0. set the directory
cd path/

#step 1. fixing palindrome 
plink2 --pfile ../oneK1K_genotype_QC1_subset --exclude-palindromic-snps --make-bed --out oneK1K_genotype_QC1_subset_nopalindromic
plink2 --pfile oneK1K_genotype_QC1_subset_nopalindromic --make-bed --out oneK1K_genotype_QC1_subset_nopalindromic

##step 2. removed the unwanted snps (allele swap and strand flip)
#run the python script to list the snps that needed to be removed (2.fix_input_imputation.py)
#remove using plink
plink --bfile oneK1K_genotype_QC1_subset_nopalindromic \
      --flip Strand-Flip.txt \
      --make-bed \
      --out TEMP

plink --bfile TEMP \
      --a2-allele Force-Allele.txt \
      --make-bed \
      --out oneK1K_genotype_QC1_subset_fixed

#step 3. saved per chromosome (format: vcf version 4.2 since the server only accept this version - not yet accept 4.3)
for chr in {1..22} X; do
	plink2 --bfile oneK1K_genotype_QC1_subset_fixed \
				--chr $chr \
	      --export vcf-4.2 bgz\
	      --out oneK1K_genotype_QC1_subset_fixed_TOPMed-chr${chr}
done

## to check the version of vcf file, you cn run:
zcat oneK1K_genotype_QC1_subset_fixed_TOPMed-chr1.vcf.gz | head -n 5
