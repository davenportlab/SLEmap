# 0. download imputation results (command is provided by the server)

#1. unzip
for chr in {1..22} X; do
unzip -P "rLC6T[R7oa&3BknQ.d" chr_${chr}.zip
done

#2. convert imputed vcf to pgen file (using plink2 per chromosome)
## for autosomal data
for chr in {1..22}; do
plink2 \
--vcf chr${chr}.dose.vcf.gz dosage=DS \
--double-id \
--make-pgen \
--out chr${chr}_imputed
done

## for chr X data (is slightly different since we need the sex information)
### create the sex information
awk 'NR>1 {print $1, $2, 2}' chrX_imputed-temporary.psam > sex.txt
#or
zcat chrX.dose.vcf.gz | \
grep "^#CHROM" | \
cut -f10- | tr '\t' '\n' | \
awk '{print $1, $1, 2}' > sex.txt
### insert the sex information and convert the file
plink2 \
--vcf chrX.dose.vcf.gz dosage=DS\
--double-id \
--update-sex sex.txt \
--make-pgen \
--out chrX_imputed

#3. merge all variants from all chromosomes 
## filter any indels (non snp)
for chr in {1..22} X; do
plink2 \
--pfile chr${chr}_imputed \
--snps-only just-acgt \
--make-pgen \
--out chr${chr}_imputed_SNPsonly
done

##change the rsID name (because multiallelic snps are splited into several biallelic snps)
#notes: previously it is not named line this, this has been changed to the correct format
for c in {2..22} X; do
echo chr${c}_imputed_SNPsonly
done > merge_list.txt

plink2 \
--pfile chr1_imputed_SNPsonly \
--pmerge-list merge_list.txt \
--set-all-var-ids 'chr@_#_$r_$a' \
--make-pgen \
--out onek1k_imputed_allchr_SNPsonly_renamed

#you can check the number of SNPs
wc -l onek1k_imputed_allchr_SNPsonly_renamed.pvar

