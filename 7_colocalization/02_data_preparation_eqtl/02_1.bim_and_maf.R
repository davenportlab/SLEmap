library(tidyverse)
library(data.table)

DIR_MAIN="/path/coloc/"
DIR_plink="/path/"

geno_bim_f=paste0(DIR_plink, "/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95.bim")
geno_maf_f="/path/coloc/02_1.load_genotypes_data/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95.allele_freqs.txt.gz"

###########################################################################
## 1. merge bim and maf
###########################################################################

print("loading bim: geno_bim")
geno_bim <- fread( geno_bim_f ) %>% as.data.frame()
colnames(geno_bim) <- c("chr", "snp", "cM", "pos", "GTF_REF", "GTF_ALT")
geno_bim$chr_corrected <- ifelse(geno_bim$chr=="PAR1" |  geno_bim$chr=="PAR2",
                                 "X", geno_bim$chr)
geno_bim$tmp_id <- paste0("chr", geno_bim$chr_corrected, "_", geno_bim$pos, "_", 
                          geno_bim$GTF_ALT, "_", geno_bim$GTF_REF)
dim(geno_bim)


print("loading maf: geno_maf")
geno_maf <- fread(geno_maf_f) %>%
  as.data.frame() 
dim(geno_maf)

print("loading maf: geno_maf2")
geno_maf2 <- fread(geno_maf_f2) %>%
  as.data.frame()
colnames(geno_maf2) <- c( colnames(geno_maf2)[2:5], "ALLELE_1:FREQ", "ALLELE_2:FREQ") 
geno_maf2 <- geno_maf2 %>%  
  mutate(., tmp_id=paste0(CHROM, "_", POS, "_",
                          gsub("(.*)\\:.*", "\\1", `ALLELE_1:FREQ`),
                          "_",
                          gsub("(.*)\\:.*", "\\1", `ALLELE_2:FREQ`)) )
dim(geno_maf2)

print("merging bim and maf")
geno_merged <- merge(geno_bim, 
                     geno_maf, by="tmp_id") 
dim(geno_merged)
all(geno_merged$snp.x == geno_merged$snp.y)
colnames(geno_merged)[3] <- "snp"

out_file=paste0(DIR_MAIN,
                "/02_1.load_genotypes_data/281_samples.bim_and_maf.txt.gz")
OUT_FILE=gzfile(out_file)
write.table(geno_merged[, c(3:4, 6:7, 9:18)], 
            OUT_FILE,
            quote=F, row.names=F, sep="\t")


###########################################################################
## 2. add MAF
###########################################################################

maf_vcf<- fread(paste0(DIR_plink,
                       "/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95.vcf.frq.gz") )
maf_vcf <- maf_vcf %>% as.data.frame()
colnames(maf_vcf) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "A1_P", "A2_P")

maf_vcf_2 <- maf_vcf %>% 
  mutate(., A1=gsub("(.*)\\:(.*)", "\\1", A1_P)) %>%
  mutate(., A2=gsub("(.*)\\:(.*)", "\\1", A2_P)) %>%
  mutate(., A1_P=gsub("(.*)\\:(.*)", "\\2", A1_P)) %>%
  mutate(., A2_P=gsub("(.*)\\:(.*)", "\\2", A2_P))
maf_vcf_2$minor_allele <- ifelse(maf_vcf_2$A1_P > maf_vcf_2$A2_P, 
                                 maf_vcf_2$A2, maf_vcf_2$A1)
maf_vcf_2$minor_allele_frq <- ifelse(maf_vcf_2$A1_P > maf_vcf_2$A2_P, 
                                     maf_vcf_2$A2_P, maf_vcf_2$A1_P)
maf_vcf_2$tmp_id <- paste0(maf_vcf_2$CHROM, "_", maf_vcf_2$POS, 
                           "_", maf_vcf_2$A1, "_", maf_vcf_2$A2)
maf_vcf_2$n <- 1:nrow(maf_vcf_2)
dim(maf_vcf_2)

maf_vcf_3 <- merge(maf_vcf_2, unique( geno_merged[,c(1,3)]), by="tmp_id")
maf_vcf_3 <- maf_vcf_3 %>% arrange(., n)
dim(maf_vcf_3)

out_file2="/path/coloc/02_1.load_genotypes_data/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95.allele_freqs.txt.gz"
OUT_FILE2=gzfile(out_file2 )
write.table(maf_vcf_3[, c(1:11,13)], OUT_FILE2,
            sep="\t", quote=F, row.names = F)



min(maf_vcf_2$N_CHR)
max(maf_vcf_2$N_CHR)

