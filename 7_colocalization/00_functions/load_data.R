# GWAS data

DIR_tensorQTL="/path/eQTLresults"
DIR_MAIN="/path/coloc/"
DIR_plink="/path"

####################################################################
## this is an output from '01_1.eQTL_data_preparations.R'
####################################################################
print("loading 'bim_and_maf'")
out_file=paste0(DIR_MAIN,
                "/02_1.load_genotypes_data/281_samples.bim_and_maf.txt.gz")

bim_and_maf <- fread(out_file) %>%
  as.data.frame()
dim(bim_and_maf)

bim_and_maf <- bim_and_maf %>%
  dplyr::mutate(chr_pos = paste0(CHROM, "_", POS) ) %>%
  dplyr::mutate(chr_pos_variants = paste0(CHROM, "_", POS, "_", GTF_ALT, "_", GTF_REF) ) 








