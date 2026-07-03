library(tidyverse)
library(data.table)

DIR_MAIN="/path/coloc"
DIR_gwas_outputs <- "/path/coloc/inputs/gwas"

GWAS=paste0(DIR_MAIN,
            "/resources/GCST90270940/GCST90270940_hg38.txt.gz")
GWAS <- fread(GWAS)

### columns
# variant_id
# chr
# pos
# other_allele
# effect_allele
# beta
# odds_ratio
# effect_allele_frequency
# standard_error
# p_value

filtered_gwas <- GWAS %>% 
  dplyr::select(., rsid, chr_hg38, pos_hg38, `other allele`, `effect allele`, 
                freq, beta, se, pval) %>% 
  as.data.frame()
colnames(filtered_gwas) <- c("variant_id", "chr", "pos", "other_allele", "effect_allele",
                              "effect_allele_frequency", "beta", "standard_error", "p_value")

out_file=GWAS=paste0(DIR_gwas_outputs,
                "/GCST90270940/GCST90270940_for_eqtl.txt.gz")
OUT_FILE=gzfile(out_file )
write.table(filtered_gwas, OUT_FILE, sep="\t", quote=F, row.names = F)


