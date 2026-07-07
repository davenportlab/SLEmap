library(tidyverse)
library(data.table)

mainDir=paste0("/path/Onek1k/Colocalisation/load_eqtls_data")

####################################################
## 1. load data: bim_and_maf
####################################################
freq <- read.table("/path/Onek1k/Genotypes/freq/onek1k_imputed_allchr_afterQC2_new_filter.afreq")
colnames(freq) <- c("CHROM","variant_id","REF","ALT","ALT_FREQS","OBS_CT")
freq <- freq %>%
  mutate(
    ALT_FREQS = as.numeric(ALT_FREQS),
    REF_FREQ  = 1 - ALT_FREQS,
    minor_allele = if_else(ALT_FREQS <= REF_FREQ, ALT, REF),
    minor_allele_frq = pmin(ALT_FREQS, REF_FREQ)
  )


####################################################
## 2. 15 cell types
####################################################
cell_types <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","Nonclassical_Monocytes","MAIT_and_GammaDelta_T_cells","All")


####################################################
## 3. merge conditional signals by cell types
####################################################
DIR_independent_nominal_p="/path/Onek1k/Colocalisation/indep_coloc"

for(cell_type in cell_types){
  print (cell_type)
  
  dir_nominal_p <- paste0(DIR_independent_nominal_p, "/", 
                          cell_type, "/nominal_p/" )
  files <- list.files(dir_nominal_p, pattern = "\\.csv$", full.names = TRUE)
  eQTL_INPUT_2 <- bind_rows(
    lapply(files, function(file) {
      read.csv(file) %>%
        mutate(source_file = basename(file))
    })
  ) %>%
    dplyr::select(., -source_file)
  dim(eQTL_INPUT_2)
  
  length(unique(eQTL_INPUT_2$phenotype_id))
  length(unique(eQTL_INPUT_2$independent_variant))
  
  eQTL_INPUT_2 <- eQTL_INPUT_2 %>%
    dplyr::mutate(eSNP_position = 
                    gsub("chr\\d+\\_(\\d+)\\_.*", "\\1", independent_variant)) %>%
    dplyr::rename(., variant_id_v2=variant_id) %>%
    dplyr::mutate(variant_id=gsub("(chr\\d+\\_\\d+)\\_.*", "\\1", variant_id_v2)) %>%
    separate_wider_delim(variant_id_v2,delim = "_",names = c("CHROM", "POS", "GTF_ALT", "GTF_REF"),cols_remove = FALSE) %>% #REF ALT is flipped on purpose due to legacy code. This is flipped back before interpretation. 
    merge(., freq[,c("variant_id","OBS_CT","minor_allele","minor_allele_frq")], by.x="variant_id_v2",by.y="variant_id")
  
  dim(eQTL_INPUT_2)
  head(eQTL_INPUT_2)
  
  OUT_FILE=gzfile( paste0( mainDir, "/", 
                           cell_type, 
                           "/cis_nominal_", cell_type, ".multiple_eSNPs.txt.gz")  )
  write.table(eQTL_INPUT_2, OUT_FILE,
              quote=F, row.names=F, sep="\t")
}
