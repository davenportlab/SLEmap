library(tidyverse)
library(data.table)

DIR_MAIN="/path/coloc/"

####################################################
## 1. load data: bim_and_maf
####################################################
source(paste0(DIR_MAIN,
              "/scripts/00_functions/load_data.R"))
print(dim(bim_and_maf))

####################################################
## 2. 15 cell types
####################################################
cell_types <- c("CD56Bright_NK_cells", "CD56Dim_NK_cells", 
                "Classical_Monocytes", 
                "CM_CD4_T_cells", "CM_CD8_T_cells",
                "Cytotoxic_CD4_T_cells", "DN_T_cells",
                "EM_CD4_T_cells", "EM_CD8_T_cells",
                "Memory_B_cells", "Naive_B_cells",
                "Naive_CD4_T_cells", "Naive_CD8_T_cells",
                "Regulatory_CD4_T_cells", "TEMRA")

####################################################
## 3. merge conditional signals by cell types
####################################################
DIR_independent_nominal_p="/path/eQTLresults/2_indep_coloc/"
mainDir=paste0(DIR_MAIN, 
               "/02_2.load_eqtls_data/")

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
    dplyr::mutate(CHROM=gsub("(chr\\d+)\\_\\d+", "\\1", variant_id)) %>%
    dplyr::mutate(POS=gsub("(chr\\d+)\\_(\\d+)", "\\2", variant_id)) %>%
    #merge(., get_N_CHR, by="variant_id_v2") %>%   # add 'N_CHR' from get_N_CHR
    merge(., bim_and_maf[, c("chr_pos_variants", "N_CHR", 
                             "GTF_REF", "GTF_ALT", "minor_allele_frq")], 
          by.x="variant_id_v2", by.y="chr_pos_variants")
  
  dim(eQTL_INPUT_2)
  head(eQTL_INPUT_2)

  OUT_FILE=gzfile( paste0( mainDir, "/", 
                           cell_type, 
                           "/cis_nominal_", cell_type, ".multiple_eSNPs.txt.gz")  )
  write.table(eQTL_INPUT_2, OUT_FILE,
              quote=F, row.names=F, sep="\t")
}

