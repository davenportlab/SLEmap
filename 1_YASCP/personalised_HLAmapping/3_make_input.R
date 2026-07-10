library(data.table)
library(dplyr)

source("/path/data_load.R")

out_hlala <- fread(paste0(DIR_wgs,
                          "/02_1.hla_typing/summary.HLA_LA/hlala.txt") ) %>%
  as.data.frame()
table(out_hlala$Locus)

genes <- c("A", "B", "C", "DQA1", "DQB1", "DRB1", "DPA1", "DPB1", "E", "F", "G" )
hlala <- out_hlala[out_hlala$Locus %in% genes, ]
table(hlala$Locus)

input_for_pers_ref <- data.frame(individual_ID=hlala$sample,
                                 HLA_allele=gsub("(.*)\\;.*", "\\1", hlala$Allele) )
length(unique(input_for_pers_ref$individual_ID))
# 311
write.table(input_for_pers_ref, 
            paste0(DIR_single_cell, 
                   "/02_1.personalized_mapping/input_HLA_alleles.311_samples.txt"),
            quote =F, row.names = F)

