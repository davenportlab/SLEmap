DIR_MAIN="/path/"
DIR_single_cell=paste0(DIR_MAIN, "/single_cells")
DIR_wgs=paste0(DIR_MAIN, "/wgs")
DIR_out=paste0(DIR_wgs, "/01_5.variants_for_eQTLs/")

meta_data <- fread("/path/SLEmap_metadata.csv") %>%
  as.data.frame()

table(meta_data$SC_included)
# drop  use 
# 13  298
table(meta_data$condition)
# HC SLE 
# 11 300

only_for_eQTLs <- filter(meta_data, 
                         SC_included=="use" & WGS_ID!="" & condition=="SLE") %>% 
  as.data.frame()
dim(only_for_eQTLs)
# [1] 281  36
write.table(only_for_eQTLs$WGS_ID,
            paste0(DIR_out, "/281_samples_for_eQTLs.txt"),
            quote=F, row.names = F, col.names = F)
