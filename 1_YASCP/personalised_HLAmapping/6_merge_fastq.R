library(data.table)
library(dplyr)

source("/path/single_cells/scripts/data_load.R")



DIR_demulti <- "/path/singlecell/yascp/slemap"
donor_wgs_f <- paste0(DIR_demulti, "/gtmatch/assignments_all_pools.tsv")
donor_wgs <- fread(donor_wgs_f) %>% as.data.frame() %>%
  mutate(pool_donor=paste0(pool, ".", donor_query)) %>% 
  select(., donor_gt, final_panel, pool_donor, pool) %>%
  filter(., donor_gt!="NONE")
dim(donor_wgs)

for(i in 51:nrow(donor_wgs)){
  print(i)
    dir_pool_donor_fastq <- paste0(DIR_demulti, 
                                   "/handover/Donor_Quantification/", donor_wgs$pool[i], 
                                   "/", donor_wgs$pool_donor[i], "_fastq")
    R1 <- list.files(path=dir_pool_donor_fastq,
                     pattern="R1",
                     recursive = T,
                     full.names = T)
    R2 <- gsub("R1", "R2", R1)
    
    DIR_OUT=paste0(DIR_pers_fastq, "/", donor_wgs$donor_gt[i], "_", donor_wgs$pool[i])
    if (!dir.exists(DIR_OUT)) {dir.create(DIR_OUT)}
    
    
    merge_R1 <- paste0("cat ", 
                       paste0(R1, collapse = " "),
                       "  > ", 
                       DIR_OUT, 
                       "/", donor_wgs$donor_gt[i], "_", donor_wgs$pool[i], "_S1_L001_R1_001.fastq.gz" )
    system(merge_R1)
    merge_R2 <- paste0("cat ", 
                       paste0(R2, collapse = " "),
                       "  > ", 
                       DIR_OUT, 
                       "/", donor_wgs$donor_gt[i], "_", donor_wgs$pool[i], "_S1_L001_R2_001.fastq.gz" )
    system(merge_R2)
}
