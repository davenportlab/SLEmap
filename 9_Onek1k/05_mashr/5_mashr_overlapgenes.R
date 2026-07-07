####get ovelapping genes tested between slemap EUR and onek1k for each celltype
library(tidyr)
library(dplyr)
library(ggplot2)

####get genes tested for each cel ltype pair and write table with overlap + make df with num of overlap
celltypes <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","All")

celltypes_df <-data.frame(celltypes=celltypes,num_slemap=NA,num_onek1k=NA,num_overlap=NA)

for(i in 1:nrow(celltypes_df)){
  celltype <- celltypes_df[i,1]
  
  if(celltype == "All"){
    slemap <- read.table("/path/ancestry_coloc/EUR/ManualPCs/All/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv",header=T,fill=T)
    onek1k <- read.table("/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv",header=T,fill=T)
  }else{
    slemap<- read.table(paste0("/path/ancestry_coloc/EUR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv"),header=T,fill=T)
    onek1k <-read.table(paste0("/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv"),header=T,fill=T)
  }
  
  celltypes_df$num_slemap[i] <- nrow(slemap)
  celltypes_df$num_onek1k[i] <- nrow(onek1k)
  
  overlap_genes <- slemap$phenotype_id[slemap$phenotype_id %in% onek1k$phenotype_id]
  celltypes_df$num_overlap[i] <- length(overlap_genes)

  write.table(overlap_genes,paste0("/path/onek1k_locus_breaker_coloc/mashr/overlap_genes/",celltype,".tsv"),row.names = F,quote=F,col.names=F)
}

write.csv(celltypes_df,"/path/onek1k_locus_breaker_coloc/mashr/numoverlap_genes.csv")

