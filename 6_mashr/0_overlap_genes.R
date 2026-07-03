####get ovelapping genes tested between celltype pairs
library(tidyr)
library(dplyr)
library(ggplot2)

####get genes tested for each cel ltype pair and write table with overlap + make df with num of overlap
celltypes <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","DN_T_cells","Memory_B_cells","Naive_B_cells","All")

pairs <-as.data.frame(t(combn(celltypes, 2)))
colnames(pairs) <- c("celltype_1", "celltype_2")
pairs$num_overlap <- NA

for(i in 1:nrow(pairs)){
  pair <- pairs[i,]
  celltype_1 <- read.table(paste0("/path/eQTLresults/",pair[1,1],"/results/TensorQTL_eQTLS/dMean__",pair[1,1],"_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv"),header=T,fill=T)
  if(pair[1,2] == "All"){
    celltype_2 <- read.table("/path/eQTLresults/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv",header=T,fill=T)
  }else{
    celltype_2 <- read.table(paste0("/path/eQTLresults/",pair[1,2],"/results/TensorQTL_eQTLS/dMean__",pair[1,2],"_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv"),header=T,fill=T)
  }
  overlap_genes <- celltype_1$phenotype_id[celltype_1$phenotype_id %in% celltype_2$phenotype_id]
  pairs$num_overlap[i] <- length(overlap_genes)
  nonoverlap_genes <- celltype_1$phenotype_id[!celltype_1$phenotype_id %in% celltype_2$phenotype_id]
  nonoverlap_genes <- append(nonoverlap_genes,celltype_2$phenotype_id[!celltype_2$phenotype_id %in% celltype_1$phenotype_id])
  pairs$num_nonoverlap[i] <- length(nonoverlap_genes)
  
  write.table(overlap_genes,paste0("/path/mashresults/overlap_genes/",pair[1,1],"_",pair[1,2],".tsv"),row.names = F,quote=F,col.names=F)
}

pairs$allgenestested <- pairs$num_overlap + pairs$num_nonoverlap
write.csv(pairs,"/path/mashresults/overlap_genes/numoverlap_genes.csv")

