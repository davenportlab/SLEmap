###make eQTL result files combined across cell types

alltested <- data.frame()
for(celltype in c("CD56Bright_NK_cells", "CD56Dim_NK_cells", "Classical_Monocytes", "CM_CD4_T_cells", "EM_CD4_T_cells", "Naive_CD4_T_cells", "Regulatory_CD4_T_cells", "CM_CD8_T_cells", "EM_CD8_T_cells", "Naive_CD8_T_cells", "TEMRA", "Memory_B_cells", "Naive_B_cells", "MAIT_and_GammaDelta_T_cells", "Nonclassical_Monocytes")){
  qval <- read.table(paste0("/path/onek1k_eQTLresults/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"),fill=T,header=T)
  qval$celltype <- celltype
  alltested <- rbind(alltested,qval)
  
}

qval <- read.table("/path/onek1k_eQTLresults/ManualPCs/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv",fill=T,header=T)
qval$celltype <- "All"
alltested <- rbind(alltested,qval)

write.csv(alltested,"/path/onek1k_eQTLresults/alltested_genes_onek1k.csv",row.names=F)

write.csv(alltested[alltested$qval < 0.05,],"/path/onek1k_eQTLresults/allsig_genes_onek1k.csv",row.names=F)

## get indep results all for onek1k
allindep <- data.frame()
for(celltype in c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","Nonclassical_Monocytes","MAIT_and_GammaDelta_T_cells")){
   indepresults <- read.table(paste0("/path/onek1k_eQTLresults/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
   indepresults$celltype <- celltype
   allindep <- rbind(allindep,indepresults)
}

indepresults <- read.table("/path/onek1k_eQTLresults/ManualPCs/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",header=T)
indepresults$celltype <- "All"
allindep <- rbind(allindep,indepresults)
write.csv(allindep, "/path/onek1k_eQTLresults/Cis_eqtls_independent_allcelltypes_wbulklike.csv")


### make summary file
### plot eQTL summary
library(ggplot2)
library(scales)

#get number of donors tested for each cell type
seurat <- readRDS("/path/onek1k_input/onek1K_eQTL_input_seurat.RDS")

numall2 <- data.frame(unique(seurat$Celltype_level1))
numall2$numdonor_over20 <- NA
numall2$numcells <- NA
for(i in unique(seurat$Celltype_level1)){
  seurat_sub <- subset(seurat,subset = Celltype_level1 == i)
  numcells <- as.data.frame(table(seurat_sub$donor_id))
  numcells <- numcells[numcells$Freq >= 20,]
  numall2$numdonor_over20[numall2$unique.seurat.Celltype_level1. == i] <- nrow(numcells)
  if(nrow(numcells) != 0){
    seurat_sub2 <- subset(seurat_sub,subset = donor_id %in% numcells$Var1)
    numall2$numcells[numall2$unique.seurat.Celltype_level1. == i] <- nrow(seurat_sub2@meta.data)
  }
}

###########
library(data.table)
celltypes <- unique(seurat$Celltype_level1)
eGenesummary <- as.data.frame(table(seurat$Celltype_level1))
allcells <- as.data.frame(table(seurat$Celltype_level0))
eGenesummary <- rbind(eGenesummary,allcells)
colnames(eGenesummary) <- c("celltype","numCells")
eGenesummary$numCells <- as.numeric(eGenesummary$numCells)
eGenesummary$numdonor <- numall2$numdonor_over20[match(eGenesummary$celltype,numall2$unique.seurat.Celltype_level1.)]
eGenesummary$numCells_ineQTL <- numall2$numcells[match(eGenesummary$celltype,numall2$unique.seurat.Celltype_level1.)]
eGenesummary$numdonor[eGenesummary$celltype == "All"] <- length(unique(seurat$donor_id))
eGenesummary$numCells_ineQTL[eGenesummary$celltype == "All"] <- nrow(seurat@meta.data)
eGenesummary <- eGenesummary[eGenesummary$numdonor > 20,]

eGenesummary$numGenes_tested <- NA
eGenesummary$numGenes_sig <- NA
eGenesummary$numGenes_conditionalsig <- NA

##make list for all genes tested - and which overlap between celltypes
genestested <- as.data.frame(matrix(ncol=2))
colnames(genestested) <- c("phenotype_id","celltype")

##make list for all genes significant - and which overlap between celltypes
genessig <- as.data.frame(matrix(ncol=6))
colnames(genessig) <- c("phenotype_id","celltype","variant_id","start_distance","af","slope")

resultsdir <- "/path/onek1k_eQTLresults/"
resultsdir_all <- "/path/onek1k_eQTLresults/"
for (i in 1:nrow(eGenesummary)){
  celltype <- eGenesummary$celltype[[i]]
  optimPCs <- eGenesummary$optimPCs[[i]]
  
  if(celltype == "All"){
    sigegene <- read.table(paste0(resultsdir_all,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv"),header=T,fill=T)
    sigegene_conditional <- read.table(paste0(resultsdir_all,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T,fill=T)
  }else{
    sigegene <- read.table(paste0(resultsdir,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Cis_eqtls_qval.tsv"),header=T,fill=T)
    sigegene_conditional <- read.table(paste0(resultsdir,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T,fill=T)
  }
  
  eGenesummary[i, 'numGenes_tested'] <- nrow(sigegene)
  sigegene$celltype <- celltype
  genestested <- rbind(genestested,sigegene[,c("phenotype_id","celltype")])
  sigegene <- sigegene[sigegene$qval < 0.05,]
  eGenesummary[i, 'numGenes_sig'] <- nrow(sigegene)
  genessig <- rbind(genessig,sigegene[,c("phenotype_id","celltype","variant_id","start_distance","af","slope")])
  eGenesummary[i, 'numGenes_conditionalsig'] <- nrow(sigegene_conditional)
}

write.csv(eGenesummary,paste0("/path/onek1k_eQTLresults/","eGenesummary_Onek1k.csv"))

