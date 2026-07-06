

alltested <- data.frame()
for(celltype in c("CD56Bright_NK_cells", "CD56Dim_NK_cells", "Classical_Monocytes", "CM_CD4_T_cells", "EM_CD4_T_cells", "Naive_CD4_T_cells", "Regulatory_CD4_T_cells", "CM_CD8_T_cells", "EM_CD8_T_cells", "Naive_CD8_T_cells", "TEMRA", "Memory_B_cells", "Naive_B_cells", "MAIT_and_GammaDelta_T_cells", "Nonclassical_Monocytes")){
  qval <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"),fill=T,header=T)
  qval$celltype <- celltype
  alltested <- rbind(alltested,qval)
  
}

qval <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv",fill=T,header=T)
qval$celltype <- "All"
alltested <- rbind(alltested,qval)

write.csv(alltested,"/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/alltested_genes_onek1k.csv",row.names=F)

write.csv(alltested[alltested$qval < 0.05,],"/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/allsig_genes_onek1k.csv",row.names=F)

###counting numbers of genes tested, egenes, indep eQTL 
alltested <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/alltested_genes_onek1k.csv")
length(unique(alltested$phenotype_id))

length(unique(alltested[alltested$qval < 0.05,]$phenotype_id))

allindependent <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Cis_eqtls_independent_allcelltypes_wbulklike.csv")

length(unique(paste0(allindependent$phenotype_id,allindependent$celltype)))



###slemap
slemap_tested <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/genestested_wbulklike.csv")
length(unique(slemap_tested$phenotype_id))

slemap_sig <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/genessig_wbulklike.csv")
length(unique(slemap_sig$phenotype_id))

slemap_indep <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/conditionaleQTL.csv")
table(slemap_indep$celltype)
table(slemap_indep$rank_new)

slemap_indep_allcells <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",sep="\t",header=T)
table(slemap_indep_allcells$rank)
head(slemap_indep)


### make summary file
### plot eQTL summary
library(ggplot2)
library(scales)

#get number of donors tested for each cell type
seurat <- readRDS("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/onek1K_eQTL_input_seurat.RDS")

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

resultsdir <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs/"
resultsdir_all <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells"
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

write.csv(eGenesummary,paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/","eGenesummary_Onek1k.csv"))

