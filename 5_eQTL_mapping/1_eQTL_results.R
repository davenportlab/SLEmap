### plot eQTL summary
library(ggplot2)
library(scales)
library(viridis)
library(schard)
library(ggrepel)
library(ggpubr)
library(data.table)

#convert anndata object to seurat
seurat = schard::h5ad2seurat("/path/testrun_input_n281.h5ad")
saveRDS(seurat,"/path/testrun_input_n281_seurat.RDS")

#get number of donors tested for each cell type
numall <- list()
for(i in 1:length(unique(seurat$Celltype_level1))){
  seurat_sub <- subset(seurat,subset = Celltype_level1 == unique(seurat$Celltype_level1)[[i]])
  numcells <- as.data.frame(table(seurat_sub$experiment_id))
  numall[[i]] <- nrow(numcells[numcells$Freq >= 20,])
}
numall2 <- data.frame(unique(seurat$Celltype_level1))
numall2$numdonor_over20 <- NA
for(i in 1:length(unique(seurat$Celltype_level1))){numall2$numdonor_over20[i] <- numall[[i]]}
write.csv(numall2,"/path/celltypes_over20cells_over100donors.csv")


## get number of cells per celltype
celltypes <- unique(seurat$Celltype_level1)
eGenesummary <- as.data.frame(table(seurat$Celltype_level1))
allcells <- as.data.frame(table(seurat$Celltype_level0))
eGenesummary <- rbind(eGenesummary,allcells)
colnames(eGenesummary) <- c("celltype","numCells")
eGenesummary$numCells <- as.numeric(eGenesummary$numCells)

##get optim PCs
resultsdir <- "/path/eQTLresults/"
resultsdir_all <- "/path/eQTLresults/"
PCs <- read.csv("/path/optim_pcs.csv")
eGenesummary$optimPCs <- PCs$PCs_use[match(eGenesummary$celltype,PCs$Celltype)]
eGenesummary$optimPCs[eGenesummary$celltype == "All"] <- 25
eGenesummary <- eGenesummary[is.na(eGenesummary$optimPCs) == F,]

bydonor <- as.data.frame(table(seurat$Celltype_level1,seurat$WGS_ID))
eGenesummary$numCells_inc_eQTL <- NA
for(i in c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","DN_T_cells","Memory_B_cells","Naive_B_cells")){
  bydonor_tmp <- bydonor[bydonor$Var1 == i,]
  donorlist <- read.table(paste0("/path/eQTLresults/",i,"/results/norm_data/dMean__",i,"_all/mappings_handeling_repeats.tsv"),header=T)
  bydonor_tmp <- bydonor_tmp[bydonor_tmp$Var2 %in% donorlist$Genotype,]
  eGenesummary$numCells_inc_eQTL[eGenesummary$celltype == i] <- sum(bydonor_tmp$Freq)
}

bydonor_allcells <- as.data.frame(table(seurat$Celltype_level0,seurat$WGS_ID))
donorlist <- read.table(paste0("/path/eQTLresults/norm_data/dMean__All_all/mappings_handeling_repeats.tsv"),header=T)
bydonor_allcells <- bydonor_allcells[bydonor_allcells$Var2 %in% donorlist$Genotype,]
eGenesummary$numCells_inc_eQTL[eGenesummary$celltype == "All"] <- sum(bydonor_allcells$Freq)

eGenesummary$numGenes_tested <- NA
eGenesummary$numGenes_sig <- NA
eGenesummary$numGenes_conditionalsig <- NA

##make list for all genes tested - and which overlap between celltypes
genestested <- as.data.frame(matrix(ncol=2))
colnames(genestested) <- c("phenotype_id","celltype")

##make list for all genes significant - and which overlap between celltypes
genessig <- as.data.frame(matrix(ncol=6))
colnames(genessig) <- c("phenotype_id","celltype","variant_id","start_distance","af","slope")

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

eGenesummary$cellgroup <- NA
eGenesummary$cellgroup[eGenesummary$celltype %in% c("CD56Bright_NK_cells","CD56Dim_NK_cells")] <- "NK"
eGenesummary$cellgroup[eGenesummary$celltype %in% c("Classical_Monocytes")] <- "Mono"
eGenesummary$cellgroup[eGenesummary$celltype %in% c("CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells")]<- "CD4_T"
eGenesummary$cellgroup[eGenesummary$celltype %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
eGenesummary$cellgroup[eGenesummary$celltype %in% c("DN_T_cells")] <- "Other_T"
eGenesummary$cellgroup[eGenesummary$celltype %in% c("Memory_B_cells","Naive_B_cells")] <- "B"
eGenesummary$cellgroup[eGenesummary$celltype %in% c("All")] <- "All"

eGenesummary$celltype_forplots <- gsub("_", " ", eGenesummary$celltype)
eGenesummary$celltype_forplots <- factor(eGenesummary$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
eGenesummary$cellgroup_forplots <- gsub("_", " ", eGenesummary$cellgroup)
eGenesummary$cellgroup_forplots <- factor(eGenesummary$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))

numdonors <- read.csv("/path/celltypes_over20cells_over100donors.csv")
eGenesummary$ndonors <- numdonors$numdonor_over20[match(eGenesummary$celltype,numdonors$unique.seurat.Celltype_level1.)]
eGenesummary$ndonors[eGenesummary$celltype == "All"] <- 281

##average number of cells per donor
eGenesummary$avg_cell_perdonor <- eGenesummary$numCells_inc_eQTL / eGenesummary$ndonors
ggplot(eGenesummary, aes(x=avg_cell_perdonor, y=prop_sigeGenes)) + geom_point(aes(color=avg_cell_perdonor),size=2)+ scale_x_continuous(labels = label_comma())+theme_classic()+geom_text_repel(aes(label = celltype_forplots),max.overlaps = 3,min.segment.length = 0.1,segment.color = NA,size=3)+xlab("Number of donors")+ylab("Proportion of significant eGenes out of genes tested")+stat_cor(method="spearman")+scale_color_viridis(option="inferno",end=0.85)
res <- cor.test(eGenesummary$avg_cell_perdonor,eGenesummary$prop_sigeGenes, method = "spearman")
res$p.value

##number of donors
ggplot(eGenesummary, aes(x=ndonors, y=prop_sigeGenes)) + geom_point(aes(color=avg_cell_perdonor),size=2)+ scale_x_continuous(labels = label_comma())+theme_classic()+geom_text_repel(aes(label = celltype_forplots),max.overlaps = 3,min.segment.length = 0.1,segment.color = NA,size=3)+xlab("Number of donors")+ylab("Proportion of significant eGenes out of genes tested")+stat_cor(method="spearman")+scale_color_viridis(option="inferno",end=0.85)
ggsave(paste0(resultsdir,"0_plots/ndonor_correlation_perc_colorcell.pdf"),width=7,height=5)

write.csv(eGenesummary,paste0(resultsdir,"1_csvfiles/eGenesummary.csv"),row.names=F)


