####collectri

library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(stringr)
library(rlist)
library(ggpubr)

##make result dir
resultdir <- "/path/motifbreakr/"
dir.create(file.path(paste0(resultdir,"decoupler")))

##get motifbreakr results - TRAF1 
motifbreakr <- read.csv(paste0(resultdir,"strong_results/TRAF1_CM_CD4_T_cells_Memory_B_cells.csv"))
cell_type <- "CM_CD4_T_cells"

##get TF info and subset to tfs of interest
net <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)
net_use <- net[net$source %in% motifbreakr$geneSymbol,]
table(net_use$source)
rm(net)

##using pseudobulk
SLEmap <- readRDS("/path/eQTL_input_n281_seurat.RDS")
SLEmap_pb <- AverageExpression(SLEmap,group.by = c("Celltype_level1", "WGS_ID"))
gene_symbol <- read.csv("/path/ensemblID_to_genesymbol.csv")
SLEmap_pb[["RNA"]]@Dimnames[[1]] <- gene_symbol$gene_symbols[match(SLEmap_pb[["RNA"]]@Dimnames[[1]],gene_symbol$X)]
saveRDS(SLEmap_pb,paste0(resultdir,"decoupler/SLEmap_pb.RDS"))
SLEmap_pb <- readRDS(paste0(resultdir,"decoupler/SLEmap_pb.RDS"))

##pb for all cells
SLEmap_pb <- AverageExpression(SLEmap,group.by = c("Celltype_level0", "WGS_ID"))
gene_symbol <- read.csv("/path/ensemblID_to_genesymbol.csv")
SLEmap_pb[["RNA"]]@Dimnames[[1]] <- gene_symbol$gene_symbols[match(SLEmap_pb[["RNA"]]@Dimnames[[1]],gene_symbol$X)]
saveRDS(SLEmap_pb,paste0(resultdir,"decoupler/SLEmap_pb_allcells.RDS"))
SLEmap_pb_allcells <- readRDS(paste0(resultdir,"decoupler/SLEmap_pb_allcells.RDS"))
SLEmap_pb_allcells[["RNA"]]@Dimnames[[2]] <- paste0("All_",SLEmap_pb_allcells[["RNA"]]@Dimnames[[2]])

SLEmap_pb <- cbind(SLEmap_pb[["RNA"]],SLEmap_pb_allcells[["RNA"]])

#check expression of TF
SLEmap_pb_TFs <- as.data.frame(t(SLEmap_pb[row.names(SLEmap_pb) %in% c(motifbreakr$geneSymbol,"TRAF1"),]))
SLEmap_pb_TFs$celltype <- sub("_.*", "", rownames(SLEmap_pb_TFs))
p1 <- ggplot(SLEmap_pb_TFs, aes(x=celltype, y=TFAP2C)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(SLEmap_pb_TFs, aes(x=celltype, y=REST)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- ggplot(SLEmap_pb_TFs, aes(x=celltype, y=TFAP2A)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 <- ggplot(SLEmap_pb_TFs, aes(x=celltype, y=PLAG1)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggarrange(p1,p2,p3,p4)

p1 <- ggplot(SLEmap_pb_TFs[SLEmap_pb_TFs$celltype %in% c("Memory-B-cells","CM-CD4-T-cells"),], aes(x=celltype, y=TFAP2C)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(SLEmap_pb_TFs[SLEmap_pb_TFs$celltype %in% c("Memory-B-cells","CM-CD4-T-cells"),], aes(x=celltype, y=REST)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- ggplot(SLEmap_pb_TFs[SLEmap_pb_TFs$celltype %in% c("Memory-B-cells","CM-CD4-T-cells"),], aes(x=celltype, y=TFAP2A)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 <- ggplot(SLEmap_pb_TFs[SLEmap_pb_TFs$celltype %in% c("Memory-B-cells","CM-CD4-T-cells"),], aes(x=celltype, y=PLAG1)) + geom_boxplot() + theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)

##run decouplr ulm
results_pb <- decoupleR::run_ulm(
  mat = SLEmap_pb,
  net = net_use,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

##explore results
results_pb$cell_type <- sub("_.*", "", results_pb$condition)
results_pb$cell_type <- gsub("-","_",results_pb$cell_type)

ggplot(results_pb[results_pb$cell_type %in% c("CM_CD4_T_cells","Memory_B_cells") & results_pb$source %in% c("PLAG1"),], aes(x=cell_type, y=score)) + geom_violin(fill="#C7DFEF") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_boxplot(width=0.1,outliers=F,color="black",fill="white")+xlab("")+ylab("PLAG1 TF activity")+ylim(0,1)
ggsave(paste0(resultdir,"decoupler/TRAF1_plot_forpaper.pdf"),width=3,height=5)

##t-test between cm cd4 t cells and memroy b cells
t.test(results_pb$score[results_pb$cell_type %in% c("CM_CD4_T_cells") & results_pb$source %in% c("PLAG1")], results_pb$score[results_pb$cell_type %in% c("Memory_B_cells") & results_pb$source %in% c("PLAG1")], var.equal = TRUE)