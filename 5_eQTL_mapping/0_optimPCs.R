### identify optim PCs 
PCs <- c(0,5,10,15,20,25,30,35,40,45,50)
celltypes <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","DN_T_cells","Memory_B_cells","Naive_B_cells")

eQTLresults <- data.frame(matrix(0, ncol = 4))
colnames(eQTLresults) <- c("Celltype","PCs","num_eGenes","optim")
celltypes <- celltypes$unique.seurat.Celltype_level1.

for (i in celltypes){
  optim <- read.table(paste0("/path/eQTLresults/TensorQTL_eQTLS/dMean__",i,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
  optim <- optim[1,1]
  for (k in PCs){
    QTL <- read.table(paste0("/path/eQTLresults/TensorQTL_eQTLS/dMean__",i,"_all/",k,"pcs/base_output/base/Cis_eqtls_qval.tsv"), fill = TRUE, header=T)
    if (k == optim){
      eQTLresults <- rbind(eQTLresults,c(i,k,nrow(QTL[QTL$qval < 0.05,]),"optim"))
    }else{
      eQTLresults <- rbind(eQTLresults,c(i,k,nrow(QTL[QTL$qval < 0.05,]),"notoptim"))
    }
  }
}

eQTLresults <- eQTLresults[-1,]
eQTLresults$num_eGenes <- as.numeric(eQTLresults$num_eGenes)
eQTLresults$PCs <- as.numeric(eQTLresults$PCs)

eQTLresults_diff <- eQTLresults %>% group_by(Celltype) %>% mutate(Diff = num_eGenes - lag(num_eGenes))
eQTLresults_diff$V5 <- NA
eQTLresults_diff$V5[eQTLresults_diff$PCs == 5] <- "0-5"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 10] <- "5-10"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 15] <- "10-15"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 20] <- "15-20"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 25] <- "20-25"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 30] <- "25-30"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 35] <- "30-35"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 40] <- "35-40"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 45] <- "40-45"
eQTLresults_diff$V5[eQTLresults_diff$PCs == 50] <- "45-50"
eQTLresults_diff$V5 <- factor(eQTLresults_diff$V5, levels=c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50"))

library(ggrepel)
for (i in celltypes){
  p1 <- ggplot(eQTLresults_diff[eQTLresults_diff$Celltype == i,], aes(x=PCs, y=num_eGenes, group=1)) +
    geom_line()+theme_classic()+ggtitle(paste0("Number of eGenes in ",i))+ylab("Number of eGenes")+xlab("Expression PCs")+ expand_limits(y = 0)
  p2 <- ggplot(eQTLresults_diff[(eQTLresults_diff$Celltype == i)&(is.na(eQTLresults_diff$Diff) == F),], aes(x=V5, y=Diff, group=1)) + geom_line()+theme_classic()+ggtitle(paste0("Difference in number of eGenes in ",i))+ylab("Change in number of eGenes")+xlab("Expression PCs")+ geom_hline(yintercept=0, linetype="dashed", color = "red")
  ggarrange(p1,p2)
  ggsave(paste0("/path/Results/3_eQTL_prep/ExpressionPC_optimization_allSNP/ExpressionPC_optimisation_",i,".pdf"),width=14,height=6)
}

##identify optim PCs in all-pbmc
PCs <- c(0,5,10,15,20,25,30,35,40,45,50)

eQTLresults <- data.frame(matrix(0, ncol = 2))
colnames(eQTLresults) <- c("PCs","num_eGenes")

for (k in PCs){
  QTL <- read.table(paste0("/eQTLresults/TensorQTL_eQTLS/dMean__All_all/",k,"pcs/base_output/base/Cis_eqtls_qval.tsv"), fill = TRUE, header=T)
  eQTLresults <- rbind(eQTLresults,c(k,nrow(QTL[QTL$qval < 0.05,])))
}

eQTLresults <- eQTLresults[-1,]
eQTLresults$num_eGenes <- as.numeric(eQTLresults$num_eGenes)
eQTLresults$PCs <- as.numeric(eQTLresults$PCs)

ggplot(eQTLresults, aes(x=PCs, y=num_eGenes, group=1)) +
  geom_line()+theme_classic()+ggtitle(paste0("Number of eGenes in ","all cells pseudobulked"))+ylab("Number of eGenes")+xlab("Expression PCs")+ expand_limits(y = 0)
