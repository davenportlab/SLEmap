
###eQTLs are run for onek1k in /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping

### 1. SET UP ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)

### 2. identify optim PCs - per cell type ####
PCs <- c(0,5,10,15,20,25,30,35,40,45,50)
celltypes <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/eb30/results/oneK1K_celltype_passfilter.txt", header = FALSE)
celltypes = celltypes$V1

#compile the results - make a df: celltype, PC, number of egenes, optim/notoptim
eQTLresults <- data.frame(matrix(0, ncol = 3))
colnames(eQTLresults) <- c("Celltype","PCs","num_eGenes")

for (i in celltypes){
  for (k in PCs){
    QTL <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_findPCs/results/TensorQTL_eQTLS/dMean__",i,"_all/",k,"pcs/base_output/base/Cis_eqtls_qval.tsv"), fill = TRUE, header=T)
    eQTLresults <- rbind(eQTLresults,c(i,k,nrow(QTL[QTL$qval < 0.05,])))
  }
}
eQTLresults <- eQTLresults[-1,]
eQTLresults$num_eGenes <- as.numeric(eQTLresults$num_eGenes)
eQTLresults$PCs <- as.numeric(eQTLresults$PCs)

#create seperate column to show the different number of genes to the previous PCs
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

#plot the PCs
pdf("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/ExpressionPC_optimisation_ALL.pdf",
    width = 14, height = 6)
for (i in celltypes){
  p1 <- ggplot(eQTLresults_diff[eQTLresults_diff$Celltype == i, ], aes(x=PCs, y=num_eGenes, group=1)) + 
    geom_line() + theme_classic() + 
    geom_point(size = 2) +
    ggtitle(paste0("Number of eGenes in ",i)) + 
    ylab("Number of eGenes") + xlab("Expression PCs") + 
    expand_limits(y = 0)
  
  p2 <- ggplot(eQTLresults_diff[(eQTLresults_diff$Celltype == i)&(is.na(eQTLresults_diff$Diff) == F),], aes(x=V5, y=Diff, group=1)) + geom_line()+theme_classic()+ggtitle(paste0("Difference in number of eGenes in ",i))+ylab("Change in number of eGenes")+xlab("Expression PCs")+ geom_hline(yintercept=0, linetype="dashed", color = "red")
  
  print(ggarrange(p1, p2, ncol = 2))
  #ggarrange(p1,p2)
  #ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/eb30/results/eQTL_mapping_findPCs/plot_evaluating_PCs/ExpressionPC_optimisation_",i,".pdf"),width=14,height=6)
}
dev.off()

#make df of the selected PCs
Celltype_nPC <- data.frame(
  Cell_type1 = celltypes,
  nPC = c(10,15,20,5,10,10,10,15,15,20,30,5,5,5,10),
  stringsAsFactors = FALSE
)
write.csv(Celltype_nPC, "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/nPC_selected_percelltype.csv", row.names = FALSE)

###compare with Eli's 
Eli_nPC <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/eb30/results/eQTL_mapping_findPCs/plot_evaluating_PCs/nPC_selected_percelltype.csv")
Celltype_nPC$Eli <- Eli_nPC$nPC[match(Celltype_nPC$Cell_type1,Eli_nPC$Cell_type1)]

#make directory of each cell type
new_dir <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs/"
setwd(new_dir)

for (i in celltypes) {
  dir.create(i)
}

### 3. identify optim PCs - for all cells ####
PCs <- c(0,5,10,15,20,25,30,35,40,45,50)

#compile the results - make a df: celltype, PC, number of egenes, optim/notoptim
eQTLresults <- data.frame(matrix(0, ncol = 3))
colnames(eQTLresults) <- c("PCs","num_eGenes","optim")
optim <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_findPCs_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))[1,1]

for (k in PCs){
  QTL <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_findPCs_allcells/results/TensorQTL_eQTLS/dMean__All_all/",k,"pcs/base_output/base/Cis_eqtls_qval.tsv"), fill = TRUE, header=T)
    if (k == optim){
      eQTLresults <- rbind(eQTLresults,c(k,nrow(QTL[QTL$qval < 0.05,]),"optim"))
    }else{
      eQTLresults <- rbind(eQTLresults,c(k,nrow(QTL[QTL$qval < 0.05,]),"notoptim"))
    }
}
eQTLresults <- eQTLresults[-1,]
eQTLresults$num_eGenes <- as.numeric(eQTLresults$num_eGenes)
eQTLresults$PCs <- as.numeric(eQTLresults$PCs)

#create seperate column to show the different number of genes to the previous PCs
eQTLresults_diff <- eQTLresults %>% mutate(Diff = num_eGenes - lag(num_eGenes))

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

#plot the PCs
p1 <- ggplot(eQTLresults_diff, aes(x=PCs, y=num_eGenes, group=1)) + 
    geom_line() + theme_classic() + 
    geom_point(aes(color = optim), size = 2) +
    scale_color_manual(values = c("notoptim" = "black", "optim" = "red")) +
    ggtitle("Number of eGenes in All") + 
    ylab("Number of eGenes") + xlab("Expression PCs") + 
    expand_limits(y = 0)
  
p2 <- ggplot(eQTLresults_diff[(is.na(eQTLresults_diff$Diff) == F),], aes(x=V5, y=Diff, group=1)) + 
  geom_line() + theme_classic() + 
  ggtitle("Difference in number of eGenes in All") + 
  ylab("Change in number of eGenes") + xlab("Expression PCs") + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")

ggarrange(p1,p2)
ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_findPCs_allcells/ExpressionPC_optimisation.pdf"),width=14,height=6)


nPC <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/nPC_selected_percelltype.csv")
nPC <- rbind(nPC,c("All",20))
write.csv(nPC,"/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/nPC_selected_percelltype.csv", row.names = FALSE)

# 
# 
# ## 
# # eQTLresults_diff_neg <- eQTLresults_diff[eQTLresults_diff$Diff < 0,]
# # eQTLresults_diff_neg  <- eQTLresults_diff_neg [!duplicated(eQTLresults_diff_neg $Celltype),]
# # eQTLresults_diff_neg$PCs_use <- eQTLresults_diff_neg$PCs - 5
# # eQTLresults_diff_neg$use <- "use"
# # eQTLresults$use <- eQTLresults_diff_neg$use[match(paste0(eQTLresults$Celltype,eQTLresults$PCs),paste0(eQTLresults_diff_neg$Celltype,eQTLresults_diff_neg$PCs_use))]
# # eQTLresults$usePC <- NA
# # eQTLresults$usePC[is.na(eQTLresults$use) == F] <- eQTLresults$PCs[is.na(eQTLresults$use) == F]
# 
# ##Using number of PCs manually chosen. 
# optim_PCs <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/3_eQTL_prep/ExpressionPC_optimization_allSNP/optim_pcs.csv")
# 
# eQTLresults$manualoptim <- NA
# eQTLresults$manualoptim[paste0(eQTLresults$Celltype,eQTLresults$PCs) %in% paste0(optim_PCs$Celltype,optim_PCs$PCs_use)] <- eQTLresults$PCs[paste0(eQTLresults$Celltype,eQTLresults$PCs) %in% paste0(optim_PCs$Celltype,optim_PCs$PCs_use)]
# 
# library(ggrepel)
# for (i in celltypes){
#   p1 <- ggplot(eQTLresults[eQTLresults$Celltype == i,], aes(x=PCs, y=num_eGenes, group=1)) +
#     geom_line()+theme_classic()+ggtitle(paste0("Number of eGenes in ",i))+ylab("Number of eGenes")+xlab("Expression PCs")+ expand_limits(y = 0)+geom_point(data = eQTLresults[(eQTLresults$Celltype == i)&(is.na(eQTLresults$manualoptim)==F),], aes(x=PCs, y=num_eGenes), size = 2,color="red")+geom_text_repel(data = eQTLresults[(eQTLresults$Celltype == i)&(is.na(eQTLresults$manualoptim)==F),], aes(x=PCs, y=num_eGenes,label=manualoptim),nudge_y=10)
#   p2 <- ggplot(eQTLresults_diff[(eQTLresults_diff$Celltype == i)&(is.na(eQTLresults_diff$Diff) == F),], aes(x=V5, y=Diff, group=1)) + geom_line()+theme_classic()+ggtitle(paste0("Difference in number of eGenes in ",i))+ylab("Number of eGenes")+xlab("Expression PCs")+ geom_hline(yintercept=0, linetype="dashed", color = "red")
#   ggarrange(p1,p2)
#   ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/3_eQTL_prep/ExpressionPC_optimization_allSNP/ExpressionPC_optimisation_",i,".pdf"),width=14,height=6)
# }
# 
# eQTLresults$Celltype_forplots <- gsub("_", " ", eQTLresults$Celltype)
# eQTLresults$Celltype_forplots <- factor(eQTLresults$Celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells"))
# ggplot(eQTLresults,aes(x=PCs, y=num_eGenes, group=Celltype_forplots,color=Celltype_forplots)) +
#   geom_line()+geom_point(data = eQTLresults[(is.na(eQTLresults$manualoptim)==F),], aes(x=PCs, y=num_eGenes), size = 1.5,color="black")+theme_classic()+ggtitle("Number of expression PCs used as covariates in the eQTL model")+ylab("Number of eGenes")+xlab("Expression PCs")
# ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/3_eQTL_prep/ExpressionPC_optimization_allSNP/ExpressionPC_optimisation_","allcelltypes",".pdf"),width=8,height=6)
# 
# 
