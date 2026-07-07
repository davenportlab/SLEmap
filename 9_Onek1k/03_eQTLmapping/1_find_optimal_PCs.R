
### 1. SET UP ####
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)

### 2. identify optim PCs - per cell type ####
PCs <- c(0,5,10,15,20,25,30,35,40,45,50)
celltypes <- read.csv("/path/oneK1K_celltype_passfilter.txt", header = FALSE) #more than 100 individuals with more than 20 cells per cell type
celltypes = celltypes$V1

#compile the results - make a df: celltype, PC, number of egenes, optim/notoptim
eQTLresults <- data.frame(matrix(0, ncol = 3))
colnames(eQTLresults) <- c("Celltype","PCs","num_eGenes")

for (i in celltypes){
  for (k in PCs){
    QTL <- read.table(paste0("/path/onek1k_eQTLresults/results/TensorQTL_eQTLS/dMean__",i,"_all/",k,"pcs/base_output/base/Cis_eqtls_qval.tsv"), fill = TRUE, header=T)
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
pdf("/path/onek1k_eQTLresults/ExpressionPC_optimisation_ALL.pdf",
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
}
dev.off()

#make df of the selected PCs
Celltype_nPC <- data.frame(
  Cell_type1 = celltypes,
  nPC = c(10,15,20,5,10,10,10,15,15,20,30,5,5,5,10),
  stringsAsFactors = FALSE
)
write.csv(Celltype_nPC, "/path/onek1k_eQTLresults/nPC_selected_percelltype.csv", row.names = FALSE)

