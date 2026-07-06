##eQTL by ancestry - identify optim PCs

eQTLresults <- data.frame(matrix(0, ncol = 5))
colnames(eQTLresults) <- c("Ancestry","Celltype","PCs","num_eGenes","optim")
celltypes <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","DN_T_cells","Memory_B_cells","Naive_B_cells")
ancestry <- c("AFR","EUR","SAS")
PCs <- c(0,5,10,15,20,25,30,35,40,45,50)

for (j in ancestry){
  basedir <- paste0("/path/ancestry_coloc/",j,"/FindPCs/results")
  for (i in celltypes){
    optim <- read.table(paste0(basedir,"/TensorQTL_eQTLS/dMean__",i,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
    optim <- optim[1,1]
    dirs <- list.dirs(paste0(basedir,"/TensorQTL_eQTLS/dMean__",i,"_all"), full.names = TRUE)
    numPC <- sum(grepl("5pc", basename(dirs))) + sum(grepl("0pc", basename(dirs)))
    for (k in PCs[1:numPC]){
      QTL <- read.table(paste0(basedir,"/TensorQTL_eQTLS/dMean__",i,"_all/",k,"pcs/base_output/base/Cis_eqtls_qval.tsv"), fill = TRUE, header=T)
      if (k == optim){
        eQTLresults <- rbind(eQTLresults,c(j,i,k,nrow(QTL[QTL$qval < 0.05,]),"optim"))
      }else{
        eQTLresults <- rbind(eQTLresults,c(j,i,k,nrow(QTL[QTL$qval < 0.05,]),"notoptim"))
      }
    }
  }
}


eQTLresults <- eQTLresults[-1,]
eQTLresults$num_eGenes <- as.numeric(eQTLresults$num_eGenes)
eQTLresults$PCs <- as.numeric(eQTLresults$PCs)

eQTLresults_diff <- eQTLresults %>% group_by(Ancestry,Celltype) %>% mutate(Diff = num_eGenes - lag(num_eGenes))
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
library(ggpubr)
library(rlist)
plotlist <- list()
for (j in ancestry){
  for (i in celltypes){
    p1 <- ggplot(eQTLresults_diff[(eQTLresults_diff$Ancestry == j)&(eQTLresults_diff$Celltype == i),], aes(x=PCs, y=num_eGenes, group=1)) +geom_line()+theme_classic()+ggtitle(paste0(j,"-Num eGenes in ",i))+ylab("Number of eGenes")+xlab("Expression PCs")+ expand_limits(y = 0)+ theme(panel.grid.major.x = element_line(color = "red",size = 0.2,linetype = 2))
    p2 <- ggplot(eQTLresults_diff[(eQTLresults_diff$Ancestry == j)&(eQTLresults_diff$Celltype == i)&(is.na(eQTLresults_diff$Diff) == F),], aes(x=V5, y=Diff, group=1)) + geom_line()+theme_classic()+ggtitle(paste0(j,"-Diff in num eGenes in ",i))+ylab("Change in number of eGenes")+xlab("Expression PCs")+ geom_hline(yintercept=0, linetype="dashed", color = "red")
    plotlist <- list.append(plotlist,ggarrange(p1,p2))
  }
}

pdf(file = "/path/ancestry_coloc/expPCs.pdf",width=10,height=4)
for(i in 1:length(plotlist)){
  print(plotlist[[i]])
}
dev.off()



###for all cells 
eQTLresults <- data.frame(matrix(0, ncol = 5))
colnames(eQTLresults) <- c("Ancestry","Celltype","PCs","num_eGenes","optim")
celltypes <- "All"
ancestry <- c("AFR","EUR","SAS")
PCs <- c(0,5,10,15,20,25,30,35,40,45,50)

for (j in ancestry){
  basedir <- paste0("/path/ancestry_coloc/",j,"/FindPCs_all/results")
  for (i in celltypes){
    optim <- read.table(paste0(basedir,"/TensorQTL_eQTLS/dMean__",i,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
    optim <- optim[1,1]
    dirs <- list.dirs(paste0(basedir,"/TensorQTL_eQTLS/dMean__",i,"_all"), full.names = TRUE)
    numPC <- sum(grepl("5pc", basename(dirs))) + sum(grepl("0pc", basename(dirs)))
    for (k in PCs[1:numPC]){
      QTL <- read.table(paste0(basedir,"/TensorQTL_eQTLS/dMean__",i,"_all/",k,"pcs/base_output/base/Cis_eqtls_qval.tsv"), fill = TRUE, header=T)
      if (k == optim){
        eQTLresults <- rbind(eQTLresults,c(j,i,k,nrow(QTL[QTL$qval < 0.05,]),"optim"))
      }else{
        eQTLresults <- rbind(eQTLresults,c(j,i,k,nrow(QTL[QTL$qval < 0.05,]),"notoptim"))
      }
    }
  }
}


eQTLresults <- eQTLresults[-1,]
eQTLresults$num_eGenes <- as.numeric(eQTLresults$num_eGenes)
eQTLresults$PCs <- as.numeric(eQTLresults$PCs)

eQTLresults_diff <- eQTLresults %>% group_by(Ancestry,Celltype) %>% mutate(Diff = num_eGenes - lag(num_eGenes))
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
library(ggpubr)
library(rlist)
plotlist <- list()
for (j in ancestry){
  for (i in celltypes){
    p1 <- ggplot(eQTLresults_diff[(eQTLresults_diff$Ancestry == j)&(eQTLresults_diff$Celltype == i),], aes(x=PCs, y=num_eGenes, group=1)) +geom_line()+theme_classic()+ggtitle(paste0(j,"-Num eGenes in ",i))+ylab("Number of eGenes")+xlab("Expression PCs")+ expand_limits(y = 0)+ theme(panel.grid.major.x = element_line(color = "red",size = 0.2,linetype = 2))
    p2 <- ggplot(eQTLresults_diff[(eQTLresults_diff$Ancestry == j)&(eQTLresults_diff$Celltype == i)&(is.na(eQTLresults_diff$Diff) == F),], aes(x=V5, y=Diff, group=1)) + geom_line()+theme_classic()+ggtitle(paste0(j,"-Diff in num eGenes in ",i))+ylab("Change in number of eGenes")+xlab("Expression PCs")+ geom_hline(yintercept=0, linetype="dashed", color = "red")
    plotlist <- list.append(plotlist,ggarrange(p1,p2))
  }
}
plotlist


### set optimPCs
optimPCs_ancestry <- expand.grid(
  celltypes = celltypes,
  ancestry = ancestry,
  stringsAsFactors = FALSE
)
optimPCs_ancestry$optimPCs <- NA

optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="Naive_B_cells"] <- 20
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="CD56Bright_NK_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="CD56Dim_NK_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="Naive_CD4_T_cells"] <- 20
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="CM_CD4_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="Naive_CD8_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="TEMRA"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="EM_CD4_T_cells"] <- 5 #
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="EM_CD8_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="Regulatory_CD4_T_cells"] <- 15 #
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="CM_CD8_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="Cytotoxic_CD4_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="DN_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="Memory_B_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="Classical_Monocytes"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="Naive_B_cells"] <- 20
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="CD56Bright_NK_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="CD56Dim_NK_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="Naive_CD4_T_cells"] <- 15
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="CM_CD4_T_cells"] <- 20
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="Naive_CD8_T_cells"] <- 15
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="TEMRA"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="EM_CD4_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="EM_CD8_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="Regulatory_CD4_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="CM_CD8_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="Cytotoxic_CD4_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="DN_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="Memory_B_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="Classical_Monocytes"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="Naive_B_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="CD56Bright_NK_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="CD56Dim_NK_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="Naive_CD4_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="CM_CD4_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="Naive_CD8_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="TEMRA"] <- 15
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="EM_CD4_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="EM_CD8_T_cells"] <- 10
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="Regulatory_CD4_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="CM_CD8_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="Cytotoxic_CD4_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="DN_T_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="Memory_B_cells"] <- 5
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="Classical_Monocytes"] <- 5

optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="AFR" & optimPCs_ancestry$celltypes=="All"] <- 20
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="SAS" & optimPCs_ancestry$celltypes=="All"] <- 15
optimPCs_ancestry$optimPCs[optimPCs_ancestry$ancestry=="EUR" & optimPCs_ancestry$celltypes=="All"] <- 20

write.csv(optimPCs_ancestry,"/path/ancestry_coloc/optimPCs.csv",row.names=F)