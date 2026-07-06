## plot genotype PCs by ancestry to decide number of PCs to include in eQTL mapping - using 2PCs for all ancestries

library(data.table)
library(ggplot2)
library(ggpubr)
library(rlist)

plotlist <- list()
for(i in c("AFR","SAS","EUR")){
  pca <- fread(paste0("/path/ancestry_coloc/",i,"/FindGenotypePCs/pca.eigenvec"))
  colnames(pca)[1:2] <- c("FID", "IID")
  eigenval <- scan(paste0("/path/ancestry_coloc/",i,"/FindGenotypePCs/pca.eigenval"))
  var <- round(eigenval / sum(eigenval) * 100, 2)
  
  p1 <- ggplot(pca, aes(PC1, PC2)) +
    geom_point(size = 2) +
    theme_classic() +
    xlab(paste0("PC1 (", var[1], "%)")) +
    ylab(paste0("PC2 (", var[2], "%)"))+ggtitle(i)
  p2 <- ggplot(pca, aes(PC3, PC4)) +
    geom_point(size = 2) +
    theme_classic() +
    xlab(paste0("PC3 (", var[3], "%)")) +
    ylab(paste0("PC4 (", var[4], "%)"))+ggtitle(i)
  p3 <- ggplot(pca, aes(PC5, PC6)) +
    geom_point(size = 2) +
    theme_classic() +
    xlab(paste0("PC5 (", var[5], "%)")) +
    ylab(paste0("PC6 (", var[6], "%)"))+ggtitle(i)
  
  pc_df <- data.frame(
    PC = seq_along(eigenval),
    VarianceExplained = eigenval / sum(eigenval) * 100
  )
  p4 <- ggplot(pc_df, aes(x = PC, y = VarianceExplained)) +
    geom_point(size = 2) +
    geom_line() +
    theme_classic() +
    xlab("Principal Component") +
    ylab("Variance explained (%)")+ggtitle(i)
  plotlist <- list.append(plotlist,ggarrange(p1,p2,p3,p4,ncol=2,nrow=2))
}

plotlist
pdf(file = "/path/ancestry_coloc/genotypePCs_maf005.pdf",width=10,height=4)
for(i in 1:length(plotlist)){
  print(plotlist[[i]])
}
dev.off()
