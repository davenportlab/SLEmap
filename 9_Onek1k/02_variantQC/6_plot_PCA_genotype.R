### 1. SET UP ####
library(data.table)
library(ggplot2)
library(ggpubr)

### 2. make a pca from the genotyping data #####
#PCA information is calcualted from plink2 command

pca <- fread("/path/oneK1K_imputation/Imputation_TOPMed_run3/find_genotype_pc/onek1k_allchr_pca.eigenvec")
colnames(pca)[1:2] <- c("FID", "IID")
eigenval <- scan("/path/oneK1K_imputation/Imputation_TOPMed_run3/find_genotype_pc/onek1k_allchr_pca.eigenval")
var <- round(eigenval / sum(eigenval) * 100, 2)
p1 <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(size = 1) +
  theme_classic() +
  xlab(paste0("PC1 (", var[1], "%)")) +
  ylab(paste0("PC2 (", var[2], "%)"))
p2 <- ggplot(pca, aes(PC3, PC4)) +
  geom_point(size = 1) +
  theme_classic() +
  xlab(paste0("PC3 (", var[3], "%)")) +
  ylab(paste0("PC4 (", var[4], "%)"))
p3 <- ggplot(pca, aes(PC5, PC6)) +
  geom_point(size = 1) +
  theme_classic() +
  xlab(paste0("PC5 (", var[5], "%)")) +
  ylab(paste0("PC6 (", var[6], "%)"))

### 2. select the right number of dim pca ####
#saved var as df > pc_df
pc_df <- data.frame(
  PC = seq_along(eigenval),
  VarianceExplained = eigenval / sum(eigenval) * 100
)

#creating elbox plot 
p4 <- ggplot(pc_df, aes(x = PC, y = VarianceExplained)) +
  geom_point(size = 2) +
  geom_line() +
  theme_classic() +
  xlab("Principal Component") +
  ylab("Variance explained (%)")

plotlist <- list(ggarrange(p1, p2, p3, p4, ncol=4))
ggarrange(plotlist = plotlist, ncol = 1)