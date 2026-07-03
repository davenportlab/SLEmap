####PCA plot for paper with ancestry colored
library(ggplot2)

covariates <- read.table("/path/eQTLresults/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_input/base_output__base/Covariates.tsv",sep="\t",row.names=1,header=T)
covariates <- as.data.frame(t(covariates))

ancestry <- read.csv("/path/KING.csv")
covariates$Ancestry <- ancestry$Ancestry[match(row.names(covariates),ancestry$WGS_ID)]
covariates$Ancestry_tidy <- covariates$Ancestry
table(covariates$Ancestry_tidy)
covariates$Ancestry_tidy[covariates$Ancestry_tidy %in% c("AFR;EUR","EAS;AMR","EAS;EUR;AMR","EUR;AFR","EUR;AMR","EUR;EAS")] <- "Unassigned"
covariates$Ancestry_tidy <- factor(covariates$Ancestry_tidy,levels=c("AFR","EUR","SAS","EAS","AMR","Unassigned"))
table(covariates$Ancestry_tidy)
colorpalette <- c(
  "#4E79A7",
  "#F28E2B",
  "#59A14F",
  "#E15759",
  "#B07AA1",
  "grey"
)

colnames(covariates)
ggplot(covariates, aes(x=`Genotype PC1`, y=`Genotype PC2`,color=Ancestry_tidy)) + geom_point() +theme_classic() +scale_color_manual(values = colorpalette)+ theme(legend.position="top")
ggsave("/path/coloc/1_plots/ancestry_pca.pdf",width=3.5,height=4)
