### identify examples of shared/unshared signals and plot with adjusted expression values for paper fig2- must use module with vcftools

##########plot adjusted expression ############
#library(biomaRt)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#ensembl <- useEnsembl("snp",dataset = "hsapiens_snp")

## set working directory somewhere with write permission - this will be where the genotype files are made
setwd("/path/SNPs")


getadjexp <- function(celltype,Gene,SNP,color="purple",minmax=NA){
  
  eQTLrun <- paste0("/path/eQTLresults/",celltype)
  
  
  ## set this to where the QTLite run results are. There should be a norm_data directory with each cell type
  norm_data_dir <- paste0(eQTLrun,"/results/norm_data/dMean__",celltype,"_all/")
  
  ## get mean expression
  mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
  
  # get SNP position
  SNPs <- setNames(data.frame(unlist(strsplit(SNP, "_"))[1], unlist(strsplit(SNP, "_"))[2], unlist(strsplit(SNP, "_"))[2]), c('chr_name', 'chrom_start','chrom_end'))
  
  ## get genotype
  system2("bash", args = c("/path/genotype.sh", SNPs$chr_name, SNPs$chrom_start,SNPs$chrom_end,SNP))
  indiv <- read.table(paste0(SNP,".012.indv"),sep="\t")
  genotype <- read.table(paste0(SNP,".012.gz"),sep="\t")
  system2("rm",args= c(paste0(SNP,"*")))
  indiv$geno <- genotype$V2
  
  ## plot (mean) gene exp
  mean_use <- as.data.frame(t(mean[Gene,]))
  mean_use$sampleID <- paste0("SLE_",row.names(mean_use) %>% strsplit( ".",fixed=T) %>%  sapply( "[", 3 ) %>% strsplit( "_",fixed=T) %>%  sapply( "[", 3 ))
  mean_use$geno <- indiv$geno[match(mean_use$sampleID, indiv$V1)]
  colnames(mean_use) <- c("Mean_exp","sampleID","Genotype")
  
  ##input covariates
  covariates <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
  
  sampleIDs <- as.character(covariates[1,])
  row.names(covariates) <- paste0(covariates$V1,covariates$V2)
  covariates <- covariates[-1,-c(1:2)]    
  colnames(covariates) <- sampleIDs[1:ncol(covariates)]
  covariates <- mutate_all(as.data.frame(t(covariates)), function(x) as.numeric(as.character(x)))
  
  mean_use <- merge(mean_use,covariates,by.x="sampleID",by.y=0)
  
  # Get optim phenotype PCs and 5 genotypePCs
  optimPCs <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
  optimPCs <- as.numeric(optimPCs[1,1])
  
  genotypePCs <- paste0("GenotypePC",1:5)
  phenotypePCs <- paste0("PhenotypePC",1:optimPCs)
  
  ##tensorQTL imputes genotypes for missing values: they use the mean genotype (not integer)
  mean_use$Genotype_imputed <- mean_use$Genotype
  mean_use$Genotype_imputed[mean_use$Genotype_imputed == -1] <- mean(mean_use$Genotype_imputed[mean_use$Genotype_imputed != -1])
  
  # get full model: the coefficient for genotypes should match the slope from tensorQTL
  null_model <- lm(mean_use[c("Mean_exp", "Genotype_imputed", genotypePCs, phenotypePCs)]) 
  
  coefs <- null_model[["coefficients"]][append(genotypePCs,phenotypePCs)]
  
  mean_use$adjusted_expr <- mean_use$Mean_exp - rowSums(sweep(mean_use[append(genotypePCs,phenotypePCs)], 2, coefs, FUN = "*"))
  beta <- as.numeric(null_model[["coefficients"]][["Genotype_imputed"]])
  
  if(color=="purple"){
    colorpalette <- c("#756BB1","#756BB1","#756BB1")
  }else if(color=="blue"){
    colorpalette <- c("#5C6BC0","#5C6BC0","#5C6BC0")
  }
  
  if(all(is.na(minmax))==T){
      # Plot adjusted expression by genotype
      p2 <- ggplot(mean_use[mean_use$Genotype != -1,], aes(y = adjusted_expr, x = Genotype, group=Genotype,color=factor(Genotype))) + scale_x_continuous(breaks = c(0, 1, 2))+scale_color_manual(values=colorpalette) + geom_jitter(size=1,width=0.25,alpha=0.7)+ geom_boxplot(outlier.shape = NA,alpha=0)+ggtitle(paste0(celltype,Gene,"-",SNP,", n=",nrow(mean_use[mean_use$Genotype != -1,])))+theme_classic()+annotate("text", x = 2, y = min(mean_use[mean_use$Genotype != -1,]$adjusted_expr), label = paste0("b=",round(beta,3)))+ylab("")+xlab("")+ theme(legend.position = "none")

  }else{
    min <- minmax[[1]]
    max <- minmax[[2]]
    p2 <- ggplot(mean_use[mean_use$Genotype != -1,], aes(y = adjusted_expr, x = Genotype, group=Genotype,color=factor(Genotype))) + scale_x_continuous(breaks = c(0, 1, 2))+scale_color_manual(values=colorpalette) + geom_jitter(size=1,width=0.25,alpha=0.7)+ geom_boxplot(outlier.shape = NA,alpha=0)+ggtitle(paste0(celltype,Gene,"-",SNP,", n=",nrow(mean_use[mean_use$Genotype != -1,])))+theme_classic()+annotate("text", x = 2, y = min, label = paste0("b=",round(beta,3)))+ylab("")+xlab("")+ theme(legend.position = "none")+ylim(min,max)
    
  }

  
  return(list(mean_use,p2))
  
}

###ENSG00000154589 - chr8_73989162_G (ref)_A (alt)
LY96_Classical_Monocytes <- getadjexp("Classical_Monocytes","ENSG00000154589","chr8_73989162_G_A","purple",minmax=c(0,0.6))
LY96_Naive_CD8_T_cells <- getadjexp("Naive_CD8_T_cells","ENSG00000154589","chr8_73989162_G_A","purple")
LY96_Allcells <- getadjexp("All","ENSG00000154589","chr8_73989162_G_A","purple",minmax=c(0,0.2))
ggarrange(LY96_Classical_Monocytes[[2]],LY96_Naive_CD8_T_cells[[2]],LY96_Allcells[[2]],nrow=1,ncol=3,common.legend = T)
ggsave("/path/eQTLresults/0_plots/LY96_example.pdf",width=6.2,height=2.7)

## check that eQTL is not sig in all cells
## check nominal p-val thresholdand whether theres a sig eQTL 
allcells_nominal_qval <- fread("/path/eQTLresults/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv")
allcells_nominal_qval[allcells_nominal_qval$phenotype_id == "ENSG00000154589",] 
allcells_nominal <- fread(paste0("grep ENSG00000154589 ","/path/eQTLresults/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.8.tsv"))
colnames(allcells_nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
allcells_nominal[allcells_nominal$variant_id == "chr8_73989162_G_A",]


###ENSG00000153283 - chr3_111535372_A (ref) _C (alt)
CD96_Naive_B_cells <- getadjexp("Naive_B_cells","ENSG00000153283","chr3_111535372_A_C","blue")
CD96_Memory_B_cells <- getadjexp("Memory_B_cells","ENSG00000153283","chr3_111535372_A_C","blue")
CD96_Allcells <- getadjexp("All","ENSG00000153283","chr3_111535372_A_C","blue",minmax=c(0,1.15))
ggarrange(CD96_Naive_B_cells[[2]],CD96_Memory_B_cells[[2]],CD96_Allcells[[2]],nrow=1,ncol=3,common.legend = T)
ggsave("/path/eQTLresults/0_plots/CD96_example.pdf",width=6.2,height=2.7)

min <- min(min(CD96_Naive_B_cells[[1]]$adjusted_expr),min(CD96_Memory_B_cells[[1]]$adjusted_expr),min(CD96_Allcells[[1]]$adjusted_expr))

CD96_Naive_B_cells <- getadjexp("Naive_B_cells","ENSG00000153283","chr3_111535372_A_C","blue",min)
CD96_Memory_B_cells <- getadjexp("Memory_B_cells","ENSG00000153283","chr3_111535372_A_C","blue",min)
CD96_Allcells <- getadjexp("All","ENSG00000153283","chr3_111535372_A_C","blue",min)
ggarrange(CD96_Naive_B_cells[[2]],CD96_Memory_B_cells[[2]],CD96_Allcells[[2]],nrow=1,ncol=3,common.legend = T)

ggsave("/path/eQTLresults/0_plots/CD96_example.pdf",width=6.3,height=3)

