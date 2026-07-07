###make eQTL plots for slemap vs onek1k example 
library(tidyverse)
library(ggpubr)

metadata <- read.csv("/path/KING.csv")

getplot <- function(celltype,Gene,SNP){
  setwd("/path/eQTL_prep/SNPs")
  if(celltype =="All"){
    eQTLrun <- "/path/eQTLmapping/"
    norm_data_dir <- paste0(eQTLrun,"/results/norm_data/dMean__",celltype,"_all/")
    mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
    covariates <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
    optimPCs <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
    
  }else{
    eQTLrun <- "/path/eQTLmapping/"
    norm_data_dir <- paste0(eQTLrun,celltype,"/results/norm_data/dMean__",celltype,"_all/")
    mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
    covariates <- read.table(paste0(eQTLrun,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
    optimPCs <- read.table(paste0(eQTLrun,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
  }
  
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
  sampleIDs <- as.character(covariates[1,])
  row.names(covariates) <- paste0(covariates$V1,covariates$V2)
  covariates <- covariates[-1,-c(1:2)]    
  colnames(covariates) <- sampleIDs[1:ncol(covariates)]
  covariates <- mutate_all(as.data.frame(t(covariates)), function(x) as.numeric(as.character(x)))
  
  mean_use <- merge(mean_use,covariates,by.x="sampleID",by.y=0)
  
  # Get optim phenotype PCs and 5 genotypePCs
  
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
  mean_use$ancestry <- metadata$Ancestry[match(mean_use$sampleID,metadata$WGS_ID)]
  mean_use$ancestry <- factor(mean_use$ancestry,levels=c("AFR","EUR","SAS"))
  
  # Plot adjusted expression by genotype
  p2 <- ggplot(mean_use[mean_use$Genotype != -1,], aes(y = adjusted_expr, x = Genotype, group=Genotype,color=factor(Genotype))) +
    geom_boxplot(outlier.shape = NA)+ scale_x_continuous(breaks = c(0, 1, 2)) +scale_color_manual(values=c("#FBB042","#FBB042","#FBB042"))+
    geom_jitter(width=0.25,size=0.9,alpha=0.7)+ggtitle(paste0("SLEmap:",Gene,"-",SNP,", n=",nrow(mean_use[mean_use$Genotype != -1,])))+theme_classic()+theme(legend.position="none")
  
  
  return(p2)
  
}

getplot_onek1k <- function(celltype,Gene,SNP){
  setwd("/path/Onek1k/Genotypes/SNPs")
  if(celltype =="All"){
    eQTLrun <- "/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/"
    norm_data_dir <- paste0(eQTLrun,"/results/norm_data/dMean__",celltype,"_all/")
    mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
    covariates <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
    optimPCs <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
    
  }else{
    eQTLrun <- "/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs/"
    norm_data_dir <- paste0(eQTLrun,celltype,"/results/norm_data/dMean__",celltype,"_all/")
    mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
    covariates <- read.table(paste0(eQTLrun,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
    optimPCs <- read.table(paste0(eQTLrun,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
  }
  
  SNPs <- setNames(data.frame(gsub("chr","",unlist(strsplit(SNP, "_"))[1]), unlist(strsplit(SNP, "_"))[2], unlist(strsplit(SNP, "_"))[2]), c('chr_name', 'chrom_start','chrom_end'))
  ## get genotype
  system2("bash", args = c("/path/Onek1k/Genotypes/SNPs/genotype.sh", SNPs$chr_name, SNPs$chrom_start,SNPs$chrom_end,SNP))
  indiv <- read.table(paste0(SNP,".012.indv"),sep="\t")
  genotype <- read.table(paste0(SNP,".012.gz"),sep="\t")
  system2("rm",args= c(paste0(SNP,"*")))
  indiv$geno <- genotype$V2
  
  ## plot (mean) gene exp
  mean_use <- as.data.frame(t(mean[Gene,]))
  mean_use$sampleID <- sub(".*_([0-9]+_[0-9]+)_\\1$", "\\1", row.names(mean_use))
  mean_use$geno <- indiv$geno[match(mean_use$sampleID, indiv$V1)]
  colnames(mean_use) <- c("Mean_exp","sampleID","Genotype")
  
  ##input covariates
  sampleIDs <- as.character(covariates[1,])
  row.names(covariates) <- paste0(covariates$V1,covariates$V2)
  covariates <- covariates[-1,-c(1:2)]    
  colnames(covariates) <- sampleIDs[1:ncol(covariates)]
  covariates <- mutate_all(as.data.frame(t(covariates)), function(x) as.numeric(as.character(x)))
  
  mean_use <- merge(mean_use,covariates,by.x="sampleID",by.y=0)
  
  # Get optim phenotype PCs and 5 genotypePCs
  
  optimPCs <- as.numeric(optimPCs[1,1])
  
  genotypePCs <- paste0("GenotypePC",1:4)
  phenotypePCs <- paste0("PhenotypePC",1:optimPCs)
  
  ##tensorQTL imputes genotypes for missing values: they use the mean genotype (not integer)
  mean_use$Genotype_imputed <- mean_use$Genotype
  mean_use$Genotype_imputed[mean_use$Genotype_imputed == -1] <- mean(mean_use$Genotype_imputed[mean_use$Genotype_imputed != -1])
  
  # get full model: the coefficient for genotypes should match the slope from tensorQTL
  null_model <- lm(mean_use[c("Mean_exp", "Genotype_imputed", genotypePCs, phenotypePCs)]) 
  coefs <- null_model[["coefficients"]][append(genotypePCs,phenotypePCs)]
  mean_use$adjusted_expr <- mean_use$Mean_exp - rowSums(sweep(mean_use[append(genotypePCs,phenotypePCs)], 2, coefs, FUN = "*"))
  
  # Plot adjusted expression by genotype
  p2 <- ggplot(mean_use[mean_use$Genotype != -1,], aes(y = adjusted_expr, x = Genotype, group=Genotype,color=factor(Genotype))) + 
    geom_boxplot(outlier.shape = NA)+ scale_x_continuous(breaks = c(0, 1, 2)) +scale_color_manual(values=c("#4C79A9","#4C79A9","#4C79A9"))+
    geom_jitter(width=0.25,size=0.9,alpha=0.7)+ggtitle(paste0("Onek1k:",Gene,"-",SNP,", n=",nrow(mean_use[mean_use$Genotype != -1,])))+theme_classic()+theme(legend.position="none")
  
  return(p2)
  
}

##get coloc results for examples
coloc <- read.csv(paste0("/path/onek1k_locus_breaker_coloc/","1_csvfiles/SLEmap_Onek1kgroup.csv"),row.names=1)


coloc[coloc$gene_symbol == "ITGAM",]
ITGAM_mono_slemap <- getplot("Classical_Monocytes","ENSG00000169896","chr16_31301932_C_T")
ITGAM_mono_onek1k <- getplot_onek1k("Classical_Monocytes","ENSG00000169896","chr16_31301932_C_T")
ggarrange(ITGAM_mono_slemap,ITGAM_mono_onek1k)
ggsave("/path/onek1k_locus_breaker_coloc/0_plots/forpaper/ITGAM_Classical_Monocytes.pdf",width=4.5,height=3)

coloc[coloc$gene_symbol == "IKZF3",]
IKZF3_naivecd8_slemap <- getplot("Naive_CD8_T_cells","ENSG00000161405","chr17_39820216_C_T")
IKZF3_naivecd8_onek1k <- getplot_onek1k("Naive_CD8_T_cells","ENSG00000161405","chr17_39820216_C_T")
ggarrange(IKZF3_naivecd8_slemap,IKZF3_naivecd8_onek1k)
ggsave("/path/onek1k_locus_breaker_coloc/0_plots/forpaper/IKZF3_Naive_CD8_T_cells.pdf",width=4.5,height=3)

ggarrange(ITGAM_mono_slemap,ITGAM_mono_onek1k,IKZF3_naivecd8_slemap,IKZF3_naivecd8_onek1k,ncol=2,nrow=2,common.legend=T)
