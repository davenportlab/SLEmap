###ploting correlation between beta of eQTL between ancestries using lead coloc snp - with locusbreaker gwas windows
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)

coloc <- read.csv(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/1_csvfiles/coloc_sigresults_GCST90270940_checksigeQTL.csv"))

coloc_ancestry <- coloc[c("cell_type","gene_id","lead_H4_variant","lead_H4_SNP.PP.H4","gene_symbol","lead_H4_variant_slope","lead_H4_variant_pval","lead_H4_variant_GTF_REF","lead_H4_variant_GTF_ALT")]
coloc_ancestry$MAIN_lead_H4_variant_slope <- NA
coloc_ancestry$MAIN_lead_H4_variant_AF <- NA
coloc_ancestry$AFR_lead_H4_variant_slope <- NA
coloc_ancestry$AFR_lead_H4_variant_AF <- NA
coloc_ancestry$EUR_lead_H4_variant_slope <- NA
coloc_ancestry$EUR_lead_H4_variant_AF <- NA
coloc_ancestry$SAS_lead_H4_variant_slope <- NA
coloc_ancestry$SAS_lead_H4_variant_AF <- NA


for(i in 1:nrow(coloc)){
  print(paste0(i,"/",nrow(coloc)))
  gene <- coloc$gene_id[i]
  gene_chr=coloc$lead_snp_chr[i]
  gene_chr <- gsub("chr", "",gene_chr)
  celltype <- coloc$cell_type[i]
  colocsnp <- paste0(coloc$lead_H4_variant[i],"_",coloc$lead_H4_variant_GTF_REF[i],"_",coloc$lead_H4_variant_GTF_ALT[i])
  
  ##find af in all
  if(celltype == "All"){
    nominal=fread(paste0("grep ",gene," ","/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  }else{
    nominal=fread(paste0("grep ",gene," ","/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  }
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  coloc_ancestry$MAIN_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  coloc_ancestry$MAIN_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  
  ###find B in AFR
  nominal<- fread(paste0("grep ",gene," ","/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  if(length(nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]) ==0){
    coloc_ancestry$AFR_lead_H4_variant_slope[i] <- NA
    coloc_ancestry$AFR_lead_H4_variant_AF[i] <- NA
  } else{
    coloc_ancestry$AFR_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
    coloc_ancestry$AFR_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  }
  
  ###find B in EUR
  nominal<- fread(paste0("grep ",gene," ","/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  if(length(nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]) ==0){
    coloc_ancestry$EUR_lead_H4_variant_slope[i] <- NA
    coloc_ancestry$EUR_lead_H4_variant_AF[i] <- NA
  }else{
    coloc_ancestry$EUR_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
    coloc_ancestry$EUR_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  }
  
  ###find B in SAS
  nominal<- fread(paste0("grep ",gene," ","/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  if(length(nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]) ==0){
    coloc_ancestry$SAS_lead_H4_variant_slope[i] <- NA
    coloc_ancestry$SAS_lead_H4_variant_AF[i] <- NA
  }else{
    coloc_ancestry$SAS_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
    coloc_ancestry$SAS_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  }
  
}

write.csv(coloc_ancestry,'/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/1_csvfiles/coloc_ancestry_slopeaf.csv',row.names=F)
coloc_ancestry <- read.csv('/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/1_csvfiles/coloc_ancestry_slopeaf.csv')

View(coloc_ancestry[coloc_ancestry$AFR_lead_H4_variant_slope * coloc_ancestry$EUR_lead_H4_variant_slope < 0,])
View(coloc_ancestry[coloc_ancestry$AFR_lead_H4_variant_slope * coloc_ancestry$SAS_lead_H4_variant_slope < 0,])
View(coloc_ancestry[coloc_ancestry$EUR_lead_H4_variant_slope * coloc_ancestry$SAS_lead_H4_variant_slope < 0,])

###plotting scatter 
library(ggplot2)
library(ggpubr)
library(viridis)

# difference
all_vals <- c(
  abs(coloc_ancestry$AFR_lead_H4_variant_AF-coloc_ancestry$EUR_lead_H4_variant_AF),
  abs(coloc_ancestry$EUR_lead_H4_variant_AF-coloc_ancestry$SAS_lead_H4_variant_AF),
  abs(coloc_ancestry$SAS_lead_H4_variant_AF-coloc_ancestry$EUR_lead_H4_variant_AF)
)


global_min <- min(all_vals, na.rm = TRUE)
global_max <- max(all_vals, na.rm = TRUE)

p1 <- ggplot(coloc_ancestry, aes(x=AFR_lead_H4_variant_slope, y=EUR_lead_H4_variant_slope,color=abs(AFR_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max))+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+xlim(-0.5,0.5)+ylim(-0.5,0.5)+ stat_cor(method = "pearson", label.x = -0.5, label.y = 0.4)
p2 <- ggplot(coloc_ancestry, aes(x=EUR_lead_H4_variant_slope, y=SAS_lead_H4_variant_slope,color=abs(EUR_lead_H4_variant_AF-SAS_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max))+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+xlim(-0.5,0.5)+ylim(-0.5,0.5)+ stat_cor(method = "pearson", label.x = -0.5, label.y = 0.4)
p3 <- ggplot(coloc_ancestry, aes(x=SAS_lead_H4_variant_slope, y=AFR_lead_H4_variant_slope,color=abs(SAS_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max))+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+xlim(-0.5,0.5)+ylim(-0.5,0.5)+ stat_cor(method = "pearson", label.x = -0.5, label.y = 0.4)
ggarrange(p1,p2,p3,common.legend = T,nrow=1,ncol=3)

cor.test(coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F,]$AFR_lead_H4_variant_slope, coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F,]$EUR_lead_H4_variant_slope, method = 'pearson')

p1 <- ggplot(coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F,], aes(x=AFR_lead_H4_variant_slope, y=EUR_lead_H4_variant_slope,color=abs(AFR_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max),name="AF difference")+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point(alpha=0.75)+ stat_cor(method = "pearson")+xlab("eQTL beta in AFR") + ylab("eQTL beta in EUR")
p2 <- ggplot(coloc_ancestry[is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F,], aes(x=EUR_lead_H4_variant_slope, y=SAS_lead_H4_variant_slope,color=abs(EUR_lead_H4_variant_AF-SAS_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max),name="AF difference")+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point(alpha=0.75)+ stat_cor(method = "pearson")+xlim(-0.3,0.3)+ylim(-0.48,0.3)+xlab("eQTL beta in EUR") + ylab("eQTL beta in SAS")
p3 <- ggplot(coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F,], aes(x=SAS_lead_H4_variant_slope, y=AFR_lead_H4_variant_slope,color=abs(SAS_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max),name="AF difference")+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point(alpha=0.75)+ stat_cor(method = "pearson")+xlab("eQTL beta in SAS") + ylab("eQTL beta in AFR")
ggarrange(p2,p1,p3,common.legend = T,nrow=1,ncol=3,legend="right")

ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/Beta_correlation_coloc.pdf",height=3.5,width=12)


###making figs of examples
##PLCL1
library(data.table)
metadata <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/slemap/phenotype/combined_data_tidied_sc_used_w_KING_28_08_25.csv")

getplot_ancestry <- function(celltype,Gene,SNP,runSNP){
  setwd("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/3_eQTL_prep/SNPs")
  if(celltype =="All"){
    eQTLrun <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/"
    norm_data_dir <- paste0(eQTLrun,"/results/norm_data/dMean__",celltype,"_all/")
    mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
    covariates <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
    optimPCs <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
    
  }else{
    eQTLrun <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/"
    norm_data_dir <- paste0(eQTLrun,celltype,"/results/norm_data/dMean__",celltype,"_all/")
    mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
    covariates <- read.table(paste0(eQTLrun,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
    optimPCs <- read.table(paste0(eQTLrun,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
  }
  
  SNPs <- setNames(data.frame(unlist(strsplit(SNP, "_"))[1], unlist(strsplit(SNP, "_"))[2], unlist(strsplit(SNP, "_"))[2]), c('chr_name', 'chrom_start','chrom_end'))
  if(runSNP == T){
    system2("bash", args = c("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/3_eQTL_prep/SNPs/genotype.sh", SNPs$chr_name, SNPs$chrom_start,SNPs$chrom_end,SNP))
  }
  ## get genotype
  indiv <- read.table(paste0(SNP,".012.indv"),sep="\t")
  genotype <- read.table(paste0(SNP,".012.gz"),sep="\t")
  #system2("rm",args= c(paste0(SNP,"*")))
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
  mean_use$Genotype <- factor(mean_use$Genotype,levels=c(0,1,2))
  
  # Plot adjusted expression by genotype
  p2 <- ggplot(mean_use[!is.na(mean_use$Genotype),],
               aes(y = adjusted_expr, x = factor(Genotype), color = ancestry)) + 
    
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
    
    geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 1), size = 0.9) +
    
    ggtitle(paste0("SLEmap:", Gene, "-", SNP, ", n=",
                   nrow(mean_use[mean_use$Genotype != -1,]))) +
    
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    
    facet_wrap(~ ancestry, nrow=1, ncol=4) +
    
    theme_classic()
  
  return(p2)
  
}

coloc <- read.csv(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/1_csvfiles/coloc_sigresults_GCST90270940_checksigeQTL.csv"))

coloc[coloc$gene_symbol == "PLCL1" & coloc$cell_type == "EM_CD8_T_cells",]
getplot_ancestry(celltype="EM_CD8_T_cells",Gene="ENSG00000115896",SNP="chr2_198035639_G_A",runSNP=F)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/forpaper/PLCL1_EM_CD8_T_cells_ancestry.pdf",width=7,height=3)

coloc[coloc$gene_symbol == "IL10RA" & coloc$cell_type == "CM_CD4_T_cells",]
getplot_ancestry(celltype="CM_CD4_T_cells",Gene="ENSG00000110324",SNP="chr11_118009424_C_T",runSNP=F)

coloc[coloc$gene_symbol == "CD58" & coloc$cell_type == "Memory_B_cells",]
getplot_ancestry(celltype="Memory_B_cells",Gene="ENSG00000116815",SNP="chr1_116529902_T_C",runSNP=F)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/forpaper/CD58_Memory_B_cells_ancestry.pdf",width=7,height=3)


slope <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/1_csvfiles/coloc_ancestry_slopeaf.csv")
slope_opposite <- slope[slope$AFR_lead_H4_variant_slope * slope$EUR_lead_H4_variant_slope < 0,]
slope_opposite <- slope_opposite[is.na(slope_opposite$cell_type) == F,]
slope_opposite[slope_opposite$gene_symbol == "PPP1R14B" & slope_opposite$cell_type == "All",]
coloc[coloc$gene_symbol == "PPP1R14B" & coloc$cell_type == "All",]
getplot_ancestry(celltype="All",Gene="ENSG00000173457",SNP="chr11_64360951_G_T",runSNP=F)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/forpaper/PPP1R14B_All_cells_ancestry.pdf",width=7,height=3)


