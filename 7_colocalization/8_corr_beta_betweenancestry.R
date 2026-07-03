###ploting correlation between beta of eQTL between ancestries 
##will try both lead coloc SNP and lead eQTL
library(ggplot2)
library(ggpubr)

GWAS_ID="GCST90270940"
DIR_MAIN="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs"

coloc <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_sigresults_",GWAS_ID,"_checksigeQTL.csv"))

coloc_ancestry <- coloc[c("cell_type","gene_id","lead_H4_variant","gene_symbol","lead_H4_variant_slope","lead_H4_variant_pval","lead_H4_variant_GTF_REF","lead_H4_variant_GTF_ALT")]
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
    nominal=fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  }else{
    nominal=fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  }
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  coloc_ancestry$MAIN_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  coloc_ancestry$MAIN_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  
  ###find B in AFR
  nominal<- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  if(length(nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]) ==0){
    coloc_ancestry$AFR_lead_H4_variant_slope[i] <- NA
    coloc_ancestry$AFR_lead_H4_variant_AF[i] <- NA
  } else{
    coloc_ancestry$AFR_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
    coloc_ancestry$AFR_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  }
  
  ###find B in EUR
  nominal<- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  if(length(nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]) ==0){
    coloc_ancestry$EUR_lead_H4_variant_slope[i] <- NA
    coloc_ancestry$EUR_lead_H4_variant_AF[i] <- NA
  }else{
    coloc_ancestry$EUR_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
    coloc_ancestry$EUR_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  }
  
  ###find B in SAS
  nominal<- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  if(length(nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]) ==0){
    coloc_ancestry$SAS_lead_H4_variant_slope[i] <- NA
    coloc_ancestry$SAS_lead_H4_variant_AF[i] <- NA
  }else{
    coloc_ancestry$SAS_lead_H4_variant_slope[i] <- nominal$slope[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
    coloc_ancestry$SAS_lead_H4_variant_AF[i] <- nominal$af[nominal$phenotype_id == gene & nominal$variant_id == colocsnp]
  }

}


write.csv(coloc_ancestry,paste0(DIR_MAIN,'/1_csvfiles/coloc_ancestry_slopeaf.csv'),row.names=F)
coloc_ancestry <- read.csv(paste0(DIR_MAIN,'/1_csvfiles/coloc_ancestry_slopeaf.csv'))

###plotting scatter 
library(ggplot2)
library(ggpubr)
library(viridis)

all_vals <- c(
  abs(coloc_ancestry$AFR_lead_H4_variant_AF-coloc_ancestry$EUR_lead_H4_variant_AF),
  abs(coloc_ancestry$EUR_lead_H4_variant_AF-coloc_ancestry$SAS_lead_H4_variant_AF),
  abs(coloc_ancestry$SAS_lead_H4_variant_AF-coloc_ancestry$EUR_lead_H4_variant_AF)
)

global_min <- min(all_vals, na.rm = TRUE)
global_max <- max(all_vals, na.rm = TRUE)

p1 <- ggplot(coloc_ancestry, aes(x=AFR_lead_H4_variant_slope, y=EUR_lead_H4_variant_slope,color=abs(AFR_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max))+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+xlim(-0.5,0.5)+ylim(-0.5,0.5)+ stat_cor(method = "pearson", label.x = -0.5, label.y = 0.4)
res <- cor.test(coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F ,]$AFR_lead_H4_variant_slope, coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F ,]$EUR_lead_H4_variant_slope, method = "pearson")
res[["p.value"]]

p2 <- ggplot(coloc_ancestry, aes(x=EUR_lead_H4_variant_slope, y=SAS_lead_H4_variant_slope,color=abs(EUR_lead_H4_variant_AF-SAS_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max))+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+xlim(-0.5,0.5)+ylim(-0.5,0.5)+ stat_cor(method = "pearson", label.x = -0.5, label.y = 0.4)
res <- cor.test(coloc_ancestry[is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F ,]$EUR_lead_H4_variant_slope, coloc_ancestry[is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F ,]$SAS_lead_H4_variant_slope, method = "pearson")
res[["p.value"]]

p3 <- ggplot(coloc_ancestry, aes(x=SAS_lead_H4_variant_slope, y=AFR_lead_H4_variant_slope,color=abs(SAS_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max))+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+xlim(-0.5,0.5)+ylim(-0.5,0.5)+ stat_cor(method = "pearson", label.x = -0.5, label.y = 0.4)
res <- cor.test(coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F ,]$AFR_lead_H4_variant_slope, coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F ,]$SAS_lead_H4_variant_slope, method = "pearson")
res[["p.value"]]

ggarrange(p1,p2,p3,common.legend = T,nrow=1,ncol=3)


p1 <- ggplot(coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F,], aes(x=AFR_lead_H4_variant_slope, y=EUR_lead_H4_variant_slope,color=abs(AFR_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max),name="AF difference")+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+ stat_cor(method = "pearson")+xlab("eQTL beta in AFR") + ylab("eQTL beta in EUR")
p2 <- ggplot(coloc_ancestry[is.na(coloc_ancestry$EUR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F,], aes(x=EUR_lead_H4_variant_slope, y=SAS_lead_H4_variant_slope,color=abs(EUR_lead_H4_variant_AF-SAS_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max),name="AF difference")+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+ stat_cor(method = "pearson")+xlim(-0.3,0.3)+ylim(-0.48,0.3)+xlab("eQTL beta in EUR") + ylab("eQTL beta in SAS")
p3 <- ggplot(coloc_ancestry[is.na(coloc_ancestry$AFR_lead_H4_variant_slope) == F & is.na(coloc_ancestry$SAS_lead_H4_variant_slope) == F,], aes(x=SAS_lead_H4_variant_slope, y=AFR_lead_H4_variant_slope,color=abs(SAS_lead_H4_variant_AF-EUR_lead_H4_variant_AF))) +theme_classic()+scale_color_viridis_c(limits = c(global_min, global_max),name="AF difference")+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")+geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + geom_point()+ stat_cor(method = "pearson")+xlab("eQTL beta in SAS") + ylab("eQTL beta in AFR")
ggarrange(p1,p2,p3,common.legend = T,nrow=1,ncol=3,legend="right")

ggsave(paste0(DIR_MAIN,"/0_plots/ancestry_beta_correlation_coloc.pdf"),height=3.5,width=12)
