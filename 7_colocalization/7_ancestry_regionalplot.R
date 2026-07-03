### regional plots for GWAS, main eQTL, eQTL by ancestry for coloc signals
library(ggplot2)
library(ggpubr)
library(rlist)
library(data.table)
library(dplyr)
library(forcats)
library(scales)

### grey dots in plot means that the LD could not be calculated because lead H4 snp had too low maf

## FOR GWAS
GWAS_input= "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/14_colocalization/inputs/gwas/"
#GWAS_ID="GCST003156"
GWAS_ID="GCST90270940"

##get maf info
maf <- fread("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/wl2/wgs/01_5.variants_for_eQTLs/2.eQTL_hwe_0.000001_maf_0.05_miss_0.95/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.afreq")
AFR_maf <- fread("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_AFR_n102_ancestrymaf_0.05.afreq")
EUR_maf <- fread("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_EUR_n81_ancestrymaf_0.05.afreq")
SAS_maf <- fread("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62_ancestrymaf_0.05.afreq")

##GET coloc results
DIR_MAIN="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs"

coloc <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_sigresults_",GWAS_ID,"_checksigeQTL.csv"))

GWAS_INPUT <- fread (paste0(GWAS_input,GWAS_ID,"/",GWAS_ID,"_for_eqtl.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::rename(., variant_id_v2=variant_id) %>%
  dplyr::mutate(variant_id = paste0(chr, "_", pos))
GWAS_INPUT <- GWAS_INPUT[!((GWAS_INPUT$chr == "chr6") & (GWAS_INPUT$pos > 25000000)&(GWAS_INPUT$pos < 34000000)),] 

GWAS_INPUT$pos <- as.numeric(GWAS_INPUT$pos)
GWAS_INPUT$p_value <- as.numeric(GWAS_INPUT$p_value)
locusbreaker <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/gwashit_locusbreaker_forinput.csv")
coloc_group <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_genegroup_",GWAS_ID,"_checksigeQTL.csv"),nrow=F)

##make output dir
if (!dir.exists(paste0(DIR_MAIN,"/0_plots/ancestry_LD_onlyoverlap"))) {
  dir.create(paste0(DIR_MAIN,"/0_plots/ancestry_LD_onlyoverlap"))
}

if (!dir.exists(paste0(DIR_MAIN,"/0_plots/ancestry_LD_allSNPs"))) {
  dir.create(paste0(DIR_MAIN,"/0_plots/ancestry_LD_allSNPs"))
}

for(i in 1:nrow(coloc)){

  colocvariant <- coloc$lead_H4_variant[i]
  coloc_genesymbol <- coloc$gene_symbol[i]
  ALT <- coloc$lead_H4_variant_GTF_ALT[i]
  REF <- coloc$lead_H4_variant_GTF_REF[i]
  print(coloc_genesymbol)
  celltype <- coloc$cell_type[i]
  print(celltype)
  gene=coloc$gene_id[i]
  gene_symbol=coloc$gene_symbol[i]
  gene_chr=coloc$lead_snp_chr[i]
  gene_chr <- gsub("chr", "",gene_chr)
  eQTL_group <- coloc_group$group[(coloc_group$gene == gene)&(coloc_group$celltype == celltype)]
  
  ##get GWAS
  GWAShit = coloc$gwas_hit[i]
  GWAShit_chr = sapply(strsplit(GWAShit, "_"), `[`, 1)
  GWAShit_pos = as.numeric(sapply(strsplit(GWAShit, "_"), `[`, 2))
  pos_start = locusbreaker$locusStart[locusbreaker$variant_id == GWAShit]
  pos_end = locusbreaker$locusEnd[locusbreaker$variant_id == GWAShit]
  
  GWAS_input_plot <- GWAS_INPUT[GWAS_INPUT$chr == GWAShit_chr,]
  GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos < pos_end,]
  GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos > pos_start,]
  
  GWAS_input_plot$color <- "notlead_colocsnp"
  GWAS_input_plot$color[GWAS_input_plot$variant_id ==colocvariant] <- "lead_colocsnp"
  GWAS_input_plot <- GWAS_input_plot %>% arrange(desc(color))
  
  ##get LD 
  system2("bash", args = c("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/20_ancestry_coloc/6_LDbyancestry.sh",paste0(colocvariant,"_",REF,"_",ALT),gene_chr, pos_start,pos_end))
  all_ld <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/0_calcLD/all/",colocvariant,"_",REF,"_",ALT,".vcor"))

  ## get eQTL
  if(celltype == "All"){
    EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
  }else{
    EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
  }
  
  independent=read.table(paste0(EQTL_input,"Cis_eqtls_independent.tsv"),header=T)
  independent_n <- nrow(independent[independent$phenotype_id == gene,])
  
  ##get nominal p-val threshold for all eQTL
  pvalthres <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/conditionaleQTL.csv")
  pvalthres <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
  pvalthres <- pvalthres[1]
  
  plot <- list()
  plot_overlap <- list()
  if(independent_n > 1){
    if(celltype == "All"){
      indeplist <- list.files(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/2_indep_coloc/",celltype,"/nominal_p/"),gene,full.names = TRUE)
      pvalthres <- read.table("//lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv",header=T,fill=T)
      pvalthres <- pvalthres$pval_nominal_threshold[pvalthres$phenotype_id == gene]
    }else{
      indeplist <- list.files(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/2_indep_coloc/",celltype,"/nominal_p/"),gene,full.names = TRUE)
      pvalthres <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/conditionaleQTL.csv")
      pvalthres <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
      pvalthres <- pvalthres[1]
    }
    
    for(m in 1:length(indeplist)){
      indepresults <- read.csv(indeplist[m],header=T)
      indepresults$POS <- as.numeric(sapply(strsplit(indepresults$variant_id, "_"), `[`, 2))
      indepresults$chr_pos <- sapply(strsplit(indepresults$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
      
      overlapping_variants <- intersect(indepresults$chr_pos,GWAS_input_plot$variant_id)
      indepresults$overlap <- "not_overlap"
      indepresults$overlap[indepresults$chr_pos %in% overlapping_variants] <- "overlap"
      indepresults$coloc_leadsnp <- "notlead_colocsnp"
      indepresults$coloc_leadsnp[indepresults$chr_pos == colocvariant] <- "lead_colocsnp"
      indepresults <- indepresults %>% arrange(desc(coloc_leadsnp))
      
      indepresults$LD <- all_ld$V7[match(indepresults$variant_id,all_ld$V6)]
      indepresults$LD <- fct_rev(cut(indepresults$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                                labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
      indepresults$LD[indepresults$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"
      
      GWAS_input_plot$overlap <- "not_overlap"
      GWAS_input_plot$overlap[GWAS_input_plot$variant_id %in% overlapping_variants] <- "overlap"
      GWAS_input_plot$coloc_leadsnp <- "notlead_colocsnp"
      GWAS_input_plot$coloc_leadsnp[GWAS_input_plot$variant_id == colocvariant] <- "lead_colocsnp"
      
      AF <- maf$ALT_FREQS[maf$ID == paste0(colocvariant,"_",REF,"_",ALT)]
      
      if(length(plot) == 0){
        plot <- list.append(plot,ggplot(GWAS_input_plot, aes(x = pos, y = -log10(p_value),shape=coloc_leadsnp)) + geom_point(aes(size=overlap,alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw()+scale_shape_manual(values=c(17,16))+ geom_hline(yintercept=5, linetype="dashed"))# to show all snps in region while having alpha higher in overlapping ones
        
        plot_overlap <- list.append(plot_overlap,ggplot(GWAS_input_plot[GWAS_input_plot$overlap == "overlap",], aes(x = pos, y = -log10(p_value),shape=coloc_leadsnp)) + geom_point(aes(size=overlap,alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw()+scale_shape_manual(values=c(17,16))+ geom_hline(yintercept=5, linetype="dashed"))# to show all snps in region while having alpha higher in overlapping ones
      }
      
      plot <- list.append(plot,ggplot(indepresults, aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol,"(indep signal ",m,") (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed")) # to show all snps in region while having alpha higher in overlapping ones
      
      plot_overlap <- list.append(plot_overlap,ggplot(indepresults[indepresults$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol,"(indep signal ",m,") (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))  + geom_hline(yintercept=-log10(pvalthres), linetype="dashed"))

    }
    
  }else{
    
    nominal <- fread(paste0("grep ",gene," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
    colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
    nominal <- nominal[nominal$phenotype_id == gene,]
    nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
    nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
    
    overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
    nominal$overlap <- "not_overlap"
    nominal$overlap[nominal$chr_pos %in% overlapping_variants] <- "overlap"
    nominal$coloc_leadsnp <- "notlead_colocsnp"
    nominal$coloc_leadsnp[nominal$chr_pos == colocvariant] <- "lead_colocsnp"
    nominal <- nominal %>% arrange(desc(coloc_leadsnp))
    
    nominal$LD <- all_ld$V7[match(nominal$variant_id,all_ld$V6)]
    nominal$LD <- fct_rev(cut(nominal$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                                   labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
    nominal$LD[nominal$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"
    
    GWAS_input_plot$overlap <- "not_overlap"
    GWAS_input_plot$overlap[GWAS_input_plot$variant_id %in% overlapping_variants] <- "overlap"
    GWAS_input_plot$coloc_leadsnp <- "notlead_colocsnp"
    GWAS_input_plot$coloc_leadsnp[GWAS_input_plot$variant_id == colocvariant] <- "lead_colocsnp"
    
    AF <- maf$ALT_FREQS[maf$ID == paste0(colocvariant,"_",REF,"_",ALT)]
    
    ##get nominal p-val threshold for all eQTL
    if(celltype == "All"){
      pvalthres <- read.table("//lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv",header=T,fill=T)
      pvalthres <- pvalthres$pval_nominal_threshold[pvalthres$phenotype_id == gene]
    }else{
      pvalthres <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/conditionaleQTL.csv")
      pvalthres <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
      pvalthres <- pvalthres[1]
    }
    
    plot <- list.append(plot,ggplot(GWAS_input_plot, aes(x = pos, y = -log10(p_value),shape=coloc_leadsnp)) + geom_point(aes(size=overlap,alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw()+scale_shape_manual(values=c(17,16))+ geom_hline(yintercept=5, linetype="dashed"))# to show all snps in region while having alpha higher in overlapping ones
    
    plot_overlap <- list.append(plot_overlap,ggplot(GWAS_input_plot[GWAS_input_plot$overlap == "overlap",], aes(x = pos, y = -log10(p_value),shape=coloc_leadsnp)) + geom_point(aes(size=overlap,alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw()+scale_shape_manual(values=c(17,16))+ geom_hline(yintercept=5, linetype="dashed"))# to show all snps in region while having alpha higher in overlapping ones
    
    plot <- list.append(plot,ggplot(nominal, aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed")) # to show all snps in region while having alpha higher in overlapping ones
    
    plot_overlap <- list.append(plot_overlap,ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed")) # to show all snps in region while having alpha higher in overlapping ones

  }
  
  ###ancestry specific
  ###AFR
  nominal<- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  nominal <- nominal[nominal$phenotype_id == gene,]
  nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
  nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
  
  overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
  nominal$overlap <- "not_overlap"
  nominal$overlap[nominal$chr_pos %in% overlapping_variants] <- "overlap"
  nominal$coloc_leadsnp <- "notlead_colocsnp"
  nominal$coloc_leadsnp[nominal$chr_pos == colocvariant] <- "lead_colocsnp"
  nominal <- nominal %>% arrange(desc(coloc_leadsnp))
  
  LDfile <- paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/0_calcLD/AFR/",colocvariant,"_",REF,"_",ALT,".vcor")
  if(file.exists(LDfile)){
    AFR_ld <- read.table(LDfile)
    nominal$LD <- AFR_ld$V7[match(nominal$variant_id,AFR_ld$V6)]
    nominal$LD <- fct_rev(cut(nominal$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                                labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
    nominal$LD[nominal$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"
  }else{
    nominal$LD <- NA
  }

  pvalthres <- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"))
  if(pvalthres$qval[pvalthres$phenotype_id == gene] < 0.05){
    pvalthres <- pvalthres$pval_nominal_threshold[pvalthres$phenotype_id == gene]
  }else{
    pvalthres <- NA
  }
  
  AF <- AFR_maf$ALT_FREQS[AFR_maf$ID == paste0(colocvariant,"_",REF,"_",ALT)]
                    
  plot <- list.append(plot,ggplot(nominal, aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("AFR: ",celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed")) # to show all snps in region while having alpha higher in overlapping ones
  
  plot_overlap <- list.append(plot_overlap,ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("AFR: ",celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed"))
  
  
  ###EUR
  nominal<- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  
  nominal <- nominal[nominal$phenotype_id == gene,]
  nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
  nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
  
  overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
  nominal$overlap <- "not_overlap"
  nominal$overlap[nominal$chr_pos %in% overlapping_variants] <- "overlap"
  nominal$coloc_leadsnp <- "notlead_colocsnp"
  nominal$coloc_leadsnp[nominal$chr_pos == colocvariant] <- "lead_colocsnp"
  nominal <- nominal %>% arrange(desc(coloc_leadsnp))
  
  LDfile <- paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/0_calcLD/EUR/",colocvariant,"_",REF,"_",ALT,".vcor")
  if(file.exists(LDfile)){
    EUR_ld <- read.table(LDfile)
    nominal$LD <- EUR_ld$V7[match(nominal$variant_id,EUR_ld$V6)]
    nominal$LD <- fct_rev(cut(nominal$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                              labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
    nominal$LD[nominal$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"
  }else{
    nominal$LD <- NA
  }
  
  pvalthres <- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"))
  if(pvalthres$qval[pvalthres$phenotype_id == gene] < 0.05){
    pvalthres <- pvalthres$pval_nominal_threshold[pvalthres$phenotype_id == gene]
  }else{
    pvalthres <- NA
  }
  
  AF <- EUR_maf$ALT_FREQS[EUR_maf$ID == paste0(colocvariant,"_",REF,"_",ALT)]
  
  plot <- list.append(plot,ggplot(nominal, aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("EUR: ",celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed")) # to show all snps in region while having alpha higher in overlapping ones
  
  plot_overlap <- list.append(plot_overlap,ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("EUR: ",celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed"))
  
  ###SAS
  nominal<- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  nominal <- nominal[nominal$phenotype_id == gene,]
  nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
  nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
  
  overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
  nominal$overlap <- "not_overlap"
  nominal$overlap[nominal$chr_pos %in% overlapping_variants] <- "overlap"
  nominal$coloc_leadsnp <- "notlead_colocsnp"
  nominal$coloc_leadsnp[nominal$chr_pos == colocvariant] <- "lead_colocsnp"
  nominal <- nominal %>% arrange(desc(coloc_leadsnp))
  
  LDfile <- paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/0_calcLD/SAS/",colocvariant,"_",REF,"_",ALT,".vcor")
  if(file.exists(LDfile)){
    SAS_ld <- read.table(LDfile)
    nominal$LD <- SAS_ld$V7[match(nominal$variant_id,SAS_ld$V6)]
    nominal$LD <- fct_rev(cut(nominal$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                              labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
    nominal$LD[nominal$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"
  }else{
    nominal$LD <- NA
  }
  
  pvalthres <- fread(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/ManualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"))
  if(pvalthres$qval[pvalthres$phenotype_id == gene] < 0.05){
    pvalthres <- pvalthres$pval_nominal_threshold[pvalthres$phenotype_id == gene]
  }else{
    pvalthres <- NA
  }
  
  AF <- SAS_maf$ALT_FREQS[SAS_maf$ID == paste0(colocvariant,"_",REF,"_",ALT)]
  
  plot <- list.append(plot,ggplot(nominal, aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("SAS: ",celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed")) # to show all snps in region while having alpha higher in overlapping ones
  
  plot_overlap <- list.append(plot_overlap,ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("SAS: ",celltype," ",coloc_genesymbol," (lead coloc SNP ALT AF=",AF,")"))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16)) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed"))  
  
  pdf(paste0(DIR_MAIN,"/0_plots/ancestry_LD_allSNPs/",gene_symbol,"_",celltype,".pdf"),height=13,width=6)
  print(do.call(ggarrange, c(plot, ncol = 1, align = "v",common.legend=T)))
  dev.off()
  
  pdf(paste0(DIR_MAIN,"/0_plots/ancestry_LD_onlyoverlap/",gene_symbol,"_",celltype,".pdf"),height=13,width=6)
  print(do.call(ggarrange, c(plot_overlap, ncol = 1, align = "v",common.legend=T)))
  dev.off()
  
  system2("rm",args= c(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/0_calcLD/*/",colocvariant,"*")))
  
}

