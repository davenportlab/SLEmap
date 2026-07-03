### examine coloc results
library(dplyr)
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggplot2)

##read in coloc results
GWAS_ID <- "GCST90270940"
DIR_MAIN="/path/coloc/outputs"

coloc_results <- read.table(paste0(DIR_MAIN,"/",GWAS_ID,"/all_checksigeQTL_checkallele/coloc_output_with_gene_name.txt"),sep="\t",fill=T,header=T)

## get only significant
sig_coloc <- coloc_results[coloc_results$PP.H4.abf >= 0.8,]

## count how many genes coloced by cell type
sig_coloc %>% group_by(cell_type) %>% tally() %>% arrange(., n)
all_cell_types <- unique(sig_coloc$cell_type)
length(unique(sig_coloc$gene_symbol))

## get cell group info
sig_coloc$cellgroup <- NA
sig_coloc$cellgroup[sig_coloc$cell_type %in% c("CD56Bright_NK_cells","CD56Dim_NK_cells")] <- "NK"
sig_coloc$cellgroup[sig_coloc$cell_type %in% c("Classical_Monocytes")] <- "Mono"
sig_coloc$cellgroup[sig_coloc$cell_type %in% c("CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells")]<- "CD4_T"
sig_coloc$cellgroup[sig_coloc$cell_type %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
sig_coloc$cellgroup[sig_coloc$cell_type %in% c("DN_T_cells")] <- "Other_T"
sig_coloc$cellgroup[sig_coloc$cell_type %in% c("Memory_B_cells","Naive_B_cells")] <- "B"
sig_coloc$cellgroup[sig_coloc$cell_type %in% c("All")] <- "All"

sig_coloc$celltype_forplots <- gsub("_", " ", sig_coloc$cell_type)
sig_coloc$celltype_forplots <- factor(sig_coloc$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
sig_coloc$cellgroup_forplots <- gsub("_", " ", sig_coloc$cellgroup)
sig_coloc$cellgroup_forplots <- factor(sig_coloc$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))

## make dataframe for plotting
dat <- sig_coloc %>% 
  mutate(., id=paste0(gene_symbol, "_", gwas_hit)) %>%
  dplyr::arrange(., gene_symbol) %>%
  dplyr::mutate(., Ord2=nrow(sig_coloc):1)

## get eQTL slope and p-value
dat$lead_H4_variant_slope <- NA
dat$lead_H4_variant_pval <- NA
dat$lead_H4_variant_GTF_REF <- NA
dat$lead_H4_variant_GTF_ALT <- NA

EQTL_input="/path/"

for(i in 1:nrow(dat)){
  celltype_use <- dat$cell_type[i]
  gene_id_use <- dat$gene_id[i]
  lead_snp_use <- dat$lead_snp[i]
  lead_H4_variant_use <- dat$lead_H4_variant[i]
  eQTL_group <- dat$eQTL_group[i]
  
  if(eQTL_group == "B"){
    nominal <- fread(paste0(EQTL_input,celltype_use,"/cis_nominal_",celltype_use,".multiple_eSNPs.txt.gz"))
    nominal <- nominal[(nominal$phenotype_id == gene_id_use) & (nominal$variant_id == lead_H4_variant_use) & (gsub("(chr.*\\_\\d+)\\_.*\\_.*", "\\1", nominal$independent_variant) == lead_snp_use),]
    dat$lead_H4_variant_slope[i] <- nominal$slope
    dat$lead_H4_variant_pval[i] <- nominal$pval_nominal
    dat$lead_H4_variant_GTF_ALT[i] <- nominal$GTF_REF ### these were flipped when coloc input files were created. REF/ALT Checked with the vcf and gtf files. the slope is not flipped. checked with nominal p-val eQTL result files
    dat$lead_H4_variant_GTF_REF[i] <- nominal$GTF_ALT ### these were flipped when coloc input files were created. REF/ALT Checked with the vcf and gtf files. the slope is not flipped. checked with nominal p-val eQTL result files
    
  }else{
    nominal <- fread(paste0(EQTL_input,celltype_use,"/cis_nominal_",celltype_use,".txt.gz"))
    nominal <- nominal[(nominal$phenotype_id == gene_id_use) & (nominal$chr_pos == lead_H4_variant_use),]
    dat$lead_H4_variant_slope[i] <- nominal$slope
    dat$lead_H4_variant_pval[i] <- nominal$pval_nominal
    dat$lead_H4_variant_GTF_ALT[i] <- nominal$GTF_REF ### these were flipped when coloc input files were created. REF/ALT Checked with the vcf and gtf files. the slope is not flipped. checked with nominal p-val eQTL result files
    dat$lead_H4_variant_GTF_REF[i] <- nominal$GTF_ALT ### these were flipped when coloc input files were created. REF/ALT Checked with the vcf and gtf files. the slope is not flipped. checked with nominal p-val eQTL result files
  }
}
rm(nominal)

### check if lead H4 variant alternative allele matches the gwas effect (risk) allele
GWAS_input= "/path/coloc/inputs/gwas/"
GWAS_INPUT <- fread (paste0(GWAS_input,GWAS_ID,"/",GWAS_ID,"_for_eqtl.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::rename(., variant_id_v2=variant_id) %>%
  dplyr::mutate(variant_id = paste0(chr, "_", pos))

GWAS_INPUT <- GWAS_INPUT[!((GWAS_INPUT$chr == "chr6") & (GWAS_INPUT$pos > 25000000)&(GWAS_INPUT$pos < 34000000)),] 
head(GWAS_INPUT)
dat$lead_H4_variant_GWAS_effect_allele <- GWAS_INPUT$other_allele[match(dat$lead_H4_variant, GWAS_INPUT$variant_id)]  ### these were flipped when coloc input files were created. Checked with the original gwas files.
dat$lead_H4_variant_GWAS_other_allele <- GWAS_INPUT$effect_allele[match(dat$lead_H4_variant, GWAS_INPUT$variant_id)]  ### these were flipped when coloc input files were created. Checked with the original gwas files.
all(dat$lead_H4_variant_GWAS_effect_allele == dat$lead_H4_variant_GTF_ALT) #should all be true, inspect if not, if the alleles are flipped, slope * -1

dat$slope_forplotting <- dat$lead_H4_variant_slope
dat$slope_forplotting[dat$lead_H4_variant_GWAS_effect_allele != dat$lead_H4_variant_GTF_ALT] <- -dat$lead_H4_variant_slope[dat$lead_H4_variant_GWAS_effect_allele != dat$lead_H4_variant_GTF_ALT]

dat$cellgroup_forplots <- factor(dat$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK"))
dat$celltype_forplots <- factor(dat$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells"))
dat$gwashit_forplot <- vapply(strsplit(dat$id, "_"), function(parts) {
  paste0(parts[2], ":", parts[3], " (", parts[1], ")")
}, character(1))

dat <- dat %>% 
  dplyr::arrange(., as.numeric(lead_snp_pos)) %>%
  dplyr::arrange(., as.numeric(gsub("chr","",lead_snp_chr))) %>%
  dplyr::mutate(., Ord1=nrow(dat):1)

ggplot(dat,aes(x=reorder(gwashit_forplot, Ord1), y = celltype_forplots, color = abs(slope_forplotting), size = -log10(lead_H4_variant_pval)))+
  theme_classic() + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_point(shape = 1,colour = "grey")+
  facet_grid(.~cellgroup_forplots, scales = "free", space = "free")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  scale_color_gradient2(low = "blue",mid = "yellow",high = "darkgreen",midpoint = 0)+
  scale_size(range = c(3, 5))+xlab("GWAS hit") +ylab("Cell type")+ coord_flip()+ggtitle(paste0("Coloc with ",GWAS_ID," (PP.H4>0.8)(dots=lead PPH4 eQTLs>0.8)"))

ggsave(paste0(DIR_MAIN,"/0_plots/colocsummary_dotplot_",GWAS_ID,"_wo_bulklike_checksigeQTL.pdf"),height=15,width=11)
write.csv(dat,paste0(DIR_MAIN,"/1_csvfiles/coloc_sigresults_",GWAS_ID,"_checksigeQTL.csv"),row.names=F)

## dotplot with PPH4
dat <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_sigresults_",GWAS_ID,"_checksigeQTL.csv"))
dat$celltype_forplots <- gsub("_", " ", dat$cell_type)
dat$celltype_forplots <- factor(dat$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
dat$cellgroup_forplots <- gsub("_", " ", dat$cellgroup)
dat$cellgroup_forplots <- factor(dat$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))

library("viridis")
ggplot(dat,aes(x=reorder(gwashit_forplot, Ord1), y = celltype_forplots, color = PP.H4.abf, size = abs(lead_H4_variant_slope)))+
  theme_classic() + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_viridis(option = "D",direction = -1)+
  geom_point(shape = 1,colour = "grey")+
  facet_grid(.~cellgroup_forplots, scales = "free", space = "free")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  scale_size(range = c(2, 6))+xlab("GWAS hit") +ylab("Cell type")+ coord_flip()+ggtitle(paste0("Coloc with ",GWAS_ID))
ggsave(paste0(DIR_MAIN,"/0_plots/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL.pdf"),height=14,width=10)