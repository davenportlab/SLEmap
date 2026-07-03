
library(data.table)
library(rlist)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(forcats)

coloc <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/1_csvfiles/coloc_sigresults_GCST90270940_checksigeQTL.csv")

### calculate proportion of coloc hits over significant independent eQTLs by cell type
eGene_summary <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/eGenesummary_wBulklike.csv")
coloc_ct <- as.data.frame(table(coloc$cell_type))
eGene_summary$n_coloc <- coloc_ct$Freq[match(eGene_summary$celltype,coloc_ct$Var1)]
eGene_summary$propcoloc <- eGene_summary$n_coloc / eGene_summary$numGenes_conditionalsig 
eGene_summary$propcoloc_eGene <- eGene_summary$n_coloc / eGene_summary$numGenes_sig

otar <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/coloc_output_otar_26.03.csv")
otar_ct <- as.data.frame(table(otar$cell_type,otar$otar_overlap))
eGene_summary$novel_coloc <- otar_ct[otar_ct$Var2 == T,]$Freq[match(eGene_summary$celltype,otar_ct$Var1)]
eGene_summary$propnovelcoloc <- eGene_summary$novel_coloc / eGene_summary$n_coloc
plot(eGene_summary$numGenes_sig,eGene_summary$n_coloc)

conditional_eQTL <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/conditionaleQTL_simple.csv")
table(conditional_eQTL$celltype)
length(unique(conditional_eQTL$phenotype_id))

#get GWAS
GWAS_input= "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/14_colocalization/inputs/gwas/"
GWAS_ID="GCST90270940"
GWAS_INPUT <- fread (paste0(GWAS_input,GWAS_ID,"/",GWAS_ID,"_for_eqtl.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::rename(., variant_id_v2=variant_id) %>%
  dplyr::mutate(variant_id = paste0(chr, "_", pos))
GWAS_INPUT <- GWAS_INPUT[!((GWAS_INPUT$chr == "chr6") & (GWAS_INPUT$pos > 25000000)&(GWAS_INPUT$pos < 34000000)),] 
GWAS_INPUT$pos <- as.numeric(GWAS_INPUT$pos)
GWAS_INPUT$p_value <- as.numeric(GWAS_INPUT$p_value)

####NFKB1 in all vs cmcd4t####
coloc_use <- coloc[coloc$gene_symbol == "NFKB1",]
coloc_use
gene_id =coloc_use$gene_id[1]
pos_start = coloc_use$locusStart[1]
pos_end =coloc_use$locusEnd[1]
gene_chr=gsub("chr", "",coloc_use$lead_snp_chr[1])

GWAS_input_plot <- GWAS_INPUT[GWAS_INPUT$chr == paste0("chr",gene_chr),]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos < pos_end,]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos > pos_start,]

##allcells - NFKB1
EQTL_input <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/"
nominal <- fread(paste0("grep ",gene_id," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
nominal$coloclead <- "NotColocLead"
nominal$celltype <- "All-PBMC"

##cm cd4 t cells- NFKB1
celltype <- "CM_CD4_T_cells"
EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
nominal2 <- fread(paste0("grep ",gene_id," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
colnames(nominal2) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
nominal2$POS <- as.numeric(sapply(strsplit(nominal2$variant_id, "_"), `[`, 2))
nominal2$chr_pos <- sapply(strsplit(nominal2$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
nominal2$coloclead <- "NotColocLead"
nominal2$coloclead[nominal2$chr_pos == coloc_use[coloc_use$cell_type == celltype,]$lead_H4_variant] <- "ColocLead"
nominal2$celltype <- celltype

nominal <- rbind(nominal,nominal2)
p3 <- ggplot(nominal[(nominal$chr_pos %in%overlapping_variants),], aes(x = POS, y = -log10(pval_nominal),color=celltype,shape=coloclead)) +geom_point(size = 1.5, alpha = 0.7)+scale_color_manual(values=c("#E76F51", "#2A9D8F")) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("NFKB1"))+scale_shape_manual(values=c("ColocLead"=17,"NotColocLead"=16))+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10))

##MANBA
coloc_use <- coloc[coloc$gene_symbol == "MANBA",]
coloc_use
gene_id =coloc_use$gene_id[1]
##allcells - MANBA
nominal <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/2_indep_coloc/All/nominal_p/ENSG00000109323_regress_chr4_102759827_chr4_102713063.csv",header=T,row.names=1)

nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
nominal$coloclead <- "NotColocLead"
nominal$coloclead[nominal$chr_pos == coloc_use[coloc_use$cell_type == celltype,]$lead_H4_variant] <- "ColocLead"
nominal$celltype <- "All-PBMC"

##cm cd4 t cells- MANBA
celltype <- "CM_CD4_T_cells"
nominal2 <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/2_indep_coloc/CM_CD4_T_cells/nominal_p/ENSG00000109323_regress_chr4_102735220.csv",header=T,row.names=1)
nominal2$POS <- as.numeric(sapply(strsplit(nominal2$variant_id, "_"), `[`, 2))
nominal2$chr_pos <- sapply(strsplit(nominal2$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
nominal2$coloclead <- "NotColocLead"
nominal2$coloclead[nominal2$chr_pos == coloc_use[coloc_use$cell_type == celltype,]$lead_H4_variant] <- "ColocLead"
nominal2$celltype <- celltype
nominal <- rbind(nominal,nominal2)
p4 <- ggplot(nominal[(nominal$chr_pos %in%overlapping_variants),], aes(x = POS, y = -log10(pval_nominal),color=celltype,shape=coloclead)) +geom_point(size = 1.5, alpha = 0.7)+scale_color_manual(values=c("#E76F51", "#2A9D8F")) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("MANBA"))+scale_shape_manual(values=c("ColocLead"=17,"NotColocLead"=16))+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10))

p1 <- ggplot(GWAS_input_plot[(GWAS_input_plot$variant_id %in% overlapping_variants),], aes(x = pos, y = -log10(p_value))) + geom_point(size = 1.3, alpha = 0.7,color="grey50") +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end)) + geom_hline(yintercept=5, linetype="dashed")+theme_classic()+ggtitle("GWAS")+scale_y_continuous(labels = function(x) formatC(x, width = 10))
p2 <- gg_genetracks(locus(xrange=c(pos_start,pos_end), seqname=paste0("chr",gene_chr),ens_db = "EnsDb.Hsapiens.v86"), filter_gene_biotype = 'protein_coding', gene_col = "#90A4AE",exon_col = "#90A4AE",exon_border = "#90A4AE")

ggarrange(p1,p3,p4,p2,ncol=1,common.legend=T)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/forpaper/NFKB1_MANBA.pdf",width=6,height=9)


#### FHL3 and C1orf122 in TEMRA and classical monocytes####
gwas_hit_use <- "chr1_37881745"
coloc_use <- coloc[coloc$gwas_hit == gwas_hit_use,]
GWAShit_chr = sapply(strsplit(gwas_hit_use, "_"), `[`, 1)
GWAShit_pos = as.numeric(sapply(strsplit(gwas_hit_use, "_"), `[`, 2))
pos_start = coloc_use$locusStart[coloc_use$gwas_hit == gwas_hit_use][1]
pos_end =coloc_use$locusEnd[coloc_use$gwas_hit == gwas_hit_use][1]

gene_symbol1="FHL3"
gene_symbol2="C1orf122"
#gene_symbol3="INPP5B"
gene1=coloc_use$gene_id[coloc_use$gene_symbol==gene_symbol1][1]
gene2=coloc_use$gene_id[coloc_use$gene_symbol==gene_symbol2][1]
#gene3=coloc_use$gene_id[coloc_use$gene_symbol==gene_symbol3][1]
gene_chr <- gsub("chr", "",coloc_use$lead_snp_chr[1])

GWAS_input_plot <- GWAS_INPUT[GWAS_INPUT$chr == GWAShit_chr,]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos < pos_end,]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos > pos_start,]


#### start with clssical monocytes
plot <- list()
for(celltype in c("TEMRA","Classical_Monocytes")){
  if(celltype == "All"){
    EQTL_input <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/"
  }else{
    EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
  }
  
  ## gene1
  nominal <- fread(paste0("grep ",gene1," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
  nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
  nominal$gene <- gene_symbol1
  nominal$coloclead <- "NotColocLead"
  nominal$coloclead[nominal$chr_pos == coloc_use[coloc_use$cell_type == celltype & coloc_use$gene_id == gene1,]$lead_H4_variant] <- "ColocLead"
  overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
  
  
  ## gene2
  nominal2 <- fread(paste0("grep ",gene2," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  colnames(nominal2) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  nominal2$POS <- as.numeric(sapply(strsplit(nominal2$variant_id, "_"), `[`, 2))
  nominal2$chr_pos <- sapply(strsplit(nominal2$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
  nominal2$gene <- gene_symbol2
  nominal2$coloclead <- "NotColocLead"
  nominal2$coloclead[nominal2$chr_pos == coloc_use[coloc_use$cell_type == celltype & coloc_use$gene_id == gene2,]$lead_H4_variant] <- "ColocLead"
  
  
  
  # ## gene3
  # nominal3 <- fread(paste0("grep ",gene3," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
  # colnames(nominal3) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
  # nominal3$POS <- as.numeric(sapply(strsplit(nominal3$variant_id, "_"), `[`, 2))
  # nominal3$chr_pos <- sapply(strsplit(nominal3$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
  # nominal3$gene <- gene_symbol3
  # nominal3$coloclead <- "NotColocLead"
  # nominal3$coloclead[nominal3$chr_pos == coloc_use[coloc_use$cell_type == celltype & coloc_use$gene_id == gene3,]$lead_H4_variant] <- "ColocLead"
  
  #pvalthres <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/1_csvfiles/nominal_pval_thres.csv")
  #pvalthres1 <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene1)]
  #pvalthres2 <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene2)]
  #ifelse(length(pvalthres1)==0,pvalthres<- pvalthres2,pvalthres<- pvalthres1)
  
  ## bind both genes
  nominal <- rbind(nominal,nominal2)
  #nominal <- rbind(nominal,nominal2,nominal3)
  plot <- list.append(plot,ggplot(nominal[(nominal$chr_pos %in%overlapping_variants),], aes(x = POS, y = -log10(pval_nominal),color=gene,shape=coloclead)) +geom_point(size = 1.5, alpha = 0.7)+scale_color_manual(values=c("#8E5BD9","#1F7A5A")) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype))+scale_shape_manual(values=c(17,16))+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10)))

}

p1 <- ggplot(GWAS_input_plot[(GWAS_input_plot$variant_id %in% overlapping_variants),], aes(x = pos, y = -log10(p_value))) + geom_point(size = 1.3, alpha = 0.7,color="grey50") +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end)) + geom_hline(yintercept=5, linetype="dashed")+theme_classic()+ggtitle("GWAS")+scale_y_continuous(labels = function(x) formatC(x, width = 10))

p2 <- gg_genetracks(locus(xrange=c(pos_start,pos_end), seqname=paste0("chr",gene_chr),ens_db = "EnsDb.Hsapiens.v86"), filter_gene_biotype = 'protein_coding', gene_col = "#90A4AE",exon_col = "#90A4AE",exon_border = "#90A4AE")+scale_y_continuous(labels = function(x) formatC(x, width = 10))

ggarrange(p1,plot[[1]],plot[[2]],p2, ncol=1,nrow=4,common.legend=T)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/forpaper/C1orf122_FHL3.pdf",width=5,height=8)

####plot traf1 - color by celltype
####NFKB1 in all vs cmcd4t####
coloc_use <- coloc[coloc$gene_symbol == "TRAF1",]
coloc_use
gene_id =coloc_use$gene_id[1]
pos_start = coloc_use$locusStart[1]
pos_end =coloc_use$locusEnd[1]
gene_chr=gsub("chr", "",coloc_use$lead_snp_chr[1])

GWAS_input_plot <- GWAS_INPUT[GWAS_INPUT$chr == paste0("chr",gene_chr),]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos < pos_end,]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos > pos_start,]

##allcells - TRAF1
EQTL_input <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/"
nominal <- fread(paste0("grep ",gene_id," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
nominal$coloclead <- "NotColocLead"
nominal$celltype <- "All-PBMC"

##cm cd4 t cells- TRAF1
celltype <- "CM_CD4_T_cells"
EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
nominal2 <- fread(paste0("grep ",gene_id," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
colnames(nominal2) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
nominal2$POS <- as.numeric(sapply(strsplit(nominal2$variant_id, "_"), `[`, 2))
nominal2$chr_pos <- sapply(strsplit(nominal2$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
nominal2$coloclead <- "NotColocLead"
nominal2$coloclead[nominal2$chr_pos == coloc_use[coloc_use$cell_type == celltype,]$lead_H4_variant] <- "ColocLead"
nominal2$celltype <- celltype

##memoryB cells- TRAF1
celltype <- "Memory_B_cells"
EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
nominal3 <- fread(paste0("grep ",gene_id," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
colnames(nominal3) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
nominal3$POS <- as.numeric(sapply(strsplit(nominal3$variant_id, "_"), `[`, 2))
nominal3$chr_pos <- sapply(strsplit(nominal3$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
nominal3$coloclead <- "NotColocLead"
nominal3$celltype <- celltype

nominal <- rbind(nominal,nominal2,nominal3)
p3 <- ggplot(nominal[(nominal$chr_pos %in%overlapping_variants),], aes(x = POS, y = -log10(pval_nominal),color=celltype,shape=coloclead)) +geom_point(size = 1.5, alpha = 0.7)+scale_color_manual(values=c("#E76F51", "#2A9D8F","#6F42C1")) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("TRAF1"))+scale_shape_manual(values=c("ColocLead"=17,"NotColocLead"=16))+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10))
p1 <- ggplot(GWAS_input_plot[(GWAS_input_plot$variant_id %in% overlapping_variants),], aes(x = pos, y = -log10(p_value))) + geom_point(size = 1.3, alpha = 0.7,color="grey50") +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end)) + geom_hline(yintercept=5, linetype="dashed")+theme_classic()+ggtitle("GWAS")+scale_y_continuous(labels = function(x) formatC(x, width = 10))
p2 <- gg_genetracks(locus(xrange=c(pos_start,pos_end), seqname=paste0("chr",gene_chr),ens_db = "EnsDb.Hsapiens.v86"), filter_gene_biotype = 'protein_coding', gene_col = "#90A4AE",exon_col = "#90A4AE",exon_border = "#90A4AE",highlight="TRAF1",highlight_col="#FBB040")
ggarrange(p1,p3,p2,ncol=1,common.legend=T)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/forpaper/TRAF1_combined.pdf",width=6,height=7)

#### plot TRAF1 in CM CD4 T, memory B - color LD

coloc_use <- coloc[coloc$gene_symbol == "TRAF1",]
coloc_use

## FOR GWAS
GWAS_input= "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/14_colocalization/inputs/gwas/"

chr = "chr9"
GWAShit = 120927961
pos_start = 120720914
pos_end =121039712
gene="ENSG00000056558"
colocvariant = "chr9_120927961"
gene_chr <- sub("chr", "", strsplit(colocvariant, "_")[[1]][1])
coloc_lead_pos <- strsplit(colocvariant, "_")[[1]][2]
REF = "G"
ALT = "A"
coloc_genesymbol = "TRAF1"

GWAS_INPUT <- fread(paste0(GWAS_input,GWAS_ID,"/",GWAS_ID,"_for_eqtl.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::rename(., variant_id_v2=variant_id) %>%
  dplyr::mutate(variant_id = paste0(chr, "_", pos))
GWAS_INPUT <- GWAS_INPUT[!((GWAS_INPUT$chr == "chr6") & (GWAS_INPUT$pos > 25000000)&(GWAS_INPUT$pos < 34000000)),] 

GWAS_INPUT$pos <- as.numeric(GWAS_INPUT$pos)
GWAS_INPUT$p_value <- as.numeric(GWAS_INPUT$p_value)

GWAS_input_plot <- GWAS_INPUT[GWAS_INPUT$chr == chr,]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos < pos_end,]
GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos > pos_start,]
ggplot(GWAS_input_plot, aes(x = pos, y = -log10(p_value))) +
  geom_point(size = 1.5, alpha = 0.7) +scale_x_continuous(labels = scales::comma)

GWAS_input_plot$color <- "nothit"
GWAS_input_plot$color[GWAS_input_plot$pos == GWAShit] <- "GWAShit"

##get LD 
system2("bash", args = c("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/20_ancestry_coloc/6.5_LDforall.sh",paste0(colocvariant,"_",REF,"_",ALT),gene_chr, pos_start,pos_end))
ld <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/0_calcLD/all/",colocvariant,"_",REF,"_",ALT,".vcor"))

## FOR EQTL

getplot <- function(celltype){
  if(celltype == "All"){
    EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
  }else{
    EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
  }
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
  
  nominal$LD <- ld$V7[match(nominal$variant_id,ld$V6)]
  nominal$LD <- fct_rev(cut(nominal$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                            labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
  nominal$LD[nominal$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"
  
  GWAS_input_plot$overlap <- "not_overlap"
  GWAS_input_plot$overlap[GWAS_input_plot$variant_id %in% overlapping_variants] <- "overlap"
  GWAS_input_plot$coloc_leadsnp <- "notlead_colocsnp"
  GWAS_input_plot$coloc_leadsnp[GWAS_input_plot$variant_id == colocvariant] <- "lead_colocsnp"
  
  ####get nominal p-val thres
  #####double check that the gene plotted has a sig eQTL, so the pvalthres is reliable
  if(celltype == "All"){
    pvalthres <- read.table("//lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv",header=T,fill=T)
    pvalthres <- pvalthres[pvalthres$qval < 0.05,]
    pvalthres <- pvalthres$pval_nominal_threshold[pvalthres$phenotype_id == gene]
  }else{
    pvalthres <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/conditionaleQTL.csv")
    pvalthres <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
    pvalthres <- pvalthres[1]
  }
  
  p1 <- ggplot(GWAS_input_plot[GWAS_input_plot$overlap == "overlap",], aes(x = pos, y = -log10(p_value),shape=coloc_leadsnp)) + geom_point(aes(size=overlap,alpha = overlap),color="grey50")+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw()+scale_shape_manual(values=c(17,16))+ geom_hline(yintercept=5, linetype="dashed")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("")+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+scale_y_continuous(labels = function(x) formatC(x, width = 10))
  if(length(pvalthres) == 0){
    p2 <- ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("")+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+scale_y_continuous(limits = c(0, 5.5),labels = function(x) formatC(x, width = 10))
  }else{
    p2 <- ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept=-log10(pvalthres), linetype="dashed")+xlab("")+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+scale_y_continuous(limits = c(0, 5.5),labels = function(x) formatC(x, width = 10))
  }
  
  
  return(list(p1,p2))
  
}
CM_CD4_T_cells <- getplot("CM_CD4_T_cells")
Memory_B_cells <- getplot("Memory_B_cells")
All <- getplot("All")

p2 <- gg_genetracks(locus(xrange=c(pos_start,pos_end), seqname=chr,ens_db = "EnsDb.Hsapiens.v86"), filter_gene_biotype = 'protein_coding', gene_col = "#90A4AE",exon_col = "#90A4AE",exon_border = "#90A4AE",highlight="TRAF1",highlight_col="#FBB040")


ggarrange(CM_CD4_T_cells[[1]],CM_CD4_T_cells[[2]],Memory_B_cells[[2]],All[[2]],p2,ncol=1,common.legend = T)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/0_plots/forpaper/TRAF1.pdf",width=5,height=10)
