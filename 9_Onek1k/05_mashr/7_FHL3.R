###FHL3 plots

###Regional plot with GWAS, SLEmap, ONek1k

library(ggplot2)
library(ggpubr)
library(rlist)
library(data.table)
library(dplyr)
library(forcats)
library(scales)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)

## FOR GWAS
GWAS_input= "/path/colocalization/inputs/gwas/"
GWAS_ID="GCST90270940"

##GET coloc results
coloc <- read.csv(paste0("/path/onek1k_locus_breaker_coloc/","1_csvfiles/SLEmap_Onek1kgroup.csv"),row.names=1)

GWAS_INPUT <- fread (paste0(GWAS_input,GWAS_ID,"/",GWAS_ID,"_for_eqtl.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::rename(., variant_id_v2=variant_id) %>%
  dplyr::mutate(variant_id = paste0(chr, "_", pos))
GWAS_INPUT <- GWAS_INPUT[!((GWAS_INPUT$chr == "chr6") & (GWAS_INPUT$pos > 25000000)&(GWAS_INPUT$pos < 34000000)),] 

GWAS_INPUT$pos <- as.numeric(GWAS_INPUT$pos)
GWAS_INPUT$p_value <- as.numeric(GWAS_INPUT$p_value)
locusbreaker <- read.csv("/path/locus_breaker_coloc/gwashit_locusbreaker_forinput.csv")


###for specific ones
row.names(coloc[coloc$cell_type == "TEMRA" & coloc$gene_symbol == "FHL3",])
i <- 8

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
system2("bash", args = c("/path/ancestry_coloc/6.5_LDforall.sh",paste0(colocvariant,"_",REF,"_",ALT),gene_chr, pos_start,pos_end))
all_ld <- read.table(paste0("/path/ancestry_coloc/0_calcLD/all/",colocvariant,"_",REF,"_",ALT,".vcor"))

EQTL_input=paste0("/path/eQTLresults/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")

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

pvalthres <- read.csv("/path/eQTLresults/1_csvfiles/conditionaleQTL.csv")
pvalthres <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
pvalthres <- pvalthres[1]

gwas_plot <- ggplot(GWAS_input_plot[GWAS_input_plot$overlap == "overlap",], aes(x = pos, y = -log10(p_value),shape=coloc_leadsnp)) + geom_point(size = 1.3, alpha = 0.7,color="grey50")+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw()+scale_shape_manual(values=c(17,16))+ geom_hline(yintercept=5, linetype="dashed")+scale_y_continuous(labels = function(x) formatC(x, width = 10))+theme_classic()# to show all snps in region while having alpha higher in overlapping ones

slemap_FHL3_plot <- ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(size = 1.3, alpha = 0.7)+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ geom_hline(yintercept=-log10(pvalthres), linetype="dashed")+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10)) # to show all snps in region while having alpha higher in overlapping ones

slemap <- nominal
slemap$dataset <- "slemap"

roof <- max(-log10(nominal$pval_nominal))
#############for onek1k #########
##get LD 
system2("bash", args = c("/path/Onek1k/9.5_LD.sh",paste0(colocvariant,"_",REF,"_",ALT),gene_chr, pos_start,pos_end))

if(file.exists(paste0("/path/Onek1k/Genotypes/calcLD/",colocvariant,"_",REF,"_",ALT,".vcor"))){
  all_ld <- read.table(paste0("/path/Onek1k/Genotypes/calcLD/",colocvariant,"_",REF,"_",ALT,".vcor"))
  
}else{
  all_ld <- read.table(paste0("/path/Onek1k/Genotypes/calcLD/",colocvariant,"_",REF,"_",ALT,".snplist"))
  all_ld$LD <- NA
  colnames(all_ld) <- c("V6","V7")
}

EQTL_input=paste0("/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")

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
nominal$LD <- as.numeric(nominal$LD)
nominal$LD <- fct_rev(cut(nominal$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                          labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
nominal$LD[nominal$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"

GWAS_input_plot$overlap <- "not_overlap"
GWAS_input_plot$overlap[GWAS_input_plot$variant_id %in% overlapping_variants] <- "overlap"
GWAS_input_plot$coloc_leadsnp <- "notlead_colocsnp"
GWAS_input_plot$coloc_leadsnp[GWAS_input_plot$variant_id == colocvariant] <- "lead_colocsnp"

pvalthres <- read.csv("/path/Onek1k/eQTLmapping/allsig_genes_onek1k.csv",header=T)
pvalthres <- pvalthres$pval_nominal_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
pvalthres <- pvalthres[1]

roof_onek1k <- max(-log10(nominal$pval_nominal))
ifelse(roof_onek1k>roof,roof <- roof_onek1k, roof <- roof)

onek1k_FHL3_plot <- ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(size = 1.3, alpha = 0.7)+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("Onek1k-",celltype," ",coloc_genesymbol))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ geom_hline(yintercept=-log10(pvalthres), linetype="dashed")+ylim(0,roof)+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10)) # to show all snps in region while having alpha higher in overlapping ones

onek1k <- nominal
onek1k$dataset <- "onek1k"

together_FHL3 <- rbind(slemap,onek1k)
together_FHL3_plot <- ggplot(together_FHL3[together_FHL3$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=dataset,shape=coloc_leadsnp))+geom_point(size = 1.3, alpha = 0.7) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol))+theme_bw()+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ylim(0,roof)+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10))+scale_color_manual(values = c("onek1k" = "#4C78A8","slemap" = "#FBB040"))
together_FHL3_plot

###for specific ones
row.names(coloc[coloc$cell_type == "Classical_Monocytes" & coloc$gene_symbol == "C1orf122",])
i <- 4

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
system2("bash", args = c("/path/ancestry_coloc/6.5_LDforall.sh",paste0(colocvariant,"_",REF,"_",ALT),gene_chr, pos_start,pos_end))
all_ld <- read.table(paste0("/path/ancestry_coloc/0_calcLD/all/",colocvariant,"_",REF,"_",ALT,".vcor"))

EQTL_input=paste0("/path/eQTLmapping/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")

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

pvalthres <- read.csv("/path/eQTLmapping/1_csvfiles/conditionaleQTL.csv")
pvalthres <- pvalthres$nominal_pvalue_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
pvalthres <- pvalthres[1]

slemap_C1orf122_plot <- ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(size = 1.3, alpha = 0.7)+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ geom_hline(yintercept=-log10(pvalthres), linetype="dashed")+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10)) # to show all snps in region while having alpha higher in overlapping ones

slemap <- nominal
slemap$dataset <- "slemap"

roof <- max(-log10(nominal$pval_nominal))
#############for onek1k #########
##get LD
system2("bash", args = c("/path/Onek1k/9.5_LD.sh",paste0(colocvariant,"_",REF,"_",ALT),gene_chr, pos_start,pos_end))

if(file.exists(paste0("/path/Onek1k/Genotypes/calcLD/",colocvariant,"_",REF,"_",ALT,".vcor"))){
  all_ld <- read.table(paste0("/path/Onek1k/Genotypes/calcLD/",colocvariant,"_",REF,"_",ALT,".vcor"))
  
}else{
  all_ld <- read.table(paste0("/path/Onek1k/Genotypes/calcLD/",colocvariant,"_",REF,"_",ALT,".snplist"))
  all_ld$LD <- NA
  colnames(all_ld) <- c("V6","V7")
}

EQTL_input=paste0("/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")

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
nominal$LD <- as.numeric(nominal$LD)
nominal$LD <- fct_rev(cut(nominal$LD,breaks = seq(0, 1, by = 0.2),include.lowest = TRUE,right = FALSE,
                          labels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")))
nominal$LD[nominal$coloc_leadsnp=="lead_colocsnp"] <- "0.8–1.0"

GWAS_input_plot$overlap <- "not_overlap"
GWAS_input_plot$overlap[GWAS_input_plot$variant_id %in% overlapping_variants] <- "overlap"
GWAS_input_plot$coloc_leadsnp <- "notlead_colocsnp"
GWAS_input_plot$coloc_leadsnp[GWAS_input_plot$variant_id == colocvariant] <- "lead_colocsnp"

pvalthres <- read.csv("/path/Onek1k/eQTLmapping/allsig_genes_onek1k.csv",header=T)
pvalthres <- pvalthres$pval_nominal_threshold[paste0(pvalthres$celltype,pvalthres$phenotype_id) == paste0(celltype,gene)]
pvalthres <- pvalthres[1]

roof_onek1k <- max(-log10(nominal$pval_nominal))
ifelse(roof_onek1k>roof,roof <- roof_onek1k, roof <- roof)

onek1k_C1orf122_plot <- ggplot(nominal[nominal$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=LD,shape=coloc_leadsnp)) +geom_point(size = 1.3, alpha = 0.7)+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0("Onek1k-",celltype," ",coloc_genesymbol))+theme_bw()+scale_color_manual(values=c("0.8–1.0"="#D33F49", "0.6–0.8"="#F79256", "0.4–0.6"="#90BE6D", "0.2–0.4"="#4E9CFF", "0–0.2"="#00117F"))+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ geom_hline(yintercept=-log10(pvalthres), linetype="dashed")+ylim(0,roof)+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10)) # to show all snps in region while having alpha higher in overlapping ones

onek1k <- nominal
onek1k$dataset <- "onek1k"

together_C1orf122 <- rbind(slemap,onek1k)
together_C1orf122_plot <- ggplot(together_C1orf122[together_C1orf122$overlap == "overlap",], aes(x = POS, y = -log10(pval_nominal),color=dataset,shape=coloc_leadsnp))+geom_point(size = 1.3, alpha = 0.7) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol))+theme_bw()+scale_shape_manual(values=c("lead_colocsnp"=17,"notlead_colocsnp"=16))+ylim(0,roof)+theme_classic()+scale_y_continuous(labels = function(x) formatC(x, width = 10))+scale_color_manual(values = c("onek1k" = "#4C78A8","slemap" = "#FBB040"))
together_C1orf122_plot

genetrack <- gg_genetracks(locus(xrange=c(pos_start,pos_end), seqname=paste0("chr",gene_chr),ens_db = "EnsDb.Hsapiens.v86"), filter_gene_biotype = 'protein_coding', gene_col = "#90A4AE",exon_col = "#90A4AE",exon_border = "#90A4AE",highlight=c("C1orf122","FHL3"),highlight_col="#FBB040")

ggarrange(gwas_plot,together_C1orf122_plot,together_FHL3_plot,genetrack,ncol=1)
ggsave("/path/onek1k_locus_breaker_coloc/0_plots/forpaper/C1orf122_onek1k_slemap.pdf",width=7,height=8.5)
ggarrange(gwas_plot,gwas_plot,slemap_C1orf122_plot,slemap_FHL3_plot,onek1k_C1orf122_plot,onek1k_FHL3_plot,genetrack,genetrack,ncol=2,nrow=4,common.legend=T)
ggsave("/path/onek1k_locus_breaker_coloc/0_plots/forpaper/C1orf122_onek1k_slemap_suppl.pdf",width=10.5,height=8.5)
