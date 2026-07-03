##compare with bulk (all-pbmc is denoted bulklike in this script)
library(ggplot2)
library(scales)
library(ggrepel)
library(ggpubr)
library(Seurat)
library(dplyr)


## look at eQTLS
optimPC <- 25
bulkresults <- read.table(paste0("/path/eQTLresults/TensorQTL_eQTLS/dMean__All_all/",optimPC,"pcs/base_output/base/Cis_eqtls_qval.tsv"),fill=T,header=T)
all_genes_inbulk <- bulkresults$phenotype_id
sig_genes_inbulk <- bulkresults[bulkresults$qval < 0.05,]$phenotype_id

## compare with single cell eQTLs
scresults = read.csv("/path/eQTLresults/1_csvfiles/genestested.csv")
all_genes_insc <- unique(scresults$phenotype_id)
all_genes <- as.data.frame(all_genes_insc)
colnames(all_genes) <- "phenotype_id"

scresults_sig = read.csv("/path/eQTLresults/1_csvfiles/genessig.csv")
sig_genes_insc <- unique(scresults_sig$phenotype_id)

all_genes$sig_in_sc <- "NOTSiginsc"
all_genes$sig_in_sc[all_genes$phenotype_id %in% sig_genes_insc] <- "Siginsc"

all_genes$sig_in_bulk <- "NOTtestedinbulk"
all_genes$sig_in_bulk[all_genes$phenotype_id %in% all_genes_inbulk] <- "Testedinbulk"
all_genes$sig_in_bulk[all_genes$phenotype_id %in% sig_genes_inbulk] <- "Siginbulk"

all_genes$combined <- paste0(all_genes$sig_in_sc,"_",all_genes$sig_in_bulk)
table(all_genes$combined)

## score cell type specificity of gene expression to compare bulk like to single cell 
seurat_obj <- readRDS("/path/testrun_input_n281_seurat.RDS")

# 1. Set cell identities (or use a custom metadata column)
Idents(seurat_obj) <- "Celltype_level1" 
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
genes_to_tau <- all_genes[all_genes$combined %in% c("NOTSiginsc_Siginbulk","Siginsc_NOTtestedinbulk","Siginsc_Testedinbulk","Siginsc_Siginbulk"),]$phenotype_id

# 2. Calculate average expression for each cluster
avg_expr <- AverageExpression(seurat_obj, return.seurat = F, features = genes_to_tau,layer='data')
avg_expr_matrix <- as.data.frame(as.matrix(avg_expr[["RNA"]]))
avg_expr_matrix <- avg_expr_matrix[,colnames(avg_expr_matrix) %in% c("CD56Bright-NK-cells","CD56Dim-NK-cells","Classical-Monocytes","CM-CD4-T-cells","Cytotoxic-CD4-T-cells","EM-CD4-T-cells","Naive-CD4-T-cells","Regulatory-CD4-T-cells","CM-CD8-T-cells","EM-CD8-T-cells", "Naive-CD8-T-cells","TEMRA","DN-T-cells","Memory-B-cells","Naive-B-cells")]

# 3. Compute Tau score for each gene
compute_tau <- function(expression_vector) {
  if (all(expression_vector == 0)) return(NA)  # avoid division by zero
  x_max <- max(expression_vector)
  tau <- sum(1 - (expression_vector / x_max)) / (length(expression_vector) - 1)
  return(tau)
}
tau_scores <- apply(avg_expr_matrix, 1, compute_tau)
tau_df <- data.frame(gene = rownames(avg_expr_matrix),tau = tau_scores)

# get gene groups
genes_siginsc_orbulk <- all_genes[all_genes$combined %in% c("NOTSiginsc_Siginbulk","Siginsc_NOTtestedinbulk","Siginsc_Testedinbulk","Siginsc_Siginbulk"),]
genes_siginsc_orbulk$group3 <- genes_siginsc_orbulk$combined
genes_siginsc_orbulk$group3[genes_siginsc_orbulk$combined %in% c("Siginsc_NOTtestedinbulk","Siginsc_Testedinbulk")] <- "Siginsc_NOTSiginbulk"
genes_siginsc_orbulk$tau <- tau_df$tau[match(genes_siginsc_orbulk$phenotype_id,tau_df$gene)]

write.csv(genes_siginsc_orbulk,"/path/eQTLresults/1_csvfiles/genes_siginsc_orbulk.csv")

genes_siginsc_orbulk$group3 <- factor(genes_siginsc_orbulk$group3, levels=c("Siginsc_NOTSiginbulk","Siginsc_Siginbulk","NOTSiginsc_Siginbulk"))
ggplot(data=genes_siginsc_orbulk, aes(x=group3, y=tau,fill=group3)) +  geom_violin(alpha=0.8)+scale_fill_manual(values=c("#FAC881", "#F3A464", "#E87443"))+theme_classic()+xlab("")+ guides(fill="none")+geom_boxplot(width=0.08, outliers = F)+ scale_x_discrete(labels= c("Only at cell type-level","Detected in both","Only in all-PBMC"))
ggsave("/path/eQTLresults/0_plots/tauscores.pdf", width=5,height=3)

# tau score to average expression correlation
seurat_obj <- readRDS("/path/testrun_input_n281_seurat.RDS")
Idents(seurat_obj) <- "Celltype_level0" 
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- subset(seurat_obj, features = genes_siginsc_orbulk$phenotype_id)
norm_mat <- GetAssayData(seurat_obj, slot = "data")
avg_exp <- rowMeans(norm_mat)
avg_exp <- as.data.frame(avg_exp)
genes_siginsc_orbulk$avg_exp <- avg_exp$avg_exp[match(genes_siginsc_orbulk$phenotype_id,row.names(avg_exp))]

ggplot(data=genes_siginsc_orbulk, aes(x=group3, y=avg_exp,fill=group3))+theme_classic()+xlab("")+ guides(fill="none")+ylab("Average expression of all-PBMC")+geom_boxplot(width=0.8, outliers = F)+ scale_x_discrete(labels= c("Only at cell\ntype-level","Detected in both","Only in all-PBMC"))+scale_fill_manual(values=c("#FAC881", "#F3A464", "#E87443"))
ggsave("/path/eQTLresults/0_plots/average_expression.pdf", width=5,height=4)

kruskal.test(avg_exp ~ combined_forplot,data = genes_siginsc_orbulk)
kruskal.test(avg_exp ~ group3,data = genes_siginsc_orbulk)

res <- wilcox.test(avg_exp ~ group3, data=genes_siginsc_orbulk[genes_siginsc_orbulk$group3 %in% c("Siginsc_Siginbulk","NOTSiginsc_Siginbulk"),])
res$p.value
res <- wilcox.test(avg_exp ~ group3, data=genes_siginsc_orbulk[genes_siginsc_orbulk$group3 %in% c("Siginsc_Siginbulk","Siginsc_NOTSiginbulk"),])
res$p.value
res <- wilcox.test(avg_exp ~ group3, data=genes_siginsc_orbulk[genes_siginsc_orbulk$group3 %in% c("NOTSiginsc_Siginbulk","Siginsc_NOTSiginbulk"),])
res$p.value

res <- wilcox.test(tau ~ group3, data=genes_siginsc_orbulk[genes_siginsc_orbulk$group3 %in% c("Siginsc_Siginbulk","NOTSiginsc_Siginbulk"),])
res$p.value
res <- wilcox.test(tau ~ group3, data=genes_siginsc_orbulk[genes_siginsc_orbulk$group3 %in% c("Siginsc_Siginbulk","Siginsc_NOTSiginbulk"),])
res$p.value
res <- wilcox.test(tau ~ group3, data=genes_siginsc_orbulk[genes_siginsc_orbulk$group3 %in% c("NOTSiginsc_Siginbulk","Siginsc_NOTSiginbulk"),])
res$p.value

### slope by celltype corr between bulk and sc
optimPC <- 25
bulkresults <- read.table(paste0("/path/eQTLresults/TensorQTL_eQTLS/dMean__All_all/",optimPC,"pcs/base_output/base/Cis_eqtls_qval.tsv"),fill=T,header=T)

pearsonr <- data.frame(matrix(ncol=5))
colnames(pearsonr) <- c("celltype",'nvariant_allsig','nvariant_samevariant','pearson_r','pval')

for(i in c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","DN_T_cells","Memory_B_cells","Naive_B_cells")){
  qval <- read.table(paste0("/path/eQTLresults/",i,"/results/TensorQTL_eQTLS/dMean__",i,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"), header=T,fill=T)
  bulkresults$sc_qval <- qval$qval[match(bulkresults$phenotype_id,qval$phenotype_id)]
  bulkresults$sc_slope <- qval$slope[match(bulkresults$phenotype_id,qval$phenotype_id)]
  bulkresults$sc_af <- qval$af[match(bulkresults$phenotype_id,qval$phenotype_id)]
  bulkresults$sc_variant <- qval$variant_id[match(bulkresults$phenotype_id,qval$phenotype_id)]
  
  bulkresults_tmp <- bulkresults[is.na(bulkresults$sc_slope) == F,]
  bulkresults_tmp$variant_info <- "Not same variant"
  bulkresults_tmp$variant_info[bulkresults_tmp$variant_id == bulkresults_tmp$sc_variant] <- "Same variant"
  
  bulkresults_tmp$group <- "NotSig"
  bulkresults_tmp$group[bulkresults_tmp$sc_qval < 0.05] <- "Sig_in_only_sc"
  bulkresults_tmp$group[bulkresults_tmp$qval < 0.05] <- "Sig_in_only_bulk"
  bulkresults_tmp$group[(bulkresults_tmp$qval < 0.05)&(bulkresults_tmp$sc_qval < 0.05)] <- "Sig_in_both"
  bulkresults_tmp <- bulkresults_tmp[bulkresults_tmp$group == "Sig_in_both",]
  
  pearson_results <- cor.test(bulkresults_tmp[bulkresults_tmp$variant_info == "Same variant",]$slope,bulkresults_tmp[bulkresults_tmp$variant_info == "Same variant",]$sc_slope, method = 'pearson')
  pearsonr <- rbind(pearsonr,c(i,nrow(bulkresults_tmp),nrow(bulkresults_tmp[bulkresults_tmp$variant_info == "Same variant",]),pearson_results[["estimate"]][["cor"]],pearson_results[["p.value"]]))
}

pearsonr <- pearsonr[-1,]
pearsonr$nvariant_samevariant <- as.numeric(pearsonr$nvariant_samevariant)
pearsonr$nvariant_allsig <- as.numeric(pearsonr$nvariant_allsig)
pearsonr$perc_overlap_variants <-  pearsonr$nvariant_samevariant / pearsonr$nvariant_allsig * 100
pearsonr <- pearsonr[,c("celltype","perc_overlap_variants","pearson_r","pval")]

eGenesummary <- read.csv("/path/eQTLresults/1_csvfiles/eGenesummary.csv")
pearsonr$numCells_inc_eQTL <- eGenesummary$numCells_inc_eQTL[match(pearsonr$celltype,eGenesummary$celltype)]
pearsonr$pearson_r <- as.numeric(pearsonr$pearson_r)

pearsonr$cellgroup <- NA
pearsonr$cellgroup[pearsonr$celltype %in% c("CD56Bright_NK_cells","CD56Dim_NK_cells")] <- "NK"
pearsonr$cellgroup[pearsonr$celltype %in% c("Classical_Monocytes")] <- "Mono"
pearsonr$cellgroup[pearsonr$celltype %in% c("CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells")]<- "CD4_T"
pearsonr$cellgroup[pearsonr$celltype %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
pearsonr$cellgroup[pearsonr$celltype %in% c("DN_T_cells")] <- "Other_T"
pearsonr$cellgroup[pearsonr$celltype %in% c("Memory_B_cells","Naive_B_cells")] <- "B"
pearsonr$celltype_forplots <- gsub("_", " ", pearsonr$celltype)
pearsonr$celltype_forplots <- factor(pearsonr$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells"))
pearsonr$cellgroup_forplots <- gsub("_", " ", pearsonr$cellgroup)
pearsonr$cellgroup_forplots <- factor(pearsonr$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK"))

library(ggrepel)
library(scales)

### slope for only those with same eSNP
ggplot(pearsonr, aes(x=numCells_inc_eQTL, y=pearson_r,label=celltype_forplots,color=cellgroup_forplots))+ geom_point(size=1,alpha=0.8)+theme_classic()+xlab("Number of cells included in eQTL mapping")+ylab("Pearson's r between slope of bulk-like eQTL and sc-eQTL")+ggtitle("")+geom_text_repel(hjust=-0.1, vjust=0, max.overlaps = 1)+ scale_x_continuous(labels = label_comma()) + ylim(0.5,1)+scale_color_brewer(palette="Dark2")
#ggsave("/path/eQTLresults/0_plots/corr_numcells_pearson.pdf", width=6,height=4)

