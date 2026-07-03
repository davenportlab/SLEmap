## XGR on bulklike vs singlecell

library(XGR)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)

#RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
RData.location <- "/software/team282/XGR/bigdata"
options(stringsAsFactors = FALSE)
info <- xRDataLoader('EpigenomeAtlas_15Segments_info',RData.location=RData.location)
GR.annotations <- paste0('EpigenomeAtlas_15Segments_', names(info))
names(GR.annotations) <- info
GR.annotations <- GR.annotations[c(29:30, 32, 34, 37:48,61)]
GR.annotations

# Get info on conditional eQTLs: background should be all tested
resultsdir <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/"
conditionalsig <- read.csv(paste0(resultsdir,"1_csvfiles/conditionaleQTL.csv"))
conditionalsig$pair <- paste0(conditionalsig$phenotype_id,"_", conditionalsig$variant_id)
conditionalsig$SNP <- sub("_[ACGT]+_[ACGT]+$", "", conditionalsig$variant_id)
conditionalsig$sc_or_bulk <- "singlecell"
colnames(conditionalsig)

# get info on bulklike eQTLs
bulklike <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",header=T)
bulklike$celltype <- "All"
bulklike$sc_or_bulk <- "All"
bulklike$pair <- paste0(bulklike$phenotype_id,"_", bulklike$variant_id)
bulklike$SNP <- sub("_[ACGT]+_[ACGT]+$", "", bulklike$variant_id)
colnames(bulklike)

conditionalsig <- rbind(conditionalsig[colnames(conditionalsig) %in% colnames(bulklike)],bulklike[colnames(bulklike) %in% colnames(conditionalsig)])

### change from hg38 to hg19
conditionalsig <- conditionalsig %>%
  separate(SNP, into = c("chr", "pos"), sep = "_", remove=F)
conditionalsig$pos <- as.numeric(conditionalsig$pos)

gr <- GRanges(seqnames = conditionalsig$chr,ranges = IRanges(start = conditionalsig$pos, end = conditionalsig$pos))

chain <- import.chain("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/3_otheresources/hg38ToHg19.over.chain")

lifted <- liftOver(gr, chain)
success_idx <- lengths(lifted) == 1  # Logical vector, TRUE for successful rows
conditionalsig_lift <- conditionalsig[success_idx, ]
lifted_flat <- unlist(lifted)
conditionalsig_lift$chr_hg19 <- as.character(seqnames(lifted_flat))
conditionalsig_lift$pos_hg19 <- start(lifted_flat)
conditionalsig_lift$SNP <- paste0(conditionalsig_lift$chr_hg19,"_",conditionalsig_lift$pos_hg19)
print(paste0("Out of ",nrow(conditionalsig),", ",nrow(conditionalsig_lift)," lifted to hg19 successively"))

###by matched cell type
RData.location = "/software/team282/XGR/bigdata"

getfunctionenrichment <- function(background,interested,annotations){
  background <- conditionalsig_lift[conditionalsig_lift$celltype == background,]
  interested <- conditionalsig_lift[conditionalsig_lift$celltype == interested,]
  b19pos.background <- xSNPlocations(background$SNP, GR.SNP = "dbSNP_Common",RData.location = RData.location)
  b19pos.intstested <- xSNPlocations(interested$SNP, GR.SNP = "dbSNP_Common",RData.location = RData.location)
  GR.annotations_specific <- GR.annotations[annotations]
  ls_df <- lapply(1:length(GR.annotations_specific), function(i){
    GR.annotation_use <- GR.annotations_specific[i]
    message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation_use),
                    as.character(Sys.time())), appendLF=T)
    df <- xGRviaGenomicAnno(data.file=b19pos.intstested, 
                            background.file = b19pos.background,
                            format.file="GRanges",
                            p.tail="one-tail",
                            GR.annotation=GR.annotation_use, 
                            RData.location=RData.location, verbose=F)
    df$group <- names(GR.annotation_use)
    return(df)
  })
  df <- do.call(rbind, ls_df)
  df$fc[which(df$adjp > 0.05)] <- NA
  for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
  df$name <- gsub(" from peripheral blood", "", df$name)
  return(df) 
}

getfunctionenrichment <- function(background,interested,annotations){
  background <- conditionalsig_lift[conditionalsig_lift$celltype == background,]
  interested <- conditionalsig_lift[conditionalsig_lift$celltype == interested,]
  b19pos.background <- xSNPlocations(background$SNP, GR.SNP = "dbSNP_Common",RData.location = RData.location)
  b19pos.intstested <- xSNPlocations(interested$SNP, GR.SNP = "dbSNP_Common",RData.location = RData.location)
  GR.annotations_specific <- GR.annotations[annotations]
  ls_df <- lapply(1:length(GR.annotations_specific), function(i){
    GR.annotation_use <- GR.annotations_specific[i]
    message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation_use),
                    as.character(Sys.time())), appendLF=T)
    df <- xGRviaGenomicAnno(data.file=b19pos.intstested, 
                            background.file = b19pos.background,
                            format.file="GRanges",
                            p.tail="one-tail",
                            GR.annotation=GR.annotation_use, 
                            RData.location=RData.location, verbose=F)
    df$group <- names(GR.annotation_use)
    return(df)
  })
  df <- do.call(rbind, ls_df)
  df$fc[which(df$adjp > 0.05)] <- NA
  for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
  df$name <- gsub(" from peripheral blood", "", df$name)
  return(df) 
}


df_list <- list()
celltype <- "CD56Bright_NK_cells"
annotations <- c("E046 (Primary Natural Killer cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "CD56Dim_NK_cells"
annotations <- c("E046 (Primary Natural Killer cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "Classical_Monocytes"
annotations <- c("E029 (Primary monocytes from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "CM_CD4_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E037 (Primary T helper memory cells from peripheral blood 2)","E040 (Primary T helper memory cells from peripheral blood 1)","E043 (Primary T helper cells from peripheral blood)","E045 (Primary T cells effector/memory enriched from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "CM_CD8_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E045 (Primary T cells effector/memory enriched from peripheral blood)","E048 (Primary T killer memory cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "Cytotoxic_CD4_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E037 (Primary T helper memory cells from peripheral blood 2)","E040 (Primary T helper memory cells from peripheral blood 1)","E043 (Primary T helper cells from peripheral blood)","E045 (Primary T cells effector/memory enriched from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "DN_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "EM_CD4_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E037 (Primary T helper memory cells from peripheral blood 2)","E040 (Primary T helper memory cells from peripheral blood 1)","E043 (Primary T helper cells from peripheral blood)","E045 (Primary T cells effector/memory enriched from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "EM_CD8_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E045 (Primary T cells effector/memory enriched from peripheral blood)","E048 (Primary T killer memory cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "Memory_B_cells"
annotations <- c("E032 (Primary B cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "Naive_B_cells"
annotations <- c("E032 (Primary B cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "Naive_CD4_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E038 (Primary T helper naive cells from peripheral blood)","E043 (Primary T helper cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "Naive_CD8_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E047 (Primary T killer naive cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "Regulatory_CD4_T_cells"
annotations <- c("E034 (Primary T cells from peripheral blood)","E043 (Primary T helper cells from peripheral blood)","E044 (Primary T regulatory cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)

celltype <- "TEMRA"
annotations <- c("E034 (Primary T cells from peripheral blood)","E045 (Primary T cells effector/memory enriched from peripheral blood)","E048 (Primary T killer memory cells from peripheral blood)")
df_list[[celltype]] <- getfunctionenrichment(background="All",interested=celltype,annotations)


library(dplyr)
result <- bind_rows(df_list, .id = "celltype")

result$group_num <- substr(result$group, 1, 4)
result$group_name <- substr(result$group, 5, nchar(result$group))

result$name <- sub(".*\\(([^)]+)\\).*", "\\1", result$name)
result$group_name <- sub(".*\\(([^)]+)\\).*", "\\1", result$group_name)
result$group_name <- sub("Primary ", "", result$group_name)
result$group_name <- sub("from peripheral blood", "", result$group_name)
result$group <- paste0(result$group_name," (",result$group_num,")")

write.csv(result,"/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/ChromHMM_sc_vs_bulk.csv")

##########plotting 
result <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/1_csvfiles/ChromHMM_sc_vs_bulk.csv")

result$celltype <- gsub("_", " ", result$celltype)
result$celltype <- factor(result$celltype, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Regulatory CD4 T cells","Cytotoxic CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","MAIT and GammaDelta T cells","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","Nonclassical Monocytes","CD56Bright NK cells","CD56Dim NK cells"))

result_sig <- result
result_sig$or[result_sig$adjp > 0.05] <- NA

ggplot(result_sig, aes(name, group)) + 
  geom_tile(aes(fill=log2(or))) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  ggtitle(paste0("cell-type-level eQTL vs. All cell eQTL"))+facet_grid(celltype~.,scales = "free_y", space = "free_y",switch = "y") + 
  scale_y_discrete(position = "right")+xlab("Chromatin state")+
  scale_fill_gradientn(colours = c("white","blue"),guide = "colorbar",na.value = "grey90")+theme(legend.position = "top",strip.text.y.left = element_text(angle = 0))+ylab("Cell type of epigenomic reference")

ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/0_plots/ChromHMM_sc_vs_bulk.pdf",width=9,height=8.5)

# result_memory_B_cells <- result[result$celltype == "Memory B cells",]
# result_memory_B_cells$adjp <- as.numeric(result_memory_B_cells$adjp)
# result_memory_B_cells_sig <- result_memory_B_cells[result_memory_B_cells$adjp < 0.05,]
# 
# ggplot(result_memory_B_cells_sig,aes(y = fct_rev(name))) +
#   theme_classic()+
#   geom_point(aes(x=or), shape=15, size=3) +
#   geom_linerange(aes(xmin=CIl, xmax=CIu)) +
#   geom_vline(xintercept = 0, linetype="dashed") +
#   labs(x="Odds Ratio", y="")

result_CD56Dim_NK_cells <- result[result$celltype == "CD56Dim NK cells",]
result_CD56Dim_NK_cells$adjp <- as.numeric(result_CD56Dim_NK_cells$adjp)
result_CD56Dim_NK_cells_sig <- result_CD56Dim_NK_cells[result_CD56Dim_NK_cells$adjp < 0.05,]

ggplot(result_CD56Dim_NK_cells_sig[result_CD56Dim_NK_cells_sig$nOverlap > 10,],aes(x = factor(name,levels=c("Genic enhancers","Enhancers","Flanking Active TSS","Active TSS")),y=or,ymin=CIl,ymax=CIu)) +
  theme_classic()+
  geom_point(size=2) +
  geom_errorbar(width=0.2)+
  geom_hline(yintercept = 0, linetype=2) +
  labs(y="Odds Ratio", x="")+coord_flip()

ggplot(result_CD56Dim_NK_cells_sig,aes(x = factor(name,levels=c("Bivalent Enhancer","Genic enhancers","Enhancers","Transcr. at gene 5' and 3'","Flanking Active TSS","Active TSS")),y=log2(or),ymin=log2(CIl),ymax=log2(CIu))) +
  theme_classic()+
  geom_point(size=2) +
  geom_errorbar(width=0.2)+
  geom_hline(yintercept = 0, linetype=2) +
  labs(y="log2(Odds Ratio)", x="")+coord_flip()

ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/0_plots/ChromHMM_sc_vs_bulk_CD56DimNKcells.pdf",width=5,height=2)

#### comparing allcells and single cell eQTL with annotation with all cells - not use

# annotations <- c("E062 (Primary mononuclear cells from peripheral blood)")
# 
# background <- conditionalsig_lift[conditionalsig_lift$sc_or_bulk == "All",]
# interested <- conditionalsig_lift[conditionalsig_lift$sc_or_bulk == "singlecell",]
# print(paste0("num background snps: ",nrow(background)))
# b19pos.background <- xSNPlocations(background$SNP, GR.SNP = "dbSNP_Common",RData.location = "/software/team282/XGR/bigdata")
# print(paste0("num interested snps: ",nrow(interested)))
# b19pos.intstested <- xSNPlocations(interested$SNP, GR.SNP = "dbSNP_Common",RData.location = "/software/team282/XGR/bigdata")
# GR.annotations_annot <- GR.annotations[annotations]
# #GR.annotations_annot <- get(load("/software/team282/XGR/bigdata/EpigenomeAtlas_15Segments_E062.RData"))
# ls_df <- lapply(1:length(GR.annotations_annot), function(i){
#   GR.annotation_use <- GR.annotations_annot[i]
#   message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation_use),
#                   as.character(Sys.time())), appendLF=T)
#   df <- xGRviaGenomicAnno(data.file=b19pos.intstested, 
#                           background.file = b19pos.background,
#                           format.file="GRanges",
#                           p.tail="two-tails",
#                           GR.annotation=GR.annotation_use, 
#                           RData.location=RData.location, verbose=F)
#   df$group <- names(GR.annotation_use)
#   return(df)
# })
# df <- do.call(rbind, ls_df)
# for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
# df$name <- gsub(" from peripheral blood", "", df$name)
# 
# p1 <- ggplot(df, aes(name, group)) + 
#   geom_tile(aes(fill=log2(fc))) + 
#   scale_fill_gradient2(low="red", mid="white", high="blue") +
#   theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
#   ggtitle(paste0("sc-eQTL vs all cells eQTL"))
# p1
# 
# df$adjp <- as.numeric(df$adjp)
# 
# ggplot(df,aes(y = fct_rev(name))) + 
#   theme_classic()+
#   geom_point(aes(x=or), shape=15, size=3) +
#   geom_linerange(aes(xmin=CIl, xmax=CIu)) +
#   geom_vline(xintercept = 0, linetype="dashed") +
#   labs(x="Odds Ratio", y="")
# 
# 
# df$fc[which(df$adjp > 0.05)] <- NA
# 
# p2 <- ggplot(df, aes(name, group)) + 
#   geom_tile(aes(fill=log2(fc))) + 
#   scale_fill_gradient2(low="red", mid="white", high="blue") +
#   theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
#   ggtitle(paste0("sc-eQTL vs all cells eQTL"))
# p2
# 
# df$adjp <- as.numeric(df$adjp)
# df_sig <- df[df$adjp < 0.05,]
# 
# ggplot(df_sig,aes(y = fct_rev(name))) + 
#   theme_classic()+
#   geom_point(aes(x=or), shape=15, size=3) +
#   geom_linerange(aes(xmin=CIl, xmax=CIu)) +
#   geom_vline(xintercept = 0, linetype="dashed") +
#   labs(x="Odds Ratio", y="")



