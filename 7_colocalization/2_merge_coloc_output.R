library(dplyr)
library(data.table)

GWAS_ID="GCST90270940"

####make output directories
if (!dir.exists(paste0(DIR_MAIN,"/0_plots"))) {
  dir.create(paste0(DIR_MAIN,"/0_plots"))
}
if (!dir.exists(paste0(DIR_MAIN,"/1_csvfiles"))) {
  dir.create(paste0(DIR_MAIN,"/1_csvfiles"))
}

###merge coloc
DIR_MAIN="/path/coloc/outputs"
files <- list.files(path=paste0(DIR_MAIN, 
                                "/",
                                GWAS_ID, "/all_checksigeQTL_checkallele/"),
                    pattern=".txt",
                    full.names = T)
all_out <- c()
for(fn in files){
  print(fn)
    f <- fread(fn, header = T) %>% 
      as.data.frame() %>%
      dplyr::filter(!is.na(nsnps))
    all_out <- rbind(all_out, f)
}

dim(all_out)

### add gene name
genenames <- read.csv("/path/ensemblID_to_genesymbol.csv")
all_out$gene_symbol <- genenames$gene_symbols[match(all_out$gene_id,genenames$X)]
all_out$PP.H4.abf <- as.numeric(all_out$PP.H4.abf)

### add coloc testing window from locus breaker
locusbreaker <- read.csv("/path/coloc/gwashit_locusbreaker_forinput.csv")
all_out$locusStart <- locusbreaker$locusStart[match(all_out$gwas_hit,locusbreaker$variant_id)]
all_out$locusEnd <- locusbreaker$locusEnd[match(all_out$gwas_hit,locusbreaker$variant_id)]

write.table(all_out,
            paste0(DIR_MAIN, "/", GWAS_ID, "/all_checksigeQTL_checkallele/coloc_output_with_gene_name.txt"),
            sep="\t", quote=F, row.names = F)

sig_coloc <- dplyr::filter(all_out, PP.H4.abf >=0.8)

dplyr::filter(all_out, PP.H4.abf >=0.8) %>% 
  dplyr::group_by(cell_type) %>% tally()

write.csv(sig_coloc,
            paste0(DIR_MAIN, "/", GWAS_ID, "/sig_checksigeQTL_checkallele/coloc_output_with_gene_name_sig.csv"), quote=F, row.names = F)

### identify if coloc genes were tested in other cell types, if there is a significant eQTL for the cell type
genes <- unique(sig_coloc$gene_id)
cell_types <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","DN_T_cells","Memory_B_cells","Naive_B_cells","All")
df <- expand.grid(gene = genes, celltype = cell_types, stringsAsFactors = FALSE)

df$genestested <- "Nottested"
genestested <- read.csv("/path/eQTLresults/1_csvfiles/genestested_wbulklike.csv")
genestested$use <- paste0(genestested$phenotype_id,genestested$celltype)
df$genestested[paste0(df$gene,df$celltype) %in% genestested$use] <- "tested"

df$genessig <- "Notsig"
genessig <- read.csv("/path/eQTLresults/1_csvfiles/genessig_wbulklike.csv")
genessig$use <- paste0(genessig$phenotype_id,genessig$celltype)
df$genessig[paste0(df$gene,df$celltype) %in% genessig$use] <- "Sig"

##check if significant eQTLs were tested for coloc in the cell type
coloc_all_files <- list.files(paste0(DIR_MAIN,"/",GWAS_ID,"/all_checksigeQTL_checkallele/"),"*.txt",full.names=T)
coloc_all_files <- coloc_all_files[!coloc_all_files == paste0(DIR_MAIN,"/",GWAS_ID,"/all_checksigeQTL_checkallele//coloc_output_with_gene_name.txt")]
coloc_all_results <- do.call(rbind, lapply(coloc_all_files, read.table, header = TRUE, sep = "\t"))
coloc_all_nottested <- coloc_all_results[coloc_all_results$status %in% c("Too small number of common SNPs between datasets","minimum GWAS p value among common SNPs is bigger than 1e-5","minimum eQTL p value among common SNPs is bigger than nominal pval threshold"),]
df$testedcoloc[paste0(df$gene,df$celltype) %in% paste0(coloc_all_nottested$gene_id,coloc_all_nottested$cell_type)] <- "Nottestcoloc"
coloc_all_tested <- coloc_all_results[coloc_all_results$status == "tested for coloc",]
df$testedcoloc[paste0(df$gene,df$celltype) %in% paste0(coloc_all_tested$gene_id,coloc_all_tested$cell_type)]<- "Testcoloc"
df$testedcoloc[is.na(df$testedcoloc) == T & df$genessig == "Sig"] <- "eQTLoutwindow"
df$testedcoloc[is.na(df$testedcoloc) == T] <- "Notsig2"

df$coloc <- "notcoloc"
df$coloc[paste0(df$celltype,df$gene) %in% paste0(sig_coloc$cell_type,sig_coloc$gene_id)] <- "coloc"

df$group <- paste0(df$genestested,"_",df$genessig,"_",df$testedcoloc,"_",df$coloc)
unique(df$group)

df$group_detail <- df$group
df$group[df$group == "Nottested_Notsig_Notsig2_notcoloc"] <- "Not tested for eQTL"
df$group[df$group == "tested_Notsig_Notsig2_notcoloc"] <- "Tested, no significant eQTL"
df$group[df$group %in% c("tested_Sig_Nottestcoloc_notcoloc","tested_Sig_eQTLoutwindow_notcoloc")] <- "Significant eQTL, not tested for colocalization"
df$group[df$group == "tested_Sig_Testcoloc_notcoloc"] <- "Tested for colocalization, not colocalized"
df$group[df$group == "tested_Sig_Testcoloc_coloc"] <- "Colocalized"
df$group <- factor(df$group,levels=c("Not tested for eQTL","Tested, no significant eQTL","Significant eQTL, not tested for colocalization","Tested for colocalization, not colocalized","Colocalized"))

df$cellgroup <- NA
df$cellgroup[df$celltype %in% c("CD56Bright_NK_cells","CD56Dim_NK_cells")] <- "NK"
df$cellgroup[df$celltype %in% c("Classical_Monocytes")] <- "Mono"
df$cellgroup[df$celltype %in% c("CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells")]<- "CD4_T"
df$cellgroup[df$celltype %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
df$cellgroup[df$celltype %in% c("DN_T_cells")] <- "Other_T"
df$cellgroup[df$celltype %in% c("Memory_B_cells","Naive_B_cells")] <- "B"
df$cellgroup[df$celltype %in% c("All")] <- "All"

df$celltype_forplots <- gsub("_", " ", df$celltype)
df$celltype_forplots <- factor(df$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
df$cellgroup_forplots <- gsub("_", " ", df$cellgroup)
df$cellgroup_forplots <- factor(df$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))

df$gene_symbol <- sig_coloc$gene_symbol[match(df$gene,sig_coloc$gene_id)]
dat <- sig_coloc %>% 
  mutate(., id=paste0(gene_symbol, "_", gwas_hit)) %>%
  dplyr::arrange(., gene_symbol) %>%
  dplyr::mutate(., Ord2=nrow(sig_coloc):1)

df$Ord2 <- dat$Ord2[match(df$gene,dat$gene_id)]

ggplot(df[(df$group != "Not tested for eQTL"),],aes(x=reorder(gene_symbol,Ord2), y = celltype_forplots, shape=group,color=group))+
  geom_point(size=3)+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  facet_grid(.~cellgroup_forplots, scales = "free", space = "free")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  scale_size(range = c(2, 6))+
  xlab("Gene") +
  ylab("Cell type")+ 
  coord_flip()+
  scale_shape_manual(values=c(4, 0,2,19))+theme(legend.position = "top")+
  ggtitle(paste0("SNPs for colocalized genes"))

ggsave(paste0(DIR_MAIN,"/0_plots/eQTLtested_shapeplot_",GWAS_ID,".pdf"),width=8,height=14)

write.csv(df,paste0(DIR_MAIN,"/1_csvfiles/coloc_genegroup.csv"))

