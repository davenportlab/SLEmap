library(dplyr)
#library(biomaRt)
library(data.table)

GWAS_ID="GCST90270940"

DIR_MAIN="/path/onek1k_locus_breaker_coloc"

files <- list.files(path=paste0(DIR_MAIN, 
                                "/",
                                GWAS_ID, "/all_checksigeQTL_checkallele/"),
                    pattern=".txt",
                    full.names = T)
all_out <- c()
everything <- c()
for(fn in files){
  print(fn)
  f <- fread(fn, header = T) %>% 
    as.data.frame()
  everything <- rbind(everything,f)
  f <- f[!is.na(f$nsnps),]
  all_out <- rbind(all_out, f)
}

dim(all_out)
dim(everything)
table(everything$status)

### add gene name
genenames <- read.csv("/path/ensemblID_to_genesymbol.csv")
all_out$gene_symbol <- genenames$gene_symbols[match(all_out$gene_id,genenames$X)]

all_out$PP.H4.abf <- as.numeric(all_out$PP.H4.abf)

write.table(all_out,
            paste0(DIR_MAIN, "/", GWAS_ID, "/all_checksigeQTL_checkallele/coloc_output_with_gene_name.txt"),
            sep="\t", quote=F, row.names = F)

sig_coloc <- dplyr::filter(all_out, PP.H4.abf >=0.8)

dim(sig_coloc) #190 20
length(unique(sig_coloc$gene_id)) #62

dim(sig_coloc[sig_coloc$cell_type != "All",]) #148
length(unique(sig_coloc[sig_coloc$cell_type != "All",]$gene_id)) #48

dplyr::filter(all_out, PP.H4.abf >=0.8) %>% 
  dplyr::group_by(cell_type) %>% tally()

write.csv(sig_coloc,
          paste0(DIR_MAIN, "/", GWAS_ID, "/sig_checksigeQTL_checkallele/coloc_output_with_gene_name_sig.csv"), quote=F, row.names = F)

####plot pph4 dot plot
sig_coloc <- read.csv(paste0(DIR_MAIN, "/", GWAS_ID, "/sig_checksigeQTL_checkallele/coloc_output_with_gene_name_sig.csv"))
dat <- sig_coloc %>% 
  mutate(., id=paste0(gene_symbol, "_", gwas_hit)) %>%
  dplyr::arrange(., gene_symbol) %>%
  dplyr::mutate(., Ord2=nrow(sig_coloc):1)
dat$gwashit_forplot <- vapply(strsplit(dat$id, "_"), function(parts) {
  paste0(parts[2], ":", parts[3], " (", parts[1], ")")
}, character(1))

dat$cellgroup <- NA
dat$cellgroup[dat$cell_type %in% c("CD56Bright_NK_cells","CD56Dim_NK_cells")] <- "NK"
dat$cellgroup[dat$cell_type %in% c("Classical_Monocytes","Nonclassical_Monocytes")] <- "Mono"
dat$cellgroup[dat$cell_type %in% c("CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells")]<- "CD4_T"
dat$cellgroup[dat$cell_type %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
dat$cellgroup[dat$cell_type %in% c("MAIT_and_GammaDelta_T_cells")] <- "Other_T"
dat$cellgroup[dat$cell_type %in% c("Memory_B_cells","Naive_B_cells")] <- "B"
dat$cellgroup[dat$cell_type %in% c("All")] <- "All"

dat$celltype_forplots <- gsub("_", " ", dat$cell_type)
dat$celltype_forplots <- factor(dat$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","MAIT and GammaDelta T cells","Naive B cells","Memory B cells","Classical Monocytes","Nonclassical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
dat$cellgroup_forplots <- gsub("_", " ", dat$cellgroup)
dat$cellgroup_forplots <- factor(dat$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))

library("viridis")
ggplot(dat,aes(x=reorder(gwashit_forplot, Ord2), y = celltype_forplots, color = PP.H4.abf, size = PP.H4.abf))+
  theme_classic() + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_viridis(option = "D",direction = -1)+
  geom_point(shape = 1,colour = "grey")+
  facet_grid(.~cellgroup_forplots, scales = "free", space = "free")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  scale_size(range = c(2, 5))+xlab("GWAS hit") +ylab("Cell type")+ coord_flip()+ggtitle(paste0("Onek1k-Coloc with ",GWAS_ID))
ggsave(paste0("/path/onek1k_locus_breaker_coloc1/0_plots/eQTLtested_dotplot_",GWAS_ID,"_PPH4.pdf"),height=14,width=10)

dat <- dat[!dat$gene_symbol %in% c("KANSL1", "KANSL1-AS1","ARL17B","FAM215B"),]
write.csv(dat,
          paste0("/path/onek1k_locus_breaker_coloc/1_csvfiles/coloc_sigresults_GCST90270940.csv"), quote=F, row.names = F)

### identify if coloc genes were tested in other cell types, if there is a significant eQTL for the cell type
genes <- unique(sig_coloc$gene_id)
cell_types <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","Nonclassical_Monocytes","MAIT_and_GammaDelta_T_cells","All")
df <- expand.grid(gene = genes, celltype = cell_types, stringsAsFactors = FALSE)

df$genestested <- "Nottested"
genestested <- read.csv("/path/Onek1k/eQTLmapping/alltested_genes_onek1k.csv",header=T)
genestested$use <- paste0(genestested$phenotype_id,genestested$celltype)
df$genestested[paste0(df$gene,df$celltype) %in% genestested$use] <- "tested"

df$genessig <- "Notsig"
genessig <- read.csv("/path/Onek1k/eQTLmapping/allsig_genes_onek1k.csv",header=T)
genessig$use <- paste0(genessig$phenotype_id,genessig$celltype)
df$genessig[paste0(df$gene,df$celltype) %in% genessig$use] <- "Sig"

##check if significant eQTLs were tested for coloc in the cell type
coloc_all_files <- list.files(paste0("/path/onek1k_locus_breaker_coloc/GCST90270940/all_checksigeQTL_checkallele/"),"*.txt",full.names=T)
coloc_all_files <- coloc_all_files[!coloc_all_files == "/path/onek1k_locus_breaker_coloc/GCST90270940/all_checksigeQTL_checkallele/coloc_output_with_gene_name.txt"]
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
df$cellgroup[df$celltype %in% c("Classical_Monocytes","Nonclassical_Monocytes")] <- "Mono"
df$cellgroup[df$celltype %in% c("CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells")]<- "CD4_T"
df$cellgroup[df$celltype %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
df$cellgroup[df$celltype %in% c("MAIT_and_GammaDelta_T_cells")] <- "Other_T"
df$cellgroup[df$celltype %in% c("Memory_B_cells","Naive_B_cells")] <- "B"
df$cellgroup[df$celltype %in% c("All")] <- "All"

df$celltype_forplots <- gsub("_", " ", df$celltype)
df$celltype_forplots <- factor(df$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","MAIT and GammaDelta T cells","Naive B cells","Memory B cells","Classical Monocytes","Nonclassical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
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

ggsave(paste0("/path/onek1k_locus_breaker_coloc/0_plots/eQTLtested_shapeplot_",GWAS_ID,".pdf"),width=8,height=14)
write.csv(df,paste0("/path/onek1k_locus_breaker_coloc/1_csvfiles/coloc_genegroup_",GWAS_ID,".csv"))




