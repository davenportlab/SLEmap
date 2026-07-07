####making plots comparing onek1k and slemap coloc results
library(ggplot2)
library(ggpubr)

GWAS_ID="GCST90270940"
OUTPUT_DIR <- "/path/onek1k_locus_breaker_coloc/"
####compare with SLEmap
onek1k<- read.csv(paste0("/path/onek1k_locus_breaker_coloc/1_csvfiles/coloc_sigresults_GCST90270940.csv"),header=T)
SLEmap <- read.csv("/path/locus_breaker_coloc/outputs/1_csvfiles/coloc_sigresults_GCST90270940.csv")
SLEmap$celltype_forplots[is.na(SLEmap$celltype_forplots) == T] <- "All"

####barplot for number of colocs per cell type
numcoloc <- as.data.frame(table(SLEmap$cell_type))
onek1k_num <- as.data.frame(table(onek1k$cell_type))
numcoloc <- merge(numcoloc,onek1k_num,by="Var1",all=T)
colnames(numcoloc) <- c("cell_type","SLEmap","Onek1k")

numcoloc$cellgroup <- NA
numcoloc$cellgroup[numcoloc$cell_type %in% c("CD56Bright_NK_cells","CD56Dim_NK_cells")] <- "NK"
numcoloc$cellgroup[numcoloc$cell_type %in% c("Classical_Monocytes","Nonclassical_Monocytes")] <- "Mono"
numcoloc$cellgroup[numcoloc$cell_type %in% c("CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","Cytotoxic_CD4_T_cells")]<- "CD4_T"
numcoloc$cellgroup[numcoloc$cell_type %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
numcoloc$cellgroup[numcoloc$cell_type %in% c("MAIT_and_GammaDelta_T_cells","DN_T_cells")] <- "Other_T"
numcoloc$cellgroup[numcoloc$cell_type %in% c("Memory_B_cells","Naive_B_cells")] <- "B"
numcoloc$cellgroup[numcoloc$cell_type %in% c("All")] <- "All"

numcoloc$celltype_forplots <- gsub("_", " ", numcoloc$cell_type)
numcoloc$celltype_forplots <- factor(numcoloc$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Regulatory CD4 T cells","Cytotoxic CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","MAIT and GammaDelta T cells","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","Nonclassical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
numcoloc$cellgroup_forplots <- gsub("_", " ", numcoloc$cellgroup)
numcoloc$cellgroup_forplots <- factor(numcoloc$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))

library(tidyr)

numcoloc_long <- numcoloc %>%
  pivot_longer(
    cols = c(SLEmap, Onek1k),
    names_to = "Study",
    values_to = "colocs"
  )

ggplot(numcoloc_long, aes(x=celltype_forplots,y=colocs,fill=Study))+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +facet_grid(~cellgroup_forplots,scales = "free_x", space = "free_x")+xlab("") +ylab("Number of colocalized GWAS hits")+scale_fill_brewer(palette="Set2")+geom_bar(stat="identity", position=position_dodge())+geom_text(aes(label=colocs), vjust=1.6,position = position_dodge(0.9), size=3.5)
ggsave(paste0(OUTPUT_DIR,"0_plots/SLEmapcomparison_barplot_allcelltypes.pdf"),height=6,width=10)

numcoloc_long_both <- numcoloc_long[!numcoloc_long$cell_type %in% c("MAIT_and_GammaDelta_T_cells","DN_T_cells","Cytotoxic_CD4_T_cells","Nonclassical_Monocytes"),]
ggplot(numcoloc_long_both, aes(x=celltype_forplots,y=colocs,fill=Study))+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +facet_grid(~cellgroup_forplots,scales = "free_x", space = "free_x")+xlab("") +ylab("Number of colocalized GWAS hits")+scale_fill_brewer(palette="Set2")+geom_bar(stat="identity", position=position_dodge())+geom_text(aes(label=colocs), vjust=1.6,position = position_dodge(0.9), size=3.5)
ggsave(paste0(OUTPUT_DIR,"0_plots/SLEmapcomparison_barplot.pdf"),height=6,width=8)

####identifying identical signals between the datasets
### figuring out whether slemap colocs are "Not tested for eQTL","Tested, no significant eQTL","Tested for colocalization, not colocalized","Colocalized" in onek1k
SLEmap$genestested_onek1k <- "Nottested"
genestested <- read.csv("/path/Onek1k/eQTLmapping/alltested_genes_onek1k.csv",header=T)
genestested$use <- paste0(genestested$phenotype_id,genestested$celltype)
SLEmap$genestested_onek1k[paste0(SLEmap$gene_id,SLEmap$cell_type) %in% genestested$use] <- "tested"

SLEmap$genessig_onek1k <- "Notsig"
genessig <- read.csv("/path/_Onek1k/eQTLmapping/allsig_genes_onek1k.csv",header=T)
genessig$use <- paste0(genessig$phenotype_id,genessig$celltype)
SLEmap$genessig_onek1k[paste0(SLEmap$gene_id,SLEmap$cell_type) %in% genessig$use] <- "Sig"

##check if significant eQTLs were tested for coloc in the cell type
coloc_all_files <- list.files(paste0("/path/onek1k_locus_breaker_coloc/GCST90270940/all_checksigeQTL_checkallele/"),"*.txt",full.names=T)
coloc_all_files <- coloc_all_files[!coloc_all_files == "/path/onek1k_locus_breaker_coloc/GCST90270940/all_checksigeQTL_checkallele/coloc_output_with_gene_name.txt"]
coloc_all_results <- do.call(rbind, lapply(coloc_all_files, read.table, header = TRUE, sep = "\t"))

SLEmap$testedcoloc_onek1k <- coloc_all_results$status[match(paste0(SLEmap$cell_type,SLEmap$gene_id,SLEmap$gwas_hit),paste0(coloc_all_results$cell_type,coloc_all_results$gene_id,coloc_all_results$gwas_hit))]
SLEmap$testedcoloc_onek1k[SLEmap$genessig_onek1k == "Sig" & SLEmap$genestested_onek1k == "tested" & is.na(SLEmap$testedcoloc_onek1k) == T] <- "lead eqtl not in window"

SLEmap$coloc_onek1k[paste0(SLEmap$gene_id,SLEmap$cell_type,SLEmap$gwas_hit) %in% paste0(onek1k$gene_id,onek1k$cell_type,onek1k$gwas_hit)] <- "coloc"

SLEmap$group_onek1k <- paste0(SLEmap$genestested_onek1k,"_",SLEmap$genessig_onek1k,"_",SLEmap$testedcoloc_onek1k,"_",SLEmap$coloc_onek1k)
unique(SLEmap$group_onek1k)

SLEmap$group_onek1k_detail <- SLEmap$group_onek1k
SLEmap$group_onek1k[SLEmap$group_onek1k == "Nottested_Notsig_NA_NA"] <- "Not tested for eQTL"
SLEmap$group_onek1k[SLEmap$group_onek1k == "tested_Notsig_NA_NA"] <- "Tested, no significant eQTL"
SLEmap$group_onek1k[SLEmap$group_onek1k %in% c("tested_Sig_lead eqtl not in window_NA","tested_Sig_minimum eQTL p value among common SNPs is bigger than nominal pval threshold_NA")] <- "Significant eQTL, not tested for colocalization"
SLEmap$group_onek1k[SLEmap$group_onek1k == "tested_Sig_tested for coloc_NA"] <- "Tested for colocalization, not colocalized"
SLEmap$group_onek1k[SLEmap$group_onek1k == "tested_Sig_tested for coloc_coloc"] <- "Colocalized"
SLEmap$group_onek1k <- factor(SLEmap$group_onek1k,levels=c("Not tested for eQTL","Tested, no significant eQTL","Significant eQTL, not tested for colocalization","Tested for colocalization, not colocalized","Colocalized"))

SLEmap$celltype_forplots <- factor(SLEmap$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Regulatory CD4 T cells","Cytotoxic CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
SLEmap$cellgroup_forplots <- gsub("_", " ", SLEmap$cellgroup)
SLEmap$cellgroup_forplots <- factor(SLEmap$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))

## dotplot with slemap and onek1k
ggplot(SLEmap[!SLEmap$cell_type %in% c("Cytotoxic_CD4_T_cells","DN_T_cells"),],aes(x=reorder(gwashit_forplot, Ord1), y = celltype_forplots))+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_point(shape = 1,colour = "black",size=4)+
  geom_point(aes(color=group_onek1k),shape=16,size=3)+
  scale_color_manual(values=c(
    "#D73027",  # red
    "#F4C430",  # yellow
    "#1A9850",  # green
    "#2C7BB6",  # blue
    "#7B3294"   # purple
  ))+
  facet_grid(.~cellgroup_forplots, scales = "free", space = "free")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+xlab("GWAS hit") +ylab("Cell type")+ coord_flip()+ggtitle(paste0("Coloc with ",GWAS_ID))
ggsave(paste0(OUTPUT_DIR,"0_plots/SLEmap_dotplot_with_Onek1k_shapeplot.pdf"),width=10,height=14)

ggplot(SLEmap[!SLEmap$cell_type %in% c("Cytotoxic_CD4_T_cells","DN_T_cells"),],aes(x=reorder(gwashit_forplot, Ord1), y = celltype_forplots))+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_point(shape = 1,colour = "black",size=4)+
  geom_point(aes(color=group_onek1k),shape=16,size=3)+
  scale_color_manual(values=c(
    "#D73027",  # red
    "#F4C430",  # yellow
    "#1A9850",  # green
    "#2C7BB6",  # blue
    "#7B3294"   # purple
  ))+
  facet_grid(cellgroup_forplots~., scales = "free", space = "free")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+xlab("GWAS hit") +ylab("Cell type")+ggtitle(paste0("Coloc with ",GWAS_ID))+theme(legend.position = "top")
ggsave(paste0(OUTPUT_DIR,"0_plots/SLEmap_dotplot_with_Onek1k_shapeplot_horizontal.pdf"),width=14,height=7)


SLEmap$group_onek1k_forplot0[SLEmap$group_onek1k == "Not tested for eQTL"] <- "Colocalized only in SLEmap"
SLEmap$group_onek1k_forplot0[SLEmap$group_onek1k == "Tested, no significant eQTL"] <- "Colocalized only in SLEmap"
SLEmap$group_onek1k_forplot0[SLEmap$group_onek1k %in% c("Significant eQTL, not tested for colocalization","Tested for colocalization, not colocalized")] <- "Colocalized only in SLEmap"
SLEmap$group_onek1k_forplot0[SLEmap$group_onek1k == "Colocalized"] <- "Colocalized in both\nSLEmap and Onek1k"

SLEmap$group_onek1k_forplot[SLEmap$group_onek1k == "Not tested for eQTL"] <- "Low expression in Onek1k"
SLEmap$group_onek1k_forplot[SLEmap$group_onek1k == "Tested, no significant eQTL"] <- "No eQTL detected in Onek1k" 
SLEmap$group_onek1k_forplot[SLEmap$group_onek1k %in% c("Significant eQTL, not tested for colocalization","Tested for colocalization, not colocalized")] <- "Colocalised only in SLEmap\n(eQTL in both)"
SLEmap$group_onek1k_forplot[SLEmap$group_onek1k == "Colocalized"] <- "Colocalized in both\nSLEmap and Onek1k"
SLEmap$group_onek1k_forplot[SLEmap$cell_type %in% c("Cytotoxic_CD4_T_cells","DN_T_cells")] <- "Celltype not tested in Onek1k"

SLEmap_onek1k_table <- as.data.frame(table(SLEmap$group_onek1k_forplot0,SLEmap$group_onek1k_forplot))
colnames(SLEmap_onek1k_table) <- c("Grouping","Onek1k analysis","Number out of SLEmap colocs")
SLEmap_onek1k_table$`Onek1k analysis` <- factor(SLEmap_onek1k_table$`Onek1k analysis`, levels=c("Low expression in Onek1k","No eQTL detected in Onek1k","Colocalised only in SLEmap\n(eQTL in both)","Colocalized in both\nSLEmap and Onek1k","Celltype not tested in Onek1k"))

SLEmap_onek1k_table$Grouping <- factor(SLEmap_onek1k_table$Grouping, levels=c("Colocalized only in SLEmap","Colocalized in both\nSLEmap and Onek1k"))

ggplot(SLEmap_onek1k_table[SLEmap_onek1k_table$`Onek1k analysis` != "Celltype not tested in Onek1k",], aes(x=`Number out of SLEmap colocs`,y=Grouping,fill=`Onek1k analysis`)) +
  geom_col() + scale_fill_manual(values=c("#808285","#D46127","#3ABCA1","#164989","grey"))+
  ## add percentage labels
  geom_text(aes(label = `Number out of SLEmap colocs`),nudge_x =0) +
  theme_classic()+theme(legend.position = "top")
ggsave(paste0(OUTPUT_DIR,"0_plots/SLEmapcolocs_onek1k_barplot.pdf"),height=2.5,width=6)

ggplot(SLEmap_onek1k_table[SLEmap_onek1k_table$`Number out of SLEmap colocs` > 0 & SLEmap_onek1k_table$`Onek1k analysis` != "Celltype not tested in Onek1k",], aes(x=`Onek1k analysis`,y=`Number out of SLEmap colocs`))+theme_classic()+geom_bar(stat="identity",fill="lightblue")+geom_text(aes(label=`Number out of SLEmap colocs`), vjust=0.2,hjust=-0.5,position = position_dodge(0.9), size=3.5)+coord_flip()+ylim(0,60)
ggsave(paste0(OUTPUT_DIR,"0_plots/SLEmapcolocs_onek1k_barplot_tidy.pdf"),height=3,width=6.5)


write.csv(SLEmap,paste0(OUTPUT_DIR,"1_csvfiles/SLEmap_Onek1kgroup.csv"))
SLEmap <- read.csv(paste0(OUTPUT_DIR,"1_csvfiles/SLEmap_Onek1kgroup.csv"),row.names=1)

SLEmap$gene_symbol[SLEmap$group_onek1k == "Tested for colocalization, not colocalized"]

###plot by cell type

SLEmap_table <- as.data.frame(table(SLEmap$cellgroup_forplots,SLEmap$celltype_forplots,SLEmap$group_onek1k_forplot))
df <- SLEmap_table[!SLEmap_table$Var2 %in% c("Cytotoxic CD4 T cells", "DN T cells"), ]
df <- df[!df$Freq == 0,]
df$Var3 <- factor(df$Var3, levels=c("Low expression in Onek1k","No eQTL detected in Onek1k","Colocalised only in SLEmap\n(eQTL in both)","Colocalized in both SLEmap and Onek1k"))
df$Var2 <- factor(df$Var2, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All"))
df$Var1 <- factor(df$Var1, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))
ggplot(data = df, aes(x = Var2, y = Freq, fill = Var3)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(~ Var1, scales = "free_x", space = "free_x")+scale_fill_manual(values=c("#808285","#D46127","#3ABCA1","#164989"))+xlab("Cell type")+ylab("Proportion")
ggsave("/path/onek1k_locus_breaker_coloc/0_plots/SLEmapcolocs_Onek1k_bycelltype.pdf",width=10,height=5)