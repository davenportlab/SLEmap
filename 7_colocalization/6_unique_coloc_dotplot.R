### examine coloc results
library(dplyr)
library(tidyverse)
library(data.table)

##read in coloc results
GWAS_ID <- "GCST90270940"
coloc_results <- read.table(paste0("/path/coloc/coloc_output_otar.csv"),sep=",",fill=T,header=T,row.names = 1)

coloc_results$otar_overlap_forplot[coloc_results$otar_overlap == T] <- "Known"
coloc_results$otar_overlap_forplot[coloc_results$otar_overlap == F] <- "Novel"
coloc_results$otar_overlap_forplot <- factor(coloc_results$otar_overlap_forplot,levels=c("Novel","Known"))
coloc_results$celltype_forplots[is.na(coloc_results$celltype_forplots) == T] <- "All"
coloc_results$cellgroup_forplots[is.na(coloc_results$cellgroup_forplots) == T] <- "All"
coloc_results$celltype_forplots <- factor(coloc_results$celltype_forplots, levels=c("All","Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells"))
coloc_results$cellgroup_forplots <- factor(coloc_results$cellgroup_forplots, levels=c("All","CD4 T","CD8 T","Other T","B","Mono","NK"))

####bar plot for coloc per cell type
coloc_results_onegene_percelltype <- coloc_results %>%
  distinct(cell_type, gene_symbol, .keep_all = TRUE)
coloc_results_onegene_percelltype <- as.data.frame(table(coloc_results_onegene_percelltype$celltype_forplots))
coloc_results_onegene_percelltype$cellgroup_forplots <- coloc_results$cellgroup_forplots[match(coloc_results_onegene_percelltype$Var1,coloc_results$celltype_forplots)]
coloc_results_onegene_percelltype$cellgroup_forplots <- factor(coloc_results_onegene_percelltype$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))
ggplot(coloc_results_onegene_percelltype, aes(x=fct_rev(Var1),y=Freq))+geom_bar(stat = "identity")+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +facet_grid(cellgroup_forplots~., scales = "free", space = "free")+xlab("") +ylab("Number of colocalised genes")+scale_fill_brewer(palette="Set2")+coord_flip()
ggsave("/path/coloc/barplot_allcolocgenes.pdf",height=5,width=5)

## dotplot with PPH4
library("viridis")

###size is pph4 and color is novel
ggplot(coloc_results, 
       aes(x = reorder(gwashit_forplot, Ord1), 
           y = celltype_forplots)) +
  geom_point(aes(size=PP.H4.abf, color=otar_overlap_forplot)) +  
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(. ~ cellgroup_forplots, scales = "free", space = "free") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("GWAS hit") +
  ylab("Cell type") +
  coord_flip()+
  scale_size(limits = c(0.8, 1),range = c(1, 4))+ 
  scale_color_manual(values = c("Novel" = "#5FA1ED","Known" = "#AEC4D6"))+
  ggtitle(paste0("Coloc with ", GWAS_ID))+ theme(panel.spacing = unit(-0.1, "lines"))

ggsave(paste0("/path/coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene_novel_color.pdf"),height=14,width=9)