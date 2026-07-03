### examine coloc results
library(dplyr)
library(tidyverse)
library(data.table)

##read in coloc results
#GWAS_ID <- "GCST003156"
GWAS_ID <- "GCST90270940"
coloc_results <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/coloc_output_otar.csv"),sep=",",fill=T,header=T,row.names = 1)

ggplot(coloc_results,aes(x=reorder(gwashit_forplot, Ord1), y = celltype_forplots, color = abs(slope_forplotting), size = -log10(lead_H4_variant_pval)))+
  theme_classic() + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_point(shape = 1,colour = "grey")+
  facet_grid(.~cellgroup_forplots, scales = "free", space = "free")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  scale_color_gradient2(low = "blue",mid = "yellow",high = "darkgreen",midpoint = 0)+
  scale_size(range = c(3, 5))+xlab("GWAS hit") +ylab("Cell type")+ coord_flip()+ggtitle(paste0("Coloc with ",GWAS_ID," (PP.H4>0.8)(dots=lead PPH4 eQTLs>0.8)"))

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
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/barplot_allcolocgenes.pdf",height=5,width=5)

## dotplot with PPH4
library("viridis")

ggplot(coloc_results, 
       aes(x = reorder(gwashit_forplot, Ord1), 
           y = celltype_forplots, 
           fill = PP.H4.abf)) +
  geom_point(shape = 21, stroke = 0,size=4) +
  scale_fill_viridis(option = "D", direction = -1) +  # for the fill
  scale_colour_identity() +                           # for literal outline colors
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(otar_overlap_forplot ~ cellgroup_forplots, scales = "free", space = "free") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("GWAS hit") +
  ylab("Cell type") +
  coord_flip() +
  ggtitle(paste0("Coloc with ", GWAS_ID))+ theme(panel.spacing = unit(-0.1, "lines"))

ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique.pdf"),height=14,width=10.5)

### by gene symbol only
ggplot(coloc_results, 
       aes(x = gene_symbol, 
           y = celltype_forplots, 
           fill = PP.H4.abf)) +
  scale_x_discrete(limits=rev)+
  geom_point(shape = 21, stroke = 0,size=4) +
  scale_fill_viridis(option = "D", direction = -1) +  # for the fill
  scale_colour_identity() +                           # for literal outline colors
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(otar_overlap_forplot ~ cellgroup_forplots, scales = "free", space = "free") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("Colocalized gene") +
  ylab("Cell type") +
  coord_flip() +
  ggtitle(paste0("Coloc with ", GWAS_ID))+ theme(panel.spacing = unit(-0.1, "lines"))

ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene.pdf"),height=14,width=9)

ggplot(coloc_results, 
       aes(x = gene_symbol, 
           y = celltype_forplots, 
           fill = PP.H4.abf)) +
  scale_x_discrete(limits=rev)+
  geom_point(shape = 21, stroke = 0,size=4) +
  scale_fill_viridis(option = "D", direction = -1) +  # for the fill
  scale_colour_identity() +                           # for literal outline colors
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(cellgroup_forplots ~ otar_overlap_forplot, scales = "free", space = "free") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("Colocalized gene") +
  ylab("Cell type") +
  ggtitle(paste0("Coloc with ", GWAS_ID))+ theme(panel.spacing = unit(-0.1, "lines"))

ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene_horizontal.pdf"),height=7,width=14)

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


ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene_novel_color.pdf"),height=14,width=9)

### gwas hit and gene separately
coloc_results$gwas_hit <- gsub("_",":",coloc_results$gwas_hit)
coloc_results$celltype_forplots <- fct_rev(coloc_results$celltype_forplots)
ggplot(coloc_results[coloc_results$otar_overlap_forplot == "Novel",], 
       aes(x = gene_symbol, 
           y = celltype_forplots)) +
  geom_point(aes(size=PP.H4.abf, color=otar_overlap_forplot)) +  
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) +
  facet_grid(cellgroup_forplots ~ reorder(gwas_hit,-Ord1), scales = "free", space = "free", switch = 'x') +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylab("Cell type") +
  scale_size(limits = c(0.8, 1),range = c(1, 4))+ 
  scale_color_manual(values = c("Novel" = "#5FA1ED","Known" = "#AEC4D6"))+
  ggtitle(paste0("Coloc with ", GWAS_ID))+ theme(strip.placement = "outside",panel.spacing = unit(0, "lines"),panel.border = element_rect(fill = NA, color = "black", linetype = "dashed"),strip.text.x.bottom = element_text(angle = 90))+scale_x_discrete(position = "top")+theme(legend.position = "top")

ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene_novel_color_gwashit_grouped_onlynovel.pdf"),height=6,width=10)

ggplot(coloc_results, 
       aes(x = gene_symbol, 
           y = celltype_forplots)) +
  geom_point(aes(size=PP.H4.abf, color=otar_overlap_forplot)) +  
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(reorder(gwas_hit,-Ord1) ~ cellgroup_forplots, scales = "free", space = "free", switch = 'y') +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("GWAS hit") +
  ylab("Cell type") +
  coord_flip()+
  scale_size(limits = c(0.8, 1),range = c(1, 4))+ 
  scale_color_manual(values = c("Novel" = "#5FA1ED","Known" = "#AEC4D6"))+
  ggtitle(paste0("Coloc with ", GWAS_ID))+theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+theme(strip.text.y.left = element_text(angle = 0))

ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene_novel_color_gwashit_grouped.pdf"),height=14,width=9)

ggplot(coloc_results, 
       aes(x = reorder(gwashit_forplot, Ord1), 
           y = celltype_forplots)) +
  scale_x_discrete(limits=rev)+
  geom_point(aes(size=PP.H4.abf, color=otar_overlap_forplot)) +  
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(cellgroup_forplots ~ ., scales = "free", space = "free") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("GWAS hit") +
  ylab("Cell type") + theme(legend.position = "top")+
  scale_size(limits = c(0.8, 1),range = c(1, 4))+ 
  scale_color_manual(values = c("Novel" = "#5FA1ED","Known" = "#AEC4D6"))+
  ggtitle(paste0("Coloc with ", GWAS_ID))+ theme(panel.spacing = unit(-0.1, "lines"))


ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene_novel_color_horizontal.pdf"),height=8,width=15)

ggplot(coloc_results, 
       aes(x = gene_symbol, 
           y = celltype_forplots)) +
  geom_point(aes(size=PP.H4.abf, color=otar_overlap_forplot)) +  
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(cellgroup_forplots ~ ., scales = "free", space = "free") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab("GWAS hit") +
  ylab("Cell type") + theme(legend.position = "top")+
  scale_size(limits = c(0.8, 1),range = c(2, 5))+ 
  scale_color_manual(values = c("Novel" = "#5FA1ED","Known" = "#AEC4D6"))+
  ggtitle(paste0("Coloc with ", GWAS_ID))+ theme(panel.spacing = unit(-0.1, "lines"))


ggsave(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/colocsummary_dotplot_",GWAS_ID,"_PPH4_checksigeQTL_unique_bygene_novel_color_horizontal_genesymbol.pdf"),height=8,width=15)
