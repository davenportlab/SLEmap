### evaluate conditionally independent eQTLs
library(scales)
library(ggrepel)
library(ggpubr)
library(reshape2)

###plot number of conditional eQTLs
eGenesummary <- read.csv("/path/eQTLresults/1_csvfiles/eGenesummary.csv")
resultsdir <- "/path/eQTLresults/"

conditionalsig <- as.data.frame(matrix(ncol=6))
colnames(conditionalsig) <- c("phenotype_id","celltype","variant_id","start_distance","af","rank")

for (i in 1:nrow(eGenesummary)){
  celltype <- eGenesummary$celltype[[i]]
  conditional <- read.table(paste0(resultsdir,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
  conditional$celltype <- celltype
  eGenesummary[i, 'conditional_1'] <- nrow(conditional[conditional$rank ==1,])
  eGenesummary[i, 'conditional_2'] <- nrow(conditional[conditional$rank ==2,])
  eGenesummary[i, 'conditional_3'] <- nrow(conditional[conditional$rank ==3,])
  eGenesummary[i, 'conditional_4'] <- nrow(conditional[conditional$rank ==4,])
  eGenesummary[i, 'conditional_5'] <- nrow(conditional[conditional$rank ==5,]) #max was 5
  conditionalsig <- rbind(conditionalsig,conditional[,c("phenotype_id","celltype","variant_id","start_distance","af","rank")])
}

conditionalsig <- conditionalsig[-1,]
table(conditionalsig$rank)

conditional_plot <- melt(eGenesummary[c("celltype","conditional_1","conditional_2","conditional_3","conditional_4","conditional_5")], id.vars = c("celltype"), variable.name = "conditional_rank")
conditional_plot$cellgroup <- eGenesummary$cellgroup[match(conditional_plot$celltype,eGenesummary$celltype)]
conditional_plot$conditional_rank <- factor(conditional_plot$conditional_rank,levels=c("conditional_5","conditional_4","conditional_3","conditional_2","conditional_1"))

conditional_plot$cellgroup <- NA
conditional_plot$cellgroup[conditional_plot$celltype %in% c("CD56Bright_NK_cells","CD56Dim_NK_cells")] <- "NK"
conditional_plot$cellgroup[conditional_plot$celltype %in% c("Classical_Monocytes")] <- "Mono"
conditional_plot$cellgroup[conditional_plot$celltype %in% c("CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells")]<- "CD4_T"
conditional_plot$cellgroup[conditional_plot$celltype %in% c("CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA")]<- "CD8_T"
conditional_plot$cellgroup[conditional_plot$celltype %in% c("DN_T_cells")] <- "Other_T"
conditional_plot$cellgroup[conditional_plot$celltype %in% c("Memory_B_cells","Naive_B_cells")] <- "B"

conditional_plot$celltype_forplots <- gsub("_", " ", conditional_plot$celltype)
conditional_plot$celltype_forplots <- factor(conditional_plot$celltype_forplots, levels=c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Cytotoxic CD4 T cells","Regulatory CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","All cells"))
conditional_plot$cellgroup_forplots <- gsub("_", " ", conditional_plot$cellgroup)
conditional_plot$cellgroup_forplots <- factor(conditional_plot$cellgroup_forplots, levels=c("CD4 T","CD8 T","Other T","B","Mono","NK","All"))


##add all cells
allcells_conditional <- read.table("/path/eQTLresults/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",header=T)
allcells_table <- as.data.frame(table(allcells_conditional$rank))
allcells_table$Freq <- as.numeric(allcells_table$Freq)
conditional_plot <- rbind(conditional_plot,c("All_cells","conditional_1",allcells_table$Freq[allcells_table$Var1 == 1],"All","All cells","All"))
conditional_plot <- rbind(conditional_plot,c("All_cells","conditional_2",allcells_table$Freq[allcells_table$Var1 == 2],"All","All cells","All"))
conditional_plot <- rbind(conditional_plot,c("All_cells","conditional_3",allcells_table$Freq[allcells_table$Var1 == 3],"All","All cells","All"))
conditional_plot <- rbind(conditional_plot,c("All_cells","conditional_4",allcells_table$Freq[allcells_table$Var1 == 4],"All","All cells","All"))
conditional_plot <- rbind(conditional_plot,c("All_cells","conditional_5",sum(allcells_table$Freq[allcells_table$Var1 %in% c(5,6,7,8,9)]),"All","All cells","All"))

conditional_plot$value <- as.numeric(conditional_plot$value)
ggplot(conditional_plot, aes(x=celltype_forplots,y=value,fill=conditional_rank))+geom_bar(stat = "identity",position="stack")+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +facet_grid(~cellgroup_forplots,scales = "free_x", space = "free_x")+xlab("") +ylab("Number of conditional eQTLs")+ scale_fill_manual(values = c('#006292', '#298299', '#51a29f', '#79c2a5', '#f7cd74'),labels=c("5+","4","3","2","1"))+ggtitle("Number of eQTLs from conditional analysis")
ggsave(paste0(resultsdir,"0_plots/numsigGene_conditional_wallcells.pdf"),width=9,height=5)


