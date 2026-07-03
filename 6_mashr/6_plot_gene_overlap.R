###make mash heatmap

order_vec <- rev(c("Classical_Monocytes","Memory_B_cells","Naive_B_cells","Naive_CD4_T_cells","Naive_CD8_T_cells","EM_CD8_T_cells","CM_CD4_T_cells","EM_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","DN_T_cells","CD56Bright_NK_cells","CD56Dim_NK_cells","Cytotoxic_CD4_T_cells","TEMRA"))

mash <- read.table("/path/mashresults/output/model_summary/pairwise_sharing.txt")
mash_clean <- mash
rownames(mash_clean) <- gsub("^merged_|\\.tsv$", "", rownames(mash_clean))
colnames(mash_clean) <- gsub("^merged_|\\.tsv$", "", colnames(mash_clean))

mash_clean <- mash_clean[order_vec, order_vec]
mash_long <- melt(
  as.matrix(mash_clean),
  varnames = c("Var1", "Var2"),
  value.name = "value"
)

mash_long$Var1 <- factor(mash_long$Var1, levels = order_vec)
mash_long$Var2 <- factor(mash_long$Var2, levels = order_vec)

plot2 <- ggplot(mash_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", name = "mash") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )+theme(legend.position = "top")

#####make heatmap of num overlap genes
pairs$overlap_prop <- pairs$num_overlap/pairs$allgenestested * 100
pairs$nonoverlap_prop <- pairs$num_nonoverlap/pairs$allgenestested * 100

pairs_full <- pairs %>% bind_rows(pairs %>% rename(celltype_1 = celltype_2,celltype_2 = celltype_1))

# prep
gene_overlap2 <- pairs_full %>%
  filter(celltype_1 != "All", celltype_2 != "All")

# factor order
gene_overlap2$celltype_1 <- factor(gene_overlap2$celltype_1, levels = order_vec)
gene_overlap2$celltype_2 <- factor(gene_overlap2$celltype_2, levels = order_vec)

# plot

plot3 <- ggplot(gene_overlap2, aes(celltype_1, celltype_2, fill =overlap_prop)) +
  geom_tile() +
  scale_fill_gradient(
    high = "red",
    low = "white",
    name = "Gene overlap_prop"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )+theme(legend.position = "top")


ggarrange(plot2,plot3)
ggsave("/path/mashresults/heatmaps_mash_geneoverlap.pdf",width=12,height=6)
