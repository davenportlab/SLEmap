####plot heatmap from mashr results
library(reshape2)
library(pheatmap)

##get cell type pairs
pairs <- read.table("/path/mashresults/pairs.tsv")

##bind all mashr results in proportion of sig hits shared
mashresults <- data.frame()
for(i in 1:nrow(pairs)){
  celltypes <- paste0(pairs$V1[i],"_",pairs$V2[i])
  results <- read.csv(paste0("/path/mashresults/output/",celltypes,"/",celltypes,"_prop_sharing.csv"))
  mashresults <- rbind(mashresults,results)
}

##prep pheatmap input
mat <- acast(mashresults, to ~ from, value.var = "prop") ### row~column
rownames(mat) <- gsub("_"," ",rownames(mat))
rownames(mat) <- gsub("cells","",rownames(mat))
rownames(mat) <- gsub("All","All-PBMC",rownames(mat))
colnames(mat) <- gsub("_"," ",colnames(mat))
colnames(mat) <- gsub("cells","",colnames(mat))
colnames(mat) <- gsub("All","All-PBMC",colnames(mat))
mat[is.na(mat)] <- 1

##plot with col cluster
ph <- pheatmap(mat, cluster_rows = F, cluster_cols = T, silent = TRUE)
ord <- ph$tree_col$order
mat_ord <- mat[ord, ord]
pdf("/path/mashresults/mash_heatmap_colcluster.pdf", width = 6.7, height = 6.6)
pheatmap(
  mat_ord,
  cluster_rows = FALSE,
  cluster_cols = T,
  color = colorRampPalette(c("white", "#1F5A9E"))(100),
  angle_col = 45, main="mash sharing: x is identification, y is replication"
)## x-axis is identification, y-axis is replication
dev.off()

##plot with row cluster
ph <- pheatmap(mat, cluster_rows = T, cluster_cols = F, silent = TRUE)
ord <- ph$tree_row$order
mat_ord <- mat[ord, ord]
pdf("/path/mashresults/mash_heatmap_rowcluster.pdf", width = 8, height = 6.6)
pheatmap(
  mat_ord,
  cluster_rows = T,
  cluster_cols = F,
  color = colorRampPalette(c("white", "#1F5A9E"))(100),
  angle_col = 45, main="mash sharing: x is identification, y is replication"
)## x-axis is identification, y-axis is replication
dev.off()

