library(STRINGdb)
library(igraph)
library(dplyr)

# ---------------------------
# 1. Initialize STRINGdb
# ---------------------------
# Human = 9606
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 700,   # high confidence (important)
  input_directory = ""
)

# ---------------------------
# 2. INPUT: your gene sets
# ---------------------------
# Replace these with your actual results per cell type
coloc <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs/1_csvfiles/coloc_sigresults_GCST90270940_checksigeQTL.csv")
gene_sets <- list(
  "B" = unique(coloc$gene_id[coloc$cellgroup== "B"]),
  "Monocyte" = unique(coloc$gene_id[coloc$cellgroup == "Mono"]),
  "CD4T" = unique(coloc$gene_id[coloc$cellgroup == "CD4_T"]),
  "CD8T" = unique(coloc$gene_id[coloc$cellgroup == "CD8_T"]),
  "NK" = unique(coloc$gene_id[coloc$cellgroup == "NK"])
)

# ---------------------------
# 3. Helper function
# ---------------------------
analyze_string_network <- function(gene_list, label) {
  
  cat("\n=========================\n")
  cat("Cell type:", label, "\n")
  cat("=========================\n")
  
  # Map genes to STRING IDs
  mapped <- string_db$map(
    data.frame(ensembl = gene_list),
    "ensembl",
    removeUnmappedRows = TRUE
  )
  
  if (nrow(mapped) < 3) {
    cat("Too few mapped genes for analysis\n")
    return(NULL)
  }
  
  # Get interactions for mapped genes
  hits <- string_db$get_interactions(mapped$STRING_id)
  
  if(nrow(hits) > 0){
    hits$celltype <- label
    return(hits)
  }
  
  # if(nrow(hits) > 0){
  #   # Build graph
  #   g <- graph_from_data_frame(
  #     hits[, c("from", "to")],
  #     directed = FALSE
  #   )
  #   
  #   # ---------------------------
  #   # Network stats
  #   # ---------------------------
  #   obs_edges <- gsize(induced_subgraph(
  #     g,
  #     vids = mapped$STRING_id
  #   ))
  #   
  #   obs_nodes <- length(mapped$STRING_id)
  #   
  #   cat("Observed nodes:", obs_nodes, "\n")
  #   cat("Observed edges:", obs_edges, "\n")
  #   
  #   # ---------------------------
  #   # Permutation test (connectivity enrichment)
  #   # ---------------------------
  #   all_ids <- string_db$proteins$protein_external_id
  #   
  #   n_iter <- 1000
  #   rand_edges <- numeric(n_iter)
  #   
  #   for (i in 1:n_iter) {
  #     rand_genes <- sample(all_ids, obs_nodes)
  #     sub_g <- induced_subgraph(g, vids = rand_genes)
  #     rand_edges[i] <- gsize(sub_g)
  #   }
  #   
  #   p_value <- mean(rand_edges >= obs_edges)
  #   
  #   cat("Connectivity enrichment p-value:", p_value, "\n")
  #   
  #   # ---------------------------
  #   # Cluster detection (modules)
  #   # ---------------------------
  #   sub_g <- induced_subgraph(g, vids = mapped$STRING_id)
  #   
  #   clusters <- cluster_louvain(sub_g)
  #   
  #   membership <- membership(clusters)
  #   
  #   cluster_table <- data.frame(
  #     gene = names(membership),
  #     cluster = membership
  #   )
  #   
  #   cat("\nClusters detected:", length(unique(membership)), "\n")
  #   
  #   # ---------------------------
  #   # Hub genes
  #   # ---------------------------
  #   deg <- degree(sub_g)
  #   hub_genes <- sort(deg, decreasing = TRUE)[1:min(10, length(deg))]
  #   
  #   hub_table <- data.frame(
  #     gene = names(hub_genes),
  #     degree = as.numeric(hub_genes)
  #   )
  #   
  #   # ---------------------------
  #   # Output
  #   # ---------------------------
  #   list(
  #     graph = sub_g,
  #     p_value = p_value,
  #     clusters = cluster_table,
  #     hubs = hub_table
  #   )
  # }

}


# ---------------------------
# 4. Run per cell type
# ---------------------------
results <- data.frame()
for(ct in names(gene_sets)){
  results <- rbind(results,analyze_string_network(gene_sets[[ct]], ct))
}
write.csv("")

results <- lapply(names(gene_sets), function(ct) {
  analyze_string_network(gene_sets[[ct]], ct)
})
ct=names(results) <- names(gene_sets)

# ---------------------------
# 5. Save outputs
# ---------------------------
for (ct in names(results)) {
  if (is.null(results[[ct]])) next
  
  write.csv(results[[ct]]$clusters,
            paste0(ct, "_STRING_clusters.csv"),
            row.names = FALSE)
  
  write.csv(results[[ct]]$hubs,
            paste0(ct, "_STRING_hubs.csv"),
            row.names = FALSE)
}

cat("\nDone.\n")
