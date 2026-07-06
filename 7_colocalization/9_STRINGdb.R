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
  score_threshold = 700,  
  input_directory = ""
)

# ---------------------------
# 2. INPUT: your gene sets
# ---------------------------
# Replace these with your actual results per cell type
coloc <- read.csv("/path/coloc/outputs/1_csvfiles/coloc_sigresults.csv")
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
