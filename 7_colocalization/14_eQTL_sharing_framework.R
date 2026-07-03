# eQTL sharing framework

setwd("/path/flanders")

library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(tibble)
library(pheatmap)
library(igraph)

# read in files
loci <- list.files("flanders_updated_outputs/results/coloc_info_tables/", full.names = T)
finemapped.loci <- lapply(loci, fread)
finemapped.loci <- rbindlist(finemapped.loci)

coloc.results <- read.delim("flanders_updated_outputs/results/coloc/coloc_run_colocalization.table.all.tsv")

coloc.info.tables <- list.files("flanders_updated_outputs/results/coloc_info_tables/", full.names = T)
coloc.info.tables <- lapply(coloc.info.tables, fread)
coloc.info <- rbindlist(coloc.info.tables)

load("results_signal_sharing.RData")
load("results_membership_df.RData")

# list of cell types
all.cell.types <- c("Naive_CD4_T_cells","CM_CD4_T_cells","EM_CD4_T_cells",
                    "Cytotoxic_CD4_T_cells","Regulatory_CD4_T_cells","Naive_CD8_T_cells",
                    "CM_CD8_T_cells","EM_CD8_T_cells","TEMRA","DN_T_cells",
                    "Naive_B_cells","Memory_B_cells","Classical_Monocytes",
                    "CD56Bright_NK_cells","CD56Dim_NK_cells", "all_cells")

# make expression matrix i.e. is each gene expressed in each cell type
expressed.genes <- list()
gene.expression <- read.table("/path/eQTLresults/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Expression_Data.sorted.bed")
expressed.genes[["all_cells"]] <- gene.expression$V4
all_genes <- gene.expression$V4

for(i in all.cell.types){
  if(i == "all_cells"){
    print(i)
  } else {
    print(i)
    gene.expression <- read.table(paste0("/path/eQTLresults/",
                                         i, "/results/TensorQTL_eQTLS/dMean__", i, "_all/OPTIM_pcs/base_output/base/Expression_Data.sorted.bed"))
    gene.expression <- gene.expression$V4
    expressed.genes[[i]] <- gene.expression
  }
}

# list of all expressed genes
all_genes <- unique(unlist(expressed.genes))

# matrix of cell type expression by cell type
gene_matrix <- sapply(expressed.genes, function(genes) all_genes %in% genes)
rownames(gene_matrix) <- all_genes

# signal presence matrix: is each eQTL signal detected in each cell type (from previous script)
# Combine list of shared signals
all_gene_signal_sharing_named <- lapply(all_gene_signal_sharing, function(df) {
  df$rownames <- rownames(df)   # Save rownames as a column
  return(df)
})

# Bind rows, filling missing columns with NA
combined <- bind_rows(all_gene_signal_sharing_named)

# restore rownames
rownames(combined) <- combined$rownames
combined$rownames <- NULL  # Drop the temporary column

# Replace NA with FALSE
signal_matrix <- combined
signal_matrix[is.na(combined)] <- FALSE

# signal info (long form)
all_genes_membership_df <- rbindlist(all_membership_df)
# add short signal name key to membership table
all_genes_membership_df$gene <- unlist(strsplit(all_genes_membership_df$node, "::"))[seq(2, by=3, length.out=nrow(all_genes_membership_df))]
all_genes_membership_df$signal <- paste0(all_genes_membership_df$gene, "_", all_genes_membership_df$group)

# key matching signal name to gene
signal_to_gene <- data.frame("Signal"=rownames(signal_matrix),
                             "Gene"=substr(rownames(signal_matrix), 1, nchar(rownames(signal_matrix))-2))

# for each gene
categorise_gene_pair <- function(gene, cellA, cellB,
                                 expr_mat, signal_matrix,
                                 signal_to_gene, coloc_table,
                                 finemapped.loci) {
  
  # Subset signals for this gene and cell types
  gene_signals <- signal_to_gene$Signal[signal_to_gene$Gene == gene]
  signals_A <- gene_signals[signal_matrix[gene_signals, cellA]]
  signals_B <- gene_signals[signal_matrix[gene_signals, cellB]]
  
  # Category 0: no fine-mapped eGene in those cell types
  if (length(signals_A) == 0 & length(signals_B) == 0) return("0")
  
  # Category 1: expression in only one
  expr_A <- expr_mat[gene, cellA]
  expr_B <- expr_mat[gene, cellB]
  if (xor(expr_A, expr_B)) return("1")
  
  # Category 2: eQTL in only one
  if ((length(signals_A) > 0 && length(signals_B) == 0) ||
      (length(signals_B) > 0 && length(signals_A) == 0)) return("2")
  
  # Category 3
  # signals are not the same
  if(length(intersect(signals_A, signals_B)) == 0) return("3")
  
  #4–5: look at shared signals from coloc
  shared.signals <- intersect(signals_A, signals_B)
  
  shared <- subset(coloc_table, (signal %in% shared.signals) &
                     (celltype == cellA | celltype == cellB))
  
  # Look at effect direction for each signal
  for (s in 1:length(shared.signals)) {
    shared.signal <- shared[shared$signal == shared.signals[s]]
    # check same SNP
    cs.snps <- unique(shared.signal$ol.cs.group)
    # this can be empty if there isn't overlap/coloc between every pair in the group
    # go back to full results and look up just that pair?
    if(is.na(cs.snps)){
      shared.signal.A <- finemapped.loci[grepl(shared.signal$node[shared.signal$celltype == cellA],
                                               finemapped.loci$credible_set_name), ]
      shared.signal.B <- finemapped.loci[grepl(shared.signal$node[shared.signal$celltype == cellB],
                                               finemapped.loci$credible_set_name), ]
      cs.snps.A <- unlist(strsplit(shared.signal.A$credible_set_snps, ","))
      cs.snps.B <- unlist(strsplit(shared.signal.B$credible_set_snps, ","))
      cs.snps <- intersect(cs.snps.A, cs.snps.B)
      
      if(is.na(cs.snps[1])) return("3") 
      
      # get rds files for those SNP
      shared.signal.A <- readRDS(shared.signal.A$path_rds)
      credible.sets.A <- lapply(names(shared.signal.A), function(signal) {
        item <- shared.signal.A[[signal]]
        if (is.list(item) && "finemapping_lABFs" %in% names(item)) {
          df <- item$finemapping_lABFs
          if (is.data.frame(df)) {
            df$signal <- signal  # add the signal name as a new column
            # df <- subset(df, is_cs == TRUE)
            return(df)
          } else {
            NULL
          }
        }
      })
      credible.sets.A <- rbindlist(credible.sets.A)
      credible.sets.A <- subset(credible.sets.A, snp %in% cs.snps)
      
      shared.signal.B <- readRDS(shared.signal.B$path_rds)
      credible.sets.B <- lapply(names(shared.signal.B), function(signal) {
        item <- shared.signal.B[[signal]]
        if (is.list(item) && "finemapping_lABFs" %in% names(item)) {
          df <- item$finemapping_lABFs
          if (is.data.frame(df)) {
            df$signal <- signal  # add the signal name as a new column
            # df <- subset(df, is_cs == TRUE)
            return(df)
          } else {
            NULL
          }
        }
      })
      credible.sets.B <- rbindlist(credible.sets.B)
      credible.sets.B <- subset(credible.sets.B, snp %in% cs.snps)
      
      betaA <- credible.sets.A$bC[credible.sets.A$snp == cs.snps[1]][1]
      betaB <- credible.sets.B$bC[credible.sets.B$snp == cs.snps[1]][1]
      
    } else {
      betaA <- shared.signal$effect[shared.signal$celltype == cellA & shared.signal$ol.cs.group == cs.snps[1]]
      betaA <- betaA[!is.na(betaA)][1]
      betaB <- shared.signal$effect[shared.signal$celltype == cellB & shared.signal$ol.cs.group == cs.snps[1]]
      betaB <- betaB[!is.na(betaB)][1]
    }
    
    if (sign(betaA) == sign(betaB)) return("4")
  }
  
  return("5")  # At least one shared signal with opposing effect
}

# Run function for each gene on this test pair of cells
celltype_pairs <- combn(all.cell.types, 2, simplify = FALSE)

# Prepare output list
results_list <- list()

for (pair in celltype_pairs) {
  cellA <- pair[1]
  cellB <- pair[2]
  
  cat("Processing:", cellA, "vs", cellB, "\n")
  
  pair_results <- sapply(all_genes, function(gene) {
    categorise_gene_pair(gene = gene, cellA = cellA, cellB = cellB,
                         expr_mat = gene_matrix,
                         signal_matrix = signal_matrix,
                         signal_to_gene = signal_to_gene,
                         coloc_table = all_genes_membership_df,
                         finemapped.loci = finemapped.loci)
  })
  
  results_list[[paste(cellA, cellB, sep = "_vs_")]] <- data.frame(
    gene = all_genes,
    celltype_A = cellA,
    celltype_B = cellB,
    category = pair_results
  )
}

# Combine into a single data.frame
all_results <- bind_rows(results_list)
save(all_results, file="eGene_sharing_categorisation.Rdata")

# sanity check
table(all_results$category, useNA = "ifany")
# 0       1       2       3       4       5 
# 1152297   25334  137482    4913   26472     142 

################################################################################
# Make a summary plot 
load("eGene_sharing_categorisation.Rdata")

# remove cases where the gene isn't expressed in either cell type
plot_data_filtered <- all_results %>%
  filter(category != "0") 
table(plot_data_filtered$category, useNA = "ifany")

# calculate proportions of each category for each pairwise comparison
plot_data_porportions <- plot_data_filtered %>%
  count(celltype_A, celltype_B, category) %>%
  group_by(celltype_A, celltype_B) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

# remove bulk-like for plotting
all.cell.types <- c("Naive_CD4_T_cells","CM_CD4_T_cells","EM_CD4_T_cells",
                    "Cytotoxic_CD4_T_cells","Regulatory_CD4_T_cells","Naive_CD8_T_cells",
                    "CM_CD8_T_cells","EM_CD8_T_cells","TEMRA","DN_T_cells",
                    "Naive_B_cells","Memory_B_cells","Classical_Monocytes",
                    "CD56Bright_NK_cells","CD56Dim_NK_cells", "all_cells")

plot_data_upper <- plot_data_porportions %>%
  mutate(
    celltype_A = factor(celltype_A, levels = all.cell.types),
    celltype_B = factor(celltype_B, levels = all.cell.types)
  )

# Name categories
plot_data_upper <- plot_data_upper %>%
  mutate(
    category = factor(category, levels = c("1", "2", "3", "4", "5"),
                      labels = c(
                        "1 - eGene expressed in only one cell type",
                        "2 - eQTL detected in only one cell type",
                        "3 - No overlapping signals between two cell types",
                        "4 - 1+ shared signal with concordant effects",
                        "5 - 1+ shared signal with opposite effects"
                      ))
  )

# # plot results

plot_data_sym <- bind_rows(
  plot_data_upper,
  plot_data_upper %>%
    rename(celltype_A = celltype_B, celltype_B = celltype_A)
)

ordered_celltypes <- c("Naive CD4 T cells","CM CD4 T cells","EM CD4 T cells","Regulatory CD4 T cells","Cytotoxic CD4 T cells","Naive CD8 T cells","CM CD8 T cells","EM CD8 T cells","TEMRA","DN T cells","Naive B cells","Memory B cells","Classical Monocytes","CD56Bright NK cells","CD56Dim NK cells","all cells")

plot_data_upper_ordered <- plot_data_sym %>%
  mutate(
    celltype_A = gsub("_"," ",plot_data_sym$celltype_A),
    celltype_B = gsub("_"," ",plot_data_sym$celltype_B),
    celltype_A = factor(celltype_A, levels = ordered_celltypes),
    celltype_B = factor(celltype_B, levels = ordered_celltypes),
    celltype_A_index = as.integer(celltype_A),
    celltype_B_index = as.integer(celltype_B)
  ) %>%
  filter(celltype_A_index <= celltype_B_index) %>%
  select(-celltype_A_index, -celltype_B_index)

# plot results
ggplot(plot_data_upper_ordered, aes(x = "", y = percentage, fill = category)) +
  geom_col(width = 0.9) +
  facet_grid(rows = vars(celltype_B), cols = vars(celltype_A), switch="both") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  labs(
    title = "eQTL Signal Sharing Between Cell Type Pairs",
    x = "",
    y = "Percentage of Genes",
    fill = "Category"
  ) +
  scale_fill_manual(values = c("#F4978E","#F6E8A6", "#8DD3C7", "#80B1D3","#FFB482"))+
  theme(
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
    #strip.text.y = element_text(angle = 90, size = 8),
    strip.text.x = element_text(angle = 90, size = 10, hjust = 1),
    # aspect.ratio = 1,
    panel.spacing = unit(0.1, "lines")
  )

ggsave("/path/flanders/flanders_grouping.pdf",height=7.5,width=9)
