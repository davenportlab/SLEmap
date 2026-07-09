# Summarise all Flanders results

###############################
# Setup
setwd("/path/flanders")

library(UpSetR)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
# library(gridGraphics)
library(tibble)
library(pheatmap)
library(igraph)

###############################

all.cell.types <- c("Naive_CD4_T_cells","CM_CD4_T_cells","EM_CD4_T_cells",
                    "Cytotoxic_CD4_T_cells","Regulatory_CD4_T_cells","Naive_CD8_T_cells",
                    "CM_CD8_T_cells","EM_CD8_T_cells","TEMRA","DN_T_cells",
                    "Naive_B_cells","Memory_B_cells","Classical_Monocytes",
                    "CD56Bright_NK_cells","CD56Dim_NK_cells", "all_cells")

# read in Flanders results files
loci <- list.files("flanders_updated_outputs/results/coloc_info_tables/", full.names = T)
finemapped.loci <- lapply(loci, fread)
finemapped.loci <- rbindlist(finemapped.loci)

coloc.results.old <- read.delim("flanders_outputs_test/results/coloc/coloc_run_colocalization.table.all.tsv")
coloc.results <- read.delim("flanders_updated_outputs/results/coloc/coloc_run_colocalization.table.all.tsv")

# this looks the same as finemapped loci
coloc.info.tables <- list.files("flanders_updated_outputs/results/coloc_info_tables/", full.names = T)
coloc.info.tables <- lapply(coloc.info.tables, fread)
coloc.info <- rbindlist(coloc.info.tables)

# read in TensorQTL results for comparison
input.info <- read.delim("slemap_flanders_input_updated.txt")
original.results <- list()
original.cond.results <- list()

for(i in 1:15){
  cell <- input.info$study_id[i]
  original.file <- input.info$input[i]
  original.results.i <- fread(original.file)
  original.results.i$celltype <- cell
  original.results[[i]] <- original.results.i
  
  cond.file <- paste0("/path/eQTLresults/",
                      cell, "/results/TensorQTL_eQTLS/dMean__",
                      cell, "_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv")
  original.cond.results.i <- fread(cond.file)
  original.cond.results.i$celltype <- cell
  original.cond.results[[i]] <- original.cond.results.i
}

original.results.i <- fread(input.info$input[16])
original.results.i$celltype <- "all_cells"
original.results[[16]] <- original.results.i

original.cond.results.i <- fread("/path/eQTLresults/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv")
original.cond.results.i$celltype <- "all_cells"
original.cond.results[[16]] <- original.cond.results.i

original.results <- rbindlist(original.results)
original.cond.results <- rbindlist(original.cond.results)

# set up lists to store per-gene outputs
genes <- unique(finemapped.loci$phenotype_id)
all_gene_signal_sharing <- list()
all_membership_df <- list()

blank <- nullGrob()

for(gi in 1:length(genes)) {
  
  gene <- genes[gi]

  pdf(paste0("/path/flanders/Per_gene_flanders_plots/", gene, ".pdf"),
      useDingbats = F, onefile = T)
  
  gene.finemapped.loci <- subset(finemapped.loci, phenotype_id == gene)
  gene.finemapped.loci$snppos <- as.numeric(unlist(strsplit(gene.finemapped.loci$snp, ":"))[seq(from=2,
                                                                                                by=4,
                                                                                                length.out=nrow(gene.finemapped.loci))])
  # plot original TensorQTL results with identified conditional signals highlighted
  gene.original <- subset(original.results, phenotype_id == gene)
  gene.cond.original <- subset(original.cond.results, phenotype_id == gene)
  
  # Create a unique key to match on
  gene.original$key <- paste(gene.original$variant_id, gene.original$celltype, sep = "_")
  gene.cond.original$key <- paste(gene.cond.original$variant_id, gene.cond.original$celltype, sep = "_")
  
  # Mark significant SNPs
  gene.original$cond_lead <- gene.original$key %in% gene.cond.original$key
  
  g0 <- ggplot(gene.original, aes(pos, -log10(pval_nominal))) +
    geom_point(aes(colour=cond_lead)) +
    scale_colour_manual(values=c("lightgrey", "darkblue")) +
    facet_wrap(~ celltype) +
    theme_bw()
  print(g0)
  
  # plot lead SNP positions from SuSiE
  g1 <- ggplot(gene.finemapped.loci, aes(snppos, -log10(top_pvalue))) +
    geom_point(aes(colour=study_id)) +
    theme_bw()
  
  grid.arrange(g1, blank, heights = c(3, 1))
  
  # get credible set for each: number of files doesn't equal number of signals
  # 1 file might have 1 signal, another 3 - need to check rank across them all
  cs.files <- unique(gene.finemapped.loci$path_rds)
  
  credible.sets <- lapply(cs.files, function(cs.file) {
    rds_obj <- readRDS(cs.file[1])
    finemapping <- lapply(names(rds_obj), function(signal) {
      item <- rds_obj[[signal]]
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
    finemapping <- rbindlist(finemapping)
  })
  
  credible.sets <- rbindlist(credible.sets)
  credible.sets$celltype <- unlist(strsplit(credible.sets$signal, "::"))[seq(from=2, by=5, length.out=nrow(credible.sets))]
  
  # for each cell type get rank (based on order of signal)
  credible.sets <- credible.sets[order(credible.sets$lABF, decreasing = T), ]
  credible.sets <- credible.sets %>%
    group_by(celltype) %>%
    mutate(
      rank = match(signal, unique(signal))
    ) %>%
    ungroup()
  
  credible.sets$is_cs[credible.sets$is_cs == TRUE] <- paste0(credible.sets$is_cs[credible.sets$is_cs == TRUE],
                                                             credible.sets$rank[credible.sets$is_cs == TRUE])
  
  g2 <- ggplot(credible.sets, aes(position, lABF)) +
    geom_point(aes(colour=is_cs)) +
    scale_colour_manual(values=c("lightgrey", "dodgerblue", "darkgreen", "red", "gold", "darkblue", "pink", "purple", "orange", "lightgreen", "lavender", "lightblue", "black")) +
    facet_wrap(~ celltype) +
    theme_bw()
  print(g2)
  
  if(nrow(gene.finemapped.loci) > 1){
    g3 <- ggplot(subset(credible.sets, is_cs != FALSE), aes(position, lABF)) +
      geom_point(aes(colour=signal)) +
      theme_bw() +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 8, byrow = TRUE))
    
    grid.arrange(g3)
  }
  
  # Check for overlap in CS
  credible.sets.full <- credible.sets
  credible.sets <- subset(credible.sets, is_cs != FALSE)
  
  # Initialize empty edge list
  edges <- c()
  
  # Get all unique nodes
  signal_names <- unique(credible.sets$signal)
  
  for (i in 1:(length(signal_names) - 1)) {
    for (j in (i + 1):length(signal_names)) {
      sig1 <- signal_names[i]
      sig2 <- signal_names[j]
      
      # Check if they share any variants
      overlap <- length(intersect(credible.sets$snp[credible.sets$signal == sig1], 
                                  credible.sets$snp[credible.sets$signal == sig2])) > 0
      
      if (overlap) {
        edges <- rbind(edges, c(sig1, sig2))
      }
    }
  }
  
  if(is.null(edges)){
    g <- make_empty_graph(n = length(signal_names), directed = FALSE)
    # Add names
    V(g)$name <- signal_names
  } else {
    # Create graph object
    g <- graph_from_data_frame(as.data.frame(edges), 
                               directed = FALSE,
                               vertices=signal_names)
  }
  
  # Plot with basic layout
  plot(g, vertex.label.cex = 0.7, vertex.size = 5, main = "Signal Overlap via Credible Sets")
  
  signal_names_short <- sub("^chr[^:]+::", "", signal_names)  # Remove the first "chrX::"
  signal_names_short <- sub("::L[0-9]+$", "", signal_names_short)     # Remove trailing "::L<number>"
  
  if(!is.null(edges) & nrow(gene.finemapped.loci) > 1){
    # then see if tested for coloc and what the results were
    gene.coloc.results <- subset(coloc.results, t1 == gene & t2 == gene)
    
    # get list of signals that can be merged based on coloc results
    gene.coloc.results.simple <- data.frame(
      "Signal_1"=sub("::L[0-9]+$", "", sub("^chr[^:]+::", "", gene.coloc.results$hit1)),
      "Signal_2"=sub("::L[0-9]+$", "", sub("^chr[^:]+::", "", gene.coloc.results$hit2)),
      "H4"=gene.coloc.results$PP.H4.abf,
      "Coloc"=gene.coloc.results$PP.H4.abf > 0.8)
    
    # Create a graph object from the data
    graph <- graph_from_data_frame(gene.coloc.results.simple[, 1:2], directed = FALSE,
                                   vertices=signal_names_short) # Use first two columns for edges
    
    # Add edge attributes based on the 'coloc' column
    E(graph)$connected <- gene.coloc.results.simple$Coloc
    
    plot(graph,
         edge.color = ifelse(E(graph)$connected, "black", "white"), # Different colors for true/false
         edge.lty = ifelse(E(graph)$connected, 1, 0), # Different line types for true/false
         vertex.label.color = "black",
         vertex.label.cex = 0.7, vertex.size = 5,
         main = "Coloc results"
    )
    
  }
  
  # Make summary table
  gene.coloc.info <- subset(coloc.info, phenotype_id == gene)
  
  # Start with all unique signals
  all_nodes <- unique(sub("::L[0-9]+$", "", sub("^chr[^:]+::", "", gene.coloc.info$credible_set_name)))
  
  # Step 2: Create graph with only connected edges, but include all nodes
  if(is.null(edges) | nrow(gene.finemapped.loci) == 1){
    g <- make_empty_graph(n = length(all_nodes), directed = FALSE)
    # Add names
    V(g)$name <- all_nodes
  } else {
    connected_edges <- gene.coloc.results.simple %>% filter(Coloc)
    g <- graph_from_data_frame(connected_edges, directed = FALSE, vertices = all_nodes)
  }
  
  # Step 3: Get connected components
  components <- components(g)
  
  # Step 4: Build data frame for group membership
  membership_df <- data.frame(
    node = names(components$membership),
    group = components$membership
  )
  
  membership_df$celltype <- unlist(strsplit(membership_df$node, "::"))[seq(from=1, by=3, length.out=nrow(membership_df))]
  
  # if any cell type has multiple signals in the same group we need to merge them 
  # (just remove the second one at this point but get average for effect size?)
  membership_df$check_dups <- paste0(membership_df$group, membership_df$celltype)
  membership_df <- membership_df[!duplicated(membership_df$check_dups), ]
  membership_df$check_dups <- NULL
  
  # Step 5: Wide format matrix
  group_matrix <- membership_df[, 2:3] %>%
    mutate(present = TRUE) %>%
    pivot_wider(names_from = celltype, values_from = present, values_fill = FALSE) %>%
    arrange(group) %>%
    select(-group)  # remove if you don't want group IDs
  
  gene_signal_sharing <- data.frame(group_matrix)
  
  # Name signals as gene and number
  rownames(gene_signal_sharing) <- paste0(gene, "_", 1:nrow(gene_signal_sharing)) # ordered by group
  all_gene_signal_sharing[[gi]] <- gene_signal_sharing
  
  
  # TO DO
  # Should I use the reduced membership df going forward?
  # Add signal rank per cell type
  
  # Alternative step 5: Wide format matrix
  group_matrix_signals <- membership_df %>%
    pivot_wider(names_from = celltype, values_from = node, values_fn = ~paste(.x, collapse = ", ")) %>%
    arrange(group)
  
  gene_signals <- data.frame(group_matrix_signals)
  
  # Name signals as gene and number
  rownames(gene_signals) <- paste0(gene, "_", 1:nrow(gene_signals))
  
  # add effect size
  gene.finemapped.loci$red.name <- paste0(gene.finemapped.loci$study_id, "::",
                                          gene.finemapped.loci$phenotype_id, "::",
                                          gene.finemapped.loci$snp)
  
  # what are the overlapping SNPs in each group (if there is any overlap) and the signal ranks?
  credible.sets <- credible.sets[order(credible.sets$lABF, decreasing=TRUE), ]
  credible.sets$signal_names_short <- sub("^chr[^:]+::", "", credible.sets$signal)  # Remove the first "chrX::"
  credible.sets$signal_names_short <- sub("::L[0-9]+$", "", credible.sets$signal_names_short)     # Remove trailing "::L<number>"
  
  membership_df$percellrank <- credible.sets$rank[match(membership_df$node,
                                                 credible.sets$signal_names_short)]
  
  membership_df$ol.cs.group <- NA
  
  for(s in 1:nrow(gene_signals)){
    signal_names <- membership_df$node[which(membership_df$group == s)]
    
    cs.snps <- credible.sets$snp[grepl(signal_names[1], credible.sets$signal)]
    
    # if more than 1 signal
    if(length(signal_names) > 1){
      for(i in 2:(length(signal_names))){
        cs.snps <- intersect(cs.snps,
                             credible.sets$snp[grepl(signal_names[i], credible.sets$signal)])
      }
    }
    
    # If there is no shared CS SNP, get frequent one to represent signal
    # need to check that it was tested in the other signals
    if(length(cs.snps) == 0){
        tested.snps <- credible.sets.full$snp[grepl(signal_names[1], credible.sets.full$signal)]
        for(i in 2:(length(signal_names))){
          tested.snps <- intersect(tested.snps,
                                   credible.sets.full$snp[grepl(signal_names[i], credible.sets.full$signal)])
        }
        cs.snps <- credible.sets$snp[credible.sets$signal_names_short %in% signal_names]
        cs.snps <- cs.snps[cs.snps %in% tested.snps]
        cs.snps <- names(table(cs.snps))[as.vector(table(cs.snps))==max(table(cs.snps))]  
      }
    
    membership_df$ol.cs.group[which(membership_df$group == s)] <- cs.snps[1]
    
    scale.limits.neg <- min(-0.1, min(credible.sets[which(credible.sets$signal_names_short %in% signal_names), "bC"]))
    scale.limits.pos <- max(0.1, max(credible.sets[which(credible.sets$signal_names_short %in% signal_names), "bC"]))
    
    scale.limits <- c((max(abs(scale.limits.neg), scale.limits.pos)*-1),
                      max(abs(scale.limits.neg), scale.limits.pos))
    
    gg <- ggplot(credible.sets[which(credible.sets$signal_names_short %in% signal_names), ], 
                 aes(celltype, snp, fill=bC)) +
      geom_tile() +
      theme_bw() + 
      scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, limits=scale.limits) +
      ggtitle(rownames(gene_signals)[s])
    print(gg)
  }
  
  # summarise direction of effect
  membership_df$effect <- NA
  membership_df$isCS <- TRUE
  membership_df_g_list <- list()
  
  for(g in 1:length(unique(membership_df$group))){
    membership_df_g <- membership_df[membership_df$group == g, ]
    
    # if there is only 1 signal in group) just look up in fine mapping
    if(nrow(membership_df_g) == 1){
      membership_df_g$effect <- gene.finemapped.loci$bC[match(membership_df_g$node,
                                                              gene.finemapped.loci$red.name)]
    }
    
    # otherwise look up CS
    for(i in 1:nrow(membership_df_g)){
      # if there is a CS SNP listed (shared across all conditions or frequent) use that
      membership_df_g$effect[i] <- credible.sets$bC[(grepl(membership_df_g$node[i], credible.sets$signal) &
                                                       grepl(membership_df_g$ol.cs.group[i], credible.sets$snp))][1]
      # if CS SNP is not shared across all, some will be NA - this is good - should look up but note that it's not in CS
      # if it's not in the credible set for that signal need to look in full outputs
      if(is.na(membership_df_g$effect[i])){
        membership_df_g$isCS[i] <- FALSE
        membership_df_g$effect[i] <- credible.sets.full$bC[(grepl(membership_df_g$node[i], credible.sets.full$signal) &
                                                              grepl(membership_df_g$ol.cs.group[i], credible.sets.full$snp))][1]
      }
    }
    membership_df_g_list[[g]] <- membership_df_g
  }
  
  membership_df_g_list_named <- lapply(membership_df_g_list, function(df) {
    df$rownames <- rownames(df)   # Save rownames as a column
    return(df)
  })
  
  # Bind rows
  membership_df <- bind_rows(membership_df_g_list_named)

  # restore rownames
  rownames(membership_df) <- membership_df$rownames
  membership_df$rownames <- NULL  # Drop the temporary column

  # TO DO
  # If there are multiple signals being combined (separated by finemapping, merged by coloc), 
  # check direction is same/set to NA?
  
  # Save original column order
  original_col_order <- colnames(membership_df)
  
  collapsed_df <- membership_df %>%
    rownames_to_column("rowname") %>%
    group_by(group, celltype) %>%
    # some signals can be NA, hopefully the first one isn't
    summarise(effect = {
        effect_no_na <- effect[!is.na(effect)]
        if (length(effect_no_na) == 0) {
          NA_real_
        } else if (all(sign(effect_no_na) == sign(effect_no_na[1]))) {
          mean(effect_no_na)
        } else {
          NA_real_
        }
      },
      isCS = any(isCS),  # custom rule for logical column
      across(
        .cols = -c(effect, isCS, rowname),
        .fns = ~ if (length(unique(.x)) == 1) unique(.x) else paste(unique(.x), collapse = ", ")), 
      rowname = paste(rowname, collapse = ", "),
      .groups = "drop") %>%
    column_to_rownames("rowname") %>%
    select(any_of(original_col_order))
  
  collapsed_df$celltype <- factor(x = collapsed_df$celltype,
                                  levels=all.cell.types)
  
  scale.limits.neg <- min(0.5, min(collapsed_df$effect))
  scale.limits.pos <- max(0.5, max(collapsed_df$effect))
  
  scale.limits <- c((max(abs(scale.limits.neg), scale.limits.pos)*-1),
                    max(abs(scale.limits.neg), scale.limits.pos))
  
  collapsed_df$isCS <- factor(collapsed_df$isCS, levels=c(TRUE, FALSE))
  
  g4 <- ggplot(collapsed_df, aes(celltype, as.factor(group))) +
    geom_tile(aes(fill=effect)) +
    geom_tile(data = subset(collapsed_df, isCS == TRUE), fill = NA, color = "grey", linewidth = 1, linetype=1) +
    geom_text(aes(label = percellrank), size = 3) +
    theme_bw() + ylab("Signal") +
    scale_x_discrete(drop=FALSE) +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
                         limits=scale.limits) +
    ggtitle(gene) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # work out how much of the plot to occupy based on n signals
  plot.size <- min(nrow(gene_signals), 3)
  grid.arrange(g4, blank, heights = c(plot.size, (3-plot.size)))
  
  # replot credible sets but numbered by shared signals
  credible.sets.full$signal_names_short <- sub("^chr[^:]+::", "", credible.sets.full$signal)  # Remove the first "chrX::"
  credible.sets.full$signal_names_short <- sub("::L[0-9]+$", "", credible.sets.full$signal_names_short)     # Remove trailing "::L<number>"
  
  credible.sets.full$shared_signal <- collapsed_df$group[match(credible.sets.full$signal_names_short,
                                                          collapsed_df$node)]
  credible.sets.full$shared_signal[which(credible.sets.full$is_cs == FALSE)] <- 0
  g5 <- ggplot(credible.sets.full, aes(position, lABF)) +
    geom_point(aes(colour=as.factor(shared_signal))) +
    scale_colour_manual(values=c("lightgrey", "dodgerblue", "darkgreen", "red", "gold", "darkblue", "pink", "purple", "orange", "lightgreen", "lavender", "lightblue", "black")) +
    facet_wrap(~ celltype) +
    theme_bw()
  print(g5)
  
  all_membership_df[[gi]] <- membership_df
  
  dev.off()
  

}

save(all_membership_df, file = "results_membership_df_ranks.RData")
