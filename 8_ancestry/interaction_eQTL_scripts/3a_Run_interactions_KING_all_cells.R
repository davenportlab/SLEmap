# Run interaction analysis for cis-eQTLs from tensorQTL

library(lme4)
library(data.table)
library(glue)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggpubr)
library(jtools)
library(lmtest)
library(dplyr)

symbols_to_names <- read.csv("ensemblID_to_genesymbol.csv")
symbols_to_names$gene_name <- symbols_to_names$gene_symbols


data_dir <- glue("Interactions_coloc_only/all_cells")
dir.create(glue("{data_dir}/output"))
dir.create(glue("{data_dir}/figures"))

figures_dir <- glue("Interactions_coloc_only/Figures/KING")
dir.create(figures_dir)

coloc_results <- read.csv(glue("{data_dir}/coloc_results.csv"))

load(glue("{data_dir}/eqtl_interact_KING_files_pass_hom_filt.rda"))
gex <- t(gex) #samples as rows, genes as columns

lead_eqtl <- sig_res[sig_res$new_rank == 1, ]
lead_eqtl <- merge(lead_eqtl, symbols_to_names[, c("gene_name", "X")], by.x = "phenotype_id", by.y = "X")

sig_res$short_variant <- sub("_[ACGT]+_[ACGT]+$", "", sig_res$variant_id)

if (nrow(pairs.int) == 0){
  next
}

DonorID <- rownames(gex)

pairs.int$gene_name <- gsub(":", ".", pairs.int$gene_name)
pairs.int$gene_name <- gsub(";", ".", pairs.int$gene_name)

irange <- 1:nrow(pairs.int)

covariates$Ancestry <- factor(covariates$Ancestry)
Ancestry <- covariates$Ancestry

results_no_interaction <- do.call(rbind, lapply(irange, function(i){
  Gene = pairs.int[i,1]
  SNP = pairs.int[i,2]
  
  sig_res_gene <- sig_res[sig_res$Gene == Gene, ]
  
  if (nrow(sig_res_gene) == 1){
    model <- lm(gex[,Gene] ~
                  geno[, as.character(SNP)] +
                  as.matrix(covariates[, !names(covariates) %in% c("Ancestry")]) + as.factor(Ancestry))
    
    summary(model)$coefficients[2,]
  } else{
    all_snps <- sig_res_gene$variant_id
    if (SNP %in% all_snps) {
      other_snps <- all_snps[all_snps != SNP]
    }
    else {
      lead_SNP = coloc_results[coloc_results$gene_symbol == Gene & coloc_results$lead_H4_variant_REF_ALT == SNP, ]$lead_snp[1]
      SNP_exclude <- sig_res[sig_res$short_variant == lead_SNP, ]$variant_id[1]
      other_snps <- all_snps[all_snps != SNP_exclude]
    }
    other_geno <- geno[, other_snps]
    model <- lm(gex[,Gene] ~
                  geno[, as.character(SNP)] +
                  other_geno +
                  as.matrix(covariates[, !names(covariates) %in% c("Ancestry")]) +  as.factor(Ancestry))
    
  }
  summary(model)$coefficients[2,]
  
}))

colnames(results_no_interaction) <- c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval")

results_no_interaction <- data.frame(
  Gene = pairs.int[irange, 1],
  SNP = pairs.int[irange, 2],
  results_no_interaction
)

results_no_interaction <- merge(results_no_interaction, qvals[, c("pval_nominal_threshold", "gene_name")], by.x = "Gene", by.y = "gene_name")

results_no_interaction$Still_sig <- results_no_interaction$eQTL_pval < results_no_interaction$pval_nominal_threshold

if (sum(results_no_interaction$Still_sig) < nrow(results_no_interaction)){
  print(glue("For all cells and discrete ancestry : {nrow(results_no_interaction[results_no_interaction$Still_sig == FALSE, ])} out of {nrow(results_no_interaction)} are no longer significant"))
}

write.csv(results_no_interaction, glue("{data_dir}/output/results_no_interaction_with_extra_covariate_KING.csv"))


results_anova <-do.call(rbind, lapply(irange, function(i){
  Gene = pairs.int[i, 1]
  SNP = pairs.int[i, 2]
  
  sig_res_gene <- sig_res[sig_res$Gene == Gene, ]
  
  if (nrow(sig_res_gene) == 1){
    
    model_interact <- lm(gex[,Gene] ~
                           geno[, as.character(SNP)] +
                           as.matrix(covariates[, !names(covariates) %in% c("Ancestry")]) +
                           Ancestry +
                           geno[, as.character(SNP)] * Ancestry)
    
    model_no_interact <- lm(gex[,Gene] ~
                              geno[, as.character(SNP)] +
                              as.matrix(covariates[, !names(covariates) %in% c("Ancestry")]) +
                              Ancestry)
    model_data <- data.frame(
      y = gex[, Gene],
      x = geno[, SNP],
      covariates
    )
    model_data <- na.omit(model_data)
    model_poly <- lm(y ~ poly(x, 2) + ., data = model_data)
    model <- lm(y ~ x + ., data = model_data)
    
    anova_out_poly <- anova(model, model_poly)
    
    anova_out <- anova(model_no_interact, model_interact)
    
    
    c(anova_out$`Pr(>F)`[2],anova_out_poly$`Pr(>F)`[2])
  } else {
    all_snps <- sig_res_gene$variant_id
    if (SNP %in% all_snps) {
      other_snps <- all_snps[all_snps != SNP]
    }
    else {
      lead_SNP = coloc_results[coloc_results$gene_symbol == Gene & coloc_results$lead_H4_variant_REF_ALT == SNP, ]$lead_snp[1]
      SNP_exclude <- sig_res[sig_res$short_variant == lead_SNP, ]$variant_id[1]
      other_snps <- all_snps[all_snps != SNP_exclude]
    }
    other_geno <- geno[, other_snps]
    
    model_interact <- lm(gex[,Gene] ~
                           geno[, as.character(SNP)] +
                           as.matrix(covariates[, !names(covariates) %in% c("Ancestry")]) +
                           Ancestry +
                           other_geno +
                           geno[, as.character(SNP)] * Ancestry)
    
    model_no_interact <- lm(gex[,Gene] ~
                              geno[, as.character(SNP)] +
                              as.matrix(covariates[, !names(covariates) %in% c("Ancestry")]) +
                              other_geno +
                              Ancestry)
    
    model_data <- data.frame(
      y = gex[, Gene],
      x = geno[, SNP],
      covariates,
      other_geno
    )
    model_data <- na.omit(model_data)
    model_poly <- lm(y ~ poly(x, 2) + ., data = model_data)
    model <- lm(y ~ x + ., data = model_data)
    
    anova_out <- anova(model_no_interact, model_interact)
    anova_out_poly <- anova(model, model_poly)
    
    c(anova_out$`Pr(>F)`[2],anova_out_poly$`Pr(>F)`[2])
  }
}))

colnames(results_anova)<-c("anova_pval", "anova_poly_pval")

results_anova <- data.frame(
  Gene = pairs.int[irange, 1],
  SNP = pairs.int[irange, 2],
  results_anova
)

#filter out results where the model with the new covariate term is no longer signficant
results_anova <- merge(results_anova, results_no_interaction[,c("Still_sig", "Gene", "SNP")])
results_anova <- results_anova[results_anova$Still_sig == TRUE, ]


#filter out non-linear main effects
results_anova$Poly_better <- results_anova$anova_poly_pval < 0.05
print(table(results_anova$Poly_better))

results_anova <- results_anova[results_anova$Poly_better == FALSE,] 

results_anova$FDR <- p.adjust(results_anova$anova_pval, method="fdr")
table(results_anova$FDR < 0.05) #1
results_anova$Sig <- results_anova$FDR < 0.05
results_anova_complete <- na.omit(results_anova)
anova_sig <- results_anova_complete[results_anova_complete$Sig == TRUE,] #1

write.csv(results_anova_complete, glue("{data_dir}/output/results_interaction_KING_with_FDR_before_filt.csv"))


if (nrow(anova_sig) > 0){
  anova_sig$lead_eqtl <- interaction(anova_sig$Gene, anova_sig$SNP) %in% interaction(lead_eqtl$gene_name, lead_eqtl$variant_id)
  
  anova_results_list <- list()
  
  for (i in 1: nrow(anova_sig)){
    
    Gene <- anova_sig[i, ]$Gene
    SNP <- anova_sig[i, ]$SNP
    
    sig_res_gene <- sig_res[sig_res$Gene == Gene, ]
    
    if (nrow(sig_res_gene) == 1){
      
      plot_df <- merge(gex[, Gene], geno[, SNP], by=0)
      colnames(plot_df) <- c("IID", "Expression", "Genotype")
      plot_df <- merge(plot_df, covariates, by.x="IID", by.y=0)
    } else {
      all_snps <- sig_res_gene$variant_id
      other_snps <- all_snps[all_snps != SNP]
      
      other_geno <- geno[, other_snps]
      other_geno <- as.data.frame(other_geno)
      other_geno$IID <- rownames(other_geno)
      
      plot_df <- merge(gex[, Gene], geno[, SNP], by=0)
      colnames(plot_df) <- c("IID", "Expression", "Genotype")
      plot_df <- merge(plot_df, other_geno, by = "IID")
      
      plot_df <- merge(plot_df, covariates, by.x="IID", by.y=0)
    }
    plot_df <- na.omit(plot_df)
    
    covariates_model <- colnames(plot_df)[4:(ncol(plot_df))]
    
    rhs <- paste(c("Genotype", "Genotype:Ancestry", covariates_model), collapse = " + ")
    model <- lm(as.formula(paste("Expression ~", rhs)), data = plot_df)
    
    fdr_val <- signif(anova_sig[i, ]$FDR, digits = 4)
    
    gg <- interact_plot(model, pred = Genotype, modx = "Ancestry",
                        plot.points = TRUE, jitter = c(0.1, 0),
                        colors=c("#F8766D", "#00BA38", "#619CFF"),
                        centered = "none",
                        partial.residuals = T, legend.main = paste0("Ancestry"))
    
    p1 <- gg + theme_bw() + ggtitle(glue("{Gene}: {SNP}, FDR: {fdr_val}")) +
      xlab(glue("Genotype")) +
      ylab(glue("Adjusted expression")) +
      scale_x_continuous(breaks=c(0,1,2)) +
      ggtitle(glue("All cells - {Gene}: {SNP}, FDR: {fdr_val}")) + labs(color = "Ancestry")
    
    ggsave(glue("{data_dir}/figures/interaction_{Gene}_{SNP}_all_cells_Ancestry_line.pdf"), p1, width = 6, height = 6)
    ggsave(glue("{figures_dir}/interaction_{Gene}_{SNP}_all_cells_Ancestry_line.pdf"), p1, width = 6, height = 6)
    
    tmp <- ggplot_build(gg)
    tmp <- tmp$data
    tmp <- tmp[[2]]
    tmp$geno <- as.factor(plot_df$Genotype)
    exp.data <- data.frame("Genotype"=tmp$geno,
                           "Genotype_numeric"=plot_df$Genotype,
                           "Exprn"=plot_df$Expression,
                           "AdjustedExprn"=tmp$y,
                           "Ancestry"=plot_df$Ancestry)
    rownames(exp.data) <- plot_df$IID
    
    p2 <- ggplot(exp.data, aes(x = Genotype, y = AdjustedExprn, color = Ancestry)) +
      geom_boxplot(outlier.shape = NA) +
      ylab(glue("Adjusted expression")) +
      geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
      theme_bw() +
      ggtitle(glue("All cells - {Gene}: {SNP}, FDR: {fdr_val}")) + labs(color = "Ancestry")
    p2
    ggsave(glue("{data_dir}/figures/interaction_{Gene}_{SNP}_all_cells_Ancestry_box.pdf"), p2, width = 6, height = 6)
    ggsave(glue("{figures_dir}/interaction_{Gene}_{SNP}_all_cells_Ancestry_box.pdf"), p2, width = 6, height = 6)
    
  }
}


KING_overall_df <- anova_sig

write.csv(KING_overall_df, "Interactions_coloc_only/all_results/all_signif_interactions_discrete_ancestry_all_cells.csv")









