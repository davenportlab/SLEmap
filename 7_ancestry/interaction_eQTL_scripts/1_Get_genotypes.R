library(ggplot2)
library(data.table)
library(glue)
library(dplyr)

celltypes_to_test <- c("Classical_Monocytes", "DN_T_cells", "Memory_B_cells", "Regulatory_CD4_T_cells", 
                       "CM_CD4_T_cells", "EM_CD4_T_cells", "Naive_B_cells", "TEMRA", "CD56Bright_NK_cells",
                       "CM_CD8_T_cells", "EM_CD8_T_cells", "Naive_CD4_T_cells", "CD56Dim_NK_cells",
                       "Cytotoxic_CD4_T_cells", "Naive_CD8_T_cells")

coloc_results <- read.csv("coloc_sigresults_GCST90270940_checksigeQTL.csv")
coloc_results$lead_H4_variant_REF_ALT <- paste(
  coloc_results$lead_H4_variant,
  coloc_results$lead_H4_variant_GTF_REF,
  coloc_results$lead_H4_variant_GTF_ALT,
  sep = "_"
)

for (celltype in celltypes_to_test){
  data_dir <- glue("Interactions_coloc_only/{celltype}")
  dir.create(file.path(data_dir))
  genotypes_dir <- glue("{data_dir}/genotypes")
  dir.create(file.path(genotypes_dir))
  coloc_results_celltype <- coloc_results[coloc_results$cell_type == celltype,]
  
  write.csv(coloc_results_celltype, glue("Interactions_coloc_only/{celltype}/coloc_results.csv"), row.names = FALSE)
  
  wgs_dir <- "WGS_data"
  
  # QTLight output directory for celltype
  qtlight_dir <- glue("5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/{celltype}")
  
  # Get normalised and pseudobulked gene expression data
  gex <- read.delim(glue("{qtlight_dir}/results/norm_data/dMean__{celltype}_all/normalised_phenotype.tsv"))
  #genes as rows, samples as cols
  keep_samples <- colnames(gex)
  keep_samples <- gsub(".*(SLE_WGS\\d+).*", "\\1", keep_samples)
  keep_samples <- gsub(".*(SLE_map\\d+).*", "\\1", keep_samples)
  colnames(gex) <- keep_samples
  
  write.table(keep_samples, glue("{data_dir}/tested_samples.txt"), quote = F, row.names = F, col.names = F)
  
  # Get genotype and phenotype PCs
  covariates <- read.delim(glue("{qtlight_dir}/results/TensorQTL_eQTLS/dMean__{celltype}_all/OPTIM_pcs/base_output/base/Covariates.tsv"), row.names = 1)
  covariates <- t(covariates) # covariates as cols, samples as rows
  colnames(covariates) <- gsub(" ", "_", colnames(covariates))
  covariates <- data.frame(covariates)
  covariates <- covariates[keep_samples, ]
  
  ggplot(data = covariates, aes(x = Genotype_PC1, y = Genotype_PC2)) + geom_point()
  
  
  # Get conditional eQTL results
  res <- read.delim(glue("{qtlight_dir}/results/TensorQTL_eQTLS/dMean__{celltype}_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"))
  sig_res <- res
  snps.sig <- unique(sig_res$variant_id)
  num_sig_snp_pairs <- nrow(snps.sig)
  
  lead_coloc_snps <- unique(coloc_results_celltype$lead_H4_variant_REF_ALT)
  
  all_snps_to_use <- unique(unlist(c(snps.sig, lead_coloc_snps)))
  
  write.table(all_snps_to_use, glue("{genotypes_dir}/conditional_eqtl_sig_snps_and_lead_coloc_snps.txt"), quote = F, row.names = F, col.names = F)
  
  preprocess_int <- glue("{data_dir}/preprocess_files.rda")
  save(list=c("keep_samples", "covariates", "lead_coloc_snps", "sig_res", "gex"), file = preprocess_int)
}

# Run plink to get geno info for main eqtl snps

#!/bin/bash

# module load HGI/softpack/users/wl2/data_QC/2
# 
# celltypes_to_test=(
#   "Classical_Monocytes" "DN_T_cells" "Memory_B_cells" "Regulatory_CD4_T_cells"
#   "CM_CD4_T_cells" "EM_CD4_T_cells" "Naive_B_cells" "TEMRA" "CD56Bright_NK_cells"
#   "CM_CD8_T_cells" "EM_CD8_T_cells" "Naive_CD4_T_cells" "CD56Dim_NK_cells"
#   "Cytotoxic_CD4_T_cells" "Naive_CD8_T_cells" "DN_T_cells"
# )
# 
# for celltype in "${celltypes_to_test[@]}"; do
# 
# plink2 \
# --vcf 281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID.vcf.gz \
# --extract Interactions_coloc_only/${celltype}/genotypes/conditional_eqtl_sig_snps_and_lead_coloc_snps.txt \
# --recode A --out Interactions_coloc_only/${celltype}/genotypes/sig_snps \
# --make-bed
# 
# done