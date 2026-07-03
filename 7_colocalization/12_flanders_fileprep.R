# SLEmap: Run Flanders to compare eQTL across cell types

library(data.table)
library(dplyr)
library(pbapply)

setwd("/path/flanders/")

# make the input parameters file for Flanders
columns.needed <- c("input", "study_id", "chr_lab", "pos_lab", "rsid_lab", "a1_lab",
                    "a0_lab", "freq_lab", "n_lab", "effect_lab", "se_lab",
                    "pvalue_lab", "type", "sdY", "s", "grch", "p_thresh1", "p_thresh2",
                    "hole", "bfile", "maf", "is_molQTL", "key", "cs_thresh", "grch_bfile",
                    # new ones
                    "process_bfile",	"per_gene_p1",	"per_gene_p2")

input.info <- matrix(NA, ncol=28, nrow=16, dimnames = list(1:16, columns.needed))

# All eQTL
input.info[, "type"] <- "quant"
input.info[, "is_molQTL"] <- 	TRUE

# Project-specific settings
input.info[, "grch"] <- 38
input.info[, "hole"] <- 	200000
input.info[, "maf"] <- 	0.05
input.info[, "cs_thresh"] <- 	0.95 # reduce to see if get fewer empty?
input.info[, "grch_bfile"] <- 38
input.info[, "process_bfile"] <- TRUE

input.info[, "p_thresh1"] <- 	1.00E-05 # to be updated per cell type and overridden by gene-specific thresholds
input.info[, "p_thresh2"] <- 	1.00E-03 # to be updated per cell type and overridden by gene-specific thresholds

# Existing columns in the tensorQTL output
input.info[, "rsid_lab"] <- "variant_id"
input.info[, "freq_lab"] <- "af"
input.info[, "effect_lab"] <- "slope"
input.info[, "se_lab"] <- "slope_se"
input.info[, "pvalue_lab"] <- "pval_nominal"
input.info[, "key"] <- 	"phenotype_id"

# Others we need to add to each results file
input.info[, "chr_lab"] <- "chr"
input.info[, "pos_lab"] <- "pos"
input.info[, "a1_lab"] <- "a1"
input.info[, "a0_lab"] <- "a0"
input.info[, "n_lab"] <- "nsamples"

# sc-eQTL results from 2 QTLight runs (by cell type, and all cells together): 

# get list of cell types
cell.results <- list.files(path="/path/mashr/input/",
                           pattern="*tsv.gz", recursive = T, full.names = T)

# For each file, read in, add additional columns and write out to flanders directory
for(i in 1:length(cell.results)){
  
  # Get cell type name from directory name and add label to flanders input file
  nominal.file.name <- cell.results[i]
  cell.name <- gsub("/path/mashr/input//merged_",
                    "", nominal.file.name)
  cell.name <- gsub(".tsv.gz", "", cell.name)
  print(cell.name)
  input.info[i, "study_id"] <- cell.name
  
  # pvalue thresholds
  qval.file.name <- paste0("/path/eQTLresults/",
                           cell.name, "/results/TensorQTL_eQTLS/dMean__", cell.name, "_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv")

  # Load significant genes (qval < 0.05)
  sig_genes <- fread(qval.file.name) %>%
    filter(qval < 0.05)
  
  # find largest significant pvalue for flat threshold
  input.info[i, "p_thresh1"] <- max(sig_genes$pval_nominal_threshold)
  # use psig threshold to get p-boundary threshold?
  input.info[i, "p_thresh2"] <- min(c(1000*(max(sig_genes$pval_nominal_threshold)), 0.001))

  # make per-gene threshold files
  pthresh1 <- sig_genes[, c("phenotype_id", "pval_nominal_threshold")]
  colnames(pthresh1) <- c("phenotype_id", "pval_thresh")
  write.table(pthresh1, paste0(cell.name, "_eqtl_pval_thresh1_per_gene.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # p-thresh2 will come from p-thresh1 with caps
  # For each gene, the variant inclusion threshold (p-thresh2) for Flanders fine-mapping was scaled from the gene-specific significance threshold (p-thresh1) estimated by eigenMT, reflecting the effective number of independent tests. We used a piecewise scaling with an upper cap: p-thresh2 = min(10 × p1, 1×10⁻⁴) for p1 < 1×10⁻⁴, p-thresh2 = min(10 × p1, 1×10⁻³) for 1×10⁻⁴ ≤ p1 < 1×10⁻³, and p-thresh2 = min(5 × p1, 1×10⁻²) for p1 ≥ 1×10⁻³. This approach allows inclusion of variants in LD with lead eQTLs while controlling region size, preserving stable fine-mapping performance.
  pthresh2 <- pthresh1 %>% mutate(p_thresh2 = ifelse(pval_thresh < 1e-4, pmin(10*pval_thresh, 1e-4), 
                              ifelse(pval_thresh < 1e-3, pmin(10*pval_thresh, 1e-3),
                                     pmin(5*pval_thresh, 1e-2)))) %>% select(-pval_thresh)
  colnames(pthresh2) <- c("phenotype_id", "pval_thresh")
  write.table(pthresh2, paste0(cell.name, "_eqtl_pval_thresh2_per_gene.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  input.info[i, "per_gene_p1"] <- paste0(cell.name, "_eqtl_pval_thresh1_per_gene.tsv")
  input.info[i, "per_gene_p2"] <- paste0(cell.name, "_eqtl_pval_thresh2_per_gene.tsv")

  
  ##### sdY ##### 
  gene_exp <- read.delim(paste0("/path/eQTLresults/",
                             cell.name, "/results/norm_data/dMean__", cell.name, "_all/dMean__", cell.name, "___phenotype_file.tsv"),
                         row.names = 1)
  
  # calculate sd for each gene 
  gene_sd <- apply(gene_exp, 1, sd)
  # Make a data frame with gene IDs and SDs
  gene_sd_df <- data.frame(
    phenotype_id = names(gene_sd),
    sdY = gene_sd
  )
  # write out all sds in the correct format: TSV table with columns phenotype_id and sdY
  write.table(gene_sd_df, paste0(cell.name, "_sdY_per_gene.tsv"), sep = "\t", 
              quote = FALSE, row.names = FALSE)
  
  input.info[i, "sdY"] <- paste0(cell.name, "_sdY_per_gene.tsv")

  # Add info to results files
  # Load nominal p-values
  nominal <- fread(nominal.file.name)  # columns: phenotype_id, variant_id, tss_distance, pval_nominal, slope, etc.
  
  # Filter for significant genes
  results <- nominal %>% filter(phenotype_id %in% sig_genes$phenotype_id)
  rm(nominal)

  # calculate number of samples (1/2 number of alleles) from the AF and MA count
  # MAF=MA count/(number of samples*2)
  # Number of samples = MA count/(MAF*2)
  results$maf <- results$af
  results$maf[which(results$maf > 0.5)] <- 1- results$maf[which(results$maf > 0.5)]
  results$nsamples <- round(results$ma_count/(results$maf*2))
  results$maf <- NULL
  
  # Split up variant info
  variant.info <- unlist(strsplit(results$variant_id, "_"))
  results$chr <- gsub("X", "23", gsub("chr", "", variant.info[seq(1, length(variant.info), by=4)]))
  results$pos <- variant.info[seq(2, length(variant.info), by=4)]
  results$a0 <- variant.info[seq(3, length(variant.info), by=4)]
  results$a1 <- variant.info[seq(4, length(variant.info), by=4)]
  rm(variant.info)
  
  # write out updated results file and put path in flanders input file
  input.info[i, "input"] <- paste0("/path/flanders/", cell.name, ".txt")
  write.table(results,
              paste0("/path/flanders/", cell.name, ".txt"), sep="\t",
              quote=F, row.names=F)
}

# add bulk
input.info[16, "study_id"] <- "all_cells"

input.info[, "bfile"] <- "/path/flanders/plink_genotypes"

# read in eQTL results file
all.results <- "/path/flanders/merged_allcells.tsv.gz"
results <- fread(all.results)

# Load significant genes (qval < 0.05)
qval.file.name <- paste0("/path/eQTLresults/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv")
sig_genes <- fread(qval.file.name) %>%
  filter(qval < 0.05)

# find largest significant pvalue for flat thresholds
input.info[16, "p_thresh1"] <- max(sig_genes$pval_nominal_threshold)
input.info[16, "p_thresh2"] <- min(c(1000*(max(sig_genes$pval_nominal_threshold)), 0.001))

range(input.info[, "p_thresh1"])
range(input.info[, "p_thresh2"])

# make per-gene threshold files
pthresh1 <- sig_genes[, c("phenotype_id", "pval_nominal_threshold")]
colnames(pthresh1) <- c("phenotype_id", "pval_thresh")
write.table(pthresh1, "all_cells_eqtl_pval_thresh1_per_gene.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# p-thresh2 will come from p-thresh1 with caps
# For each gene, the variant inclusion threshold (p-thresh2) for Flanders fine-mapping was scaled from the gene-specific significance threshold (p-thresh1) estimated by eigenMT, reflecting the effective number of independent tests. We used a piecewise scaling with an upper cap: p-thresh2 = min(10 × p1, 1×10⁻⁴) for p1 < 1×10⁻⁴, p-thresh2 = min(10 × p1, 1×10⁻³) for 1×10⁻⁴ ≤ p1 < 1×10⁻³, and p-thresh2 = min(5 × p1, 1×10⁻²) for p1 ≥ 1×10⁻³. This approach allows inclusion of variants in LD with lead eQTLs while controlling region size, preserving stable fine-mapping performance.
pthresh2 <- pthresh1 %>% mutate(p_thresh2 = ifelse(pval_thresh < 1e-4, pmin(10*pval_thresh, 1e-4), 
                                                   ifelse(pval_thresh < 1e-3, pmin(10*pval_thresh, 1e-3),
                                                          pmin(5*pval_thresh, 1e-2)))) %>% select(-pval_thresh)
colnames(pthresh2) <- c("phenotype_id", "pval_thresh")
write.table(pthresh2, "all_cells_eqtl_pval_thresh2_per_gene.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

input.info[16, "per_gene_p1"] <- "all_cells_eqtl_pval_thresh1_per_gene.tsv"
input.info[16, "per_gene_p2"] <- "all_cells_eqtl_pval_thresh2_per_gene.tsv"

##### sdY
gene_exp <- read.delim("/path/eQTLresults/results/norm_data/dMean__All_all/dMean__All___phenotype_file.tsv",
                              row.names = 1)

# calculate sd for each gene 
gene_sd <- apply(gene_exp, 1, sd)
# Make a data frame with gene IDs and SDs
gene_sd_df <- data.frame(
  phenotype_id = names(gene_sd),
  sdY = gene_sd
)
# write out all sds in the correct format: TSV table with columns phenotype_id and sdY
write.table(gene_sd_df, "all_cells_sdY_per_gene.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)

input.info[16, "sdY"] <- "all_cells_sdY_per_gene.tsv"

# Filter for significant genes
results <- results %>% filter(phenotype_id %in% sig_genes$phenotype_id)

# calculate number of samples (1/2 number of alleles) from the AF and MA count
# MAF=MA count/(number of samples*2)
# Number of samples = MA count/(MAF*2)
results$maf <- results$af
results$maf[which(results$maf > 0.5)] <- 1- results$maf[which(results$maf > 0.5)]
results$nsamples <- round(results$ma_count/(results$maf*2))
results$maf <- NULL
  
# Split up variant info
variant.info <- unlist(strsplit(results$variant_id, "_"))
results$chr <- gsub("X", "23", gsub("chr", "", variant.info[seq(1, length(variant.info), by=4)]))
results$pos <- variant.info[seq(2, length(variant.info), by=4)]
results$a0 <- variant.info[seq(3, length(variant.info), by=4)]
results$a1 <- variant.info[seq(4, length(variant.info), by=4)]
rm(variant.info)

# write out updated results file and put path in flanders input file
input.info[16, "input"] <- paste0("/path/flanders/all_cells.txt")
write.table(results, "/path/flanders/all_cells.txt", sep="\t",
            quote=F, row.names=F)

# Write out input file
write.table(input.info, 
            "/path/flanders/slemap_flanders_input_updated.txt", 
            sep="\t", row.names=F, quote=F)

###############
# to run
# flanders --summarystats_input slemap_flanders_input_updated.txt --run_liftover FALSE --susie_max_iter 1000 --pph4_threshold 0.8 --pph3_threshold 0.8 --chromosomes 1-23 --publish_susie TRUE --outdir flanders_updated_outputs