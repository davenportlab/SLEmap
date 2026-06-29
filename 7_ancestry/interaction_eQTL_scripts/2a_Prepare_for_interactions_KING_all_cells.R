# Prepare files for mapping eQTL interactions using output from QTLight pipeline
# main eQTL mapping
# library(gdata)
library(ggplot2)
library(data.table)
library(glue)
library(dplyr)


symbols_to_names <- read.csv("ensemblID_to_genesymbol.csv")
symbols_to_names$gene_name <- symbols_to_names$gene_symbols

wgs_dir <- "WGS_data"

data_dir <- glue("Interactions_coloc_only/all_cells")

qtlight_dir <- glue("5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/")

load(glue("{data_dir}/preprocess_files.rda"))

coloc_results <- read.csv(glue("{data_dir}/coloc_results.csv"))

lead_coloc_snps <- gsub(";", ".", lead_coloc_snps)

king_ancestry <- read.table(glue("{wgs_dir}/King_out/king_InferredAncestry.txt"), header = T)
covariates <- merge(covariates, king_ancestry[, c("IID", "Ancestry")], by.x = 0, by.y = "IID")

ggplot(data = covariates, aes(x = Genotype_PC1, y = Genotype_PC2, colour = Ancestry)) + geom_point()
ggplot(data = covariates, aes(x = Genotype_PC3, y = Genotype_PC4, colour = Ancestry)) + geom_point()

covariates <- covariates[covariates$Ancestry %in% c("EUR", "AFR", "SAS"),]
rownames(covariates) <- covariates$Row.names
covariates$Row.names <- NULL

# Filter GEX data to indivudals in those ancestries

gex <- gex[colnames(gex) %in% rownames(covariates)]

keep_samples <- colnames(gex)  

geno <- data.frame(fread(glue("{data_dir}/genotypes/sig_snps.raw"),
                         sep="\t", drop = c(1,3:6)))
rownames(geno) <- geno$IID
geno <- geno[keep_samples, ]
colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno))-2)
geno[, 1] <- NULL
geno <- 2 - geno #plink --extract A is giving 0 as the ALT and 2 as the REF, we want 0 as REF and 2 as ALT

geno <- as.matrix(geno)  
# Switch from ensembl ID to gene names

sig_res$Gene <- symbols_to_names$gene_name[match(sig_res$phenotype_id, symbols_to_names$X)]

new_rownames <- symbols_to_names$gene_name[match(rownames(gex), symbols_to_names$X)]

rownames(gex) <- new_rownames

keep_snps <- list()

for (i in 1:length(lead_coloc_snps)){
  
  if (length(unique(keep_samples[which(covariates$Ancestry == "AFR"
                                       & geno[, lead_coloc_snps[[i]]] == 2)])) > 1 &
      length(unique(keep_samples[which(covariates$Ancestry == "EUR" 
                                       & geno[, lead_coloc_snps[[i]]] == 2)])) > 1 &
      length(unique(keep_samples[which(covariates$Ancestry == "SAS" 
                                       & geno[, lead_coloc_snps[[i]]] == 2)])) > 1 &
      length(unique(keep_samples[which(covariates$Ancestry == "AFR" 
                                       & geno[, lead_coloc_snps[[i]]] == 0)])) > 1 &
      length(unique(keep_samples[which(covariates$Ancestry == "EUR" 
                                       & geno[, lead_coloc_snps[[i]]] == 0)])) > 1 &
      length(unique(keep_samples[which(covariates$Ancestry == "SAS" 
                                       & geno[, lead_coloc_snps[[i]]] == 0)])) > 1){
    keep_snps[[length(keep_snps) + 1]] <- lead_coloc_snps[[i]]
  }
}
tested <- length(keep_snps)

print(glue("{length(keep_snps)} out of {length(lead_coloc_snps)} passed filter for testing"))

pass_hom_eqtl <- coloc_results[coloc_results$lead_H4_variant_REF_ALT %in% keep_snps, ]
pass_hom_eqtl$gene_name <- pass_hom_eqtl$gene_symbol


pairs.int <- pass_hom_eqtl[, c("gene_name", "lead_H4_variant_REF_ALT")]
pairs.int$gene_name <- as.character(pairs.int$gene_name)
pairs.int$lead_H4_variant_REF_ALT <- as.character(pairs.int$lead_H4_variant_REF_ALT)
dim(pairs.int)

qvals <- read.delim(glue("{qtlight_dir}/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"))
qvals <- qvals[qvals$phenotype_id %in% pass_hom_eqtl$gene_id, ]
qvals <- merge(qvals, symbols_to_names[, c("gene_name", "X")], by.x = "phenotype_id", by.y = "X")

int_file <- glue("{data_dir}/eqtl_interact_KING_files_pass_hom_filt.rda")
save(list=c("gex", "geno", "covariates", "pairs.int", "pass_hom_eqtl", "sig_res", "qvals"), file = int_file)

