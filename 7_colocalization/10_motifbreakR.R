###running motifbreakR for snps where coloc is only in slemap (but eQTL in both slemap and onek1k) to see if the eSNV disrupts different tf binding sites. 

library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
library(data.table)
library(tidyverse)

resultdir <- "/path/motifbreakr/"
dir.create(file.path(resultdir))
dir.create(file.path(paste0(resultdir,"bedfiles")))
dir.create(file.path(paste0(resultdir,"strong_results")))

####identify eSNVs to test 
coloc <- read.csv(paste0("/path/coloc/","1_csvfiles/coloc_sigresults.csv"),row.names=1)

####running for specific cases1 - TRAF1 in CD4 T cells####
coloc_use <- coloc[coloc$cell_type == "CM_CD4_T_cells" & coloc$gene_symbol == "TRAF1",]
coloc_use
lead_h4variant <- paste0(coloc_use$lead_H4_variant,"_",coloc_use$lead_H4_variant_GTF_REF,"_",coloc_use$lead_H4_variant_GTF_ALT)
gene_id <- coloc_use$gene_id


##make bed file
vars <- c(lead_h4variant)
parts <- do.call(rbind, strsplit(vars, "_"))
bed <- data.frame(
  chrom = parts[,1],
  start = as.integer(parts[,2]) - 1,
  end = as.integer(parts[,2]),
  name = paste(parts[,1], parts[,2], parts[,3], parts[,4], sep=":"),
  score = 0,
  strand = ".",
  stringsAsFactors = FALSE
)
write.table(bed, paste0(resultdir,"bedfiles/TRAF1.bed"), sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE)

snps.mb.frombed <- snps.from.file(file = paste0(resultdir,"bedfiles/TRAF1.bed"),
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed")

##run motif breakr
data("encodemotif")
results <- motifbreakR(snpList = snps.mb.frombed, filterp = TRUE,
                       pwmList = encodemotif,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())

strong_res <- results[mcols(results)$effect == "strong"]
table(strong_res$SNP_id)
results_df <- as.data.frame(strong_res,row.names=NULL)
results_df$cell_type[results_df$SNP_id == gsub("_",":",lead_h4variant)] <- "CM_CD4_T_cells"
#results_df$cell_type[results_df$SNP_id == gsub("_",":",Bcell_variant)] <- "Memory_B_cells"
results_df[] <- lapply(results_df, function(x) {
  if (is.list(x)) sapply(x, toString) else x
})
write.csv(results_df,paste0(resultdir,"strong_results/TRAF1_CM_CD4_T_cells.csv"),row.names=F)
View(results_df[,c("SNP_id","geneSymbol","dataSource","effect","alleleEffectSize","cell_type")])

pdf("/path/coloc/outputs/0_plots/TRAF1_motifbreakr_plot.pdf", width = 6, height = 7.7)
plotMB(results = results, rsid = gsub("_",":",lead_h4variant), effect = "strong")
dev.off()