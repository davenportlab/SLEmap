## this is to see if our coloc results are novel/known by comparing with otar db
## known coloc = lead coloc snp from otar is in gwas loci (defined by locus breaker) 
## run /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/19_unique_coloc/1_otar_coloc.ipynb in local first
library(readr)
library(stringr)
DIR_MAIN="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs"

##get otar colocs 
otar <- read_csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/coloc_alreadyknown_otar_uniquevariants.csv")
otar$chr <- paste0("chr", sub("_.*", "",otar$variantId))
otar$pos <- sub("^[^_]+_([^_]+)_.*", "\\1", otar$variantId)

#get my colocs
coloc <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_sigresults_GCST90270940_checksigeQTL.csv"))

#compare
window <- 500000

for (i in seq_len(nrow(coloc))) {
  chr_i <- coloc$lead_snp_chr[i]
  pos_i <- coloc$lead_snp_pos[i]
  
  hits <- otar$chr == chr_i &
    otar$pos >= (pos_i - window) &
    otar$pos <= (pos_i + window)
  
  coloc$otar_leadsnp_overlap[i] <- any(hits)
  
  ids <- unique(otar$projectId[
    otar$chr == chr_i &
      otar$pos >= (pos_i - window) &
      otar$pos <= (pos_i + window)
  ])
  
  coloc$otar_projectIds[i] <- if (length(ids)) paste(ids, collapse = ",") else NA
  
  genes <- unique(otar$geneId[
    otar$chr == chr_i &
      otar$pos >= (pos_i - window) &
      otar$pos <= (pos_i + window)
  ])
  
  coloc$otar_geneIds[i] <- if (length( genes)) paste( genes, collapse = ",") else NA
}

coloc$gene_overlap <- mapply(function(gene, gene_list) {
  gene %in% str_split(gene_list, ",")[[1]]
}, coloc$gene_id, coloc$otar_geneIds)

coloc$otar_overlap <- coloc$leadsnp_overlap & coloc$gene_overlap

write.csv(coloc,"/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/coloc_output_otar.csv")

###
coloc <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/19_unique_coloc/coloc_output_otar.csv")
otar_overlap_gene <- as.data.frame(table(coloc$gene_id,coloc$otar_overlap))
colnames(otar_overlap_gene) <- c("gene","overlap","numcoloc")
nrow(otar_overlap_gene[otar_overlap_gene$overlap == F,]) ##allgenes
nrow(otar_overlap_gene[(otar_overlap_gene$overlap == F)&(otar_overlap_gene$numcoloc >0),]) ##novel genes

coloc_filtered <- coloc[!coloc$gene_symbol %in% c("KANSL1","KANSL1-AS1","ARL17B"),]
otar_overlap_gene <- as.data.frame(table(coloc_filtered$gene_id,coloc_filtered$otar_overlap))
colnames(otar_overlap_gene) <- c("gene","overlap","numcoloc")
nrow(otar_overlap_gene[otar_overlap_gene$overlap == F,]) ##allgenes
nrow(otar_overlap_gene[(otar_overlap_gene$overlap == F)&(otar_overlap_gene$numcoloc >0),]) ##novel genes


