## this is to see if our coloc results are novel/known by comparing with otar db
## run 03_OTAR/1_otar_coloc.ipynb first
library(readr)
library(stringr)
DIR_MAIN="/path/coloc/outputs"

##get otar colocs 
otar <- read_csv("/path/03_OTAR/coloc_alreadyknown_otar_uniquevariants.csv")
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

write.csv(coloc,"/path/coloc/coloc_output_otar.csv")


