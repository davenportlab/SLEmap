library(data.table)
library(dplyr)

DIR_gwas_resources <- "/path/coloc/resources/"
f <- fread(paste0(DIR_gwas_resources,
                  "/GCST90270940/GWAS_results.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::select(., chr, `pos(hg19)`) %>%
  dplyr::rename(., start=`pos(hg19)`) %>%
  dplyr::mutate(end=start) %>%
  dplyr::mutate(chr=paste0("chr", chr))
  
write.table(f, 
            paste0(DIR_gwas_resources,
                   "/GCST90270940/input_hg19.bed"),
            sep="\t", quote=F, col.names = F, row.names = F)


