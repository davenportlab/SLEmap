library(tidyverse)
library(data.table)

DIR_MAIN="/path/coloc"

# gwas data
GWAS_INPUT=fread(paste0(DIR_MAIN,
                        "/resources/GCST90270940/GWAS_results.txt.gz")) %>%
  as.data.frame() %>% dplyr::mutate(hg19_id=paste0("chr", chr, "_", `pos(hg19)`))

input_hg19=fread(paste0(DIR_MAIN,
                        "/resources/GCST90270940/input_hg19.bed")) %>% 
  as.data.frame() %>% dplyr::select(., V1, V2)
colnames(input_hg19) <- c("chr_hg19", "pos_hg19")
dim(input_hg19)

output_hg38=fread(paste0(DIR_MAIN,
                         "/resources/GCST90270940/output_hg38.bed")) %>% 
  as.data.frame() %>% dplyr::select(., V1, V2)
colnames(output_hg38) <- c("chr_hg38", "pos_hg38")
dim(output_hg38)

unMapped=fread(paste0(DIR_MAIN,
                      "/resources/GCST90270940/unMapped.bed")) %>% as.data.frame()

# filter unapped positions
input_hg19_ft <- input_hg19 %>% 
  dplyr::filter(., !(chr_hg19=="chr6" & pos_hg19==51097410) )

# combine hg19 and hg38
hg19_hg38 <- data.frame(input_hg19_ft[, 1:2], output_hg38[,1:2]) %>%
  dplyr::mutate(hg19_id=paste0(chr_hg19, "_", pos_hg19)) %>%
  dplyr::select(chr_hg38, pos_hg38, hg19_id)
dim(hg19_hg38)

GWAS_hg38 <- merge(GWAS_INPUT, hg19_hg38, by="hg19_id", all=TRUE) %>%
  dplyr::arrange(., chr, pos_hg38)

out_file=paste0(DIR_MAIN,
                "/resources/GCST90270940/GCST90270940_hg38.txt.gz")
OUT_FILE=gzfile(out_file )
write.table(GWAS_hg38, OUT_FILE, sep="\t", quote=F, row.names = F)

