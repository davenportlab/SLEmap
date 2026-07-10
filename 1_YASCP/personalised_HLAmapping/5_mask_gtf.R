# primary gtf

setwd("/path/single_cells/scripts/02_1.personalized_mapping")

gtf.primary <- fread(ref_genes) %>% 
  as.data.frame()
gtf.primary <- gtf.primary[order(gtf.primary$`##date: 2019-09-05`, gtf.primary$V4 ),]
write.table(gtf.primary,
            paste0(DIR_single_cell,
                   "/02_1.personalized_mapping/pers_refs/masked/tmp.sorted.GRCh38.primary.gtf"),
            sep="\t", row.names = F, quote=F)

gtf.primary <- rtracklayer::import(paste0(DIR_single_cell,
                                          "/02_1.personalized_mapping/pers_refs/masked/tmp.sorted.GRCh38.primary.gtf"))
hla_genes_to_mask <- read.table(paste0(DIR_single_cell,
                                       "/02_1.personalized_mapping/pers_refs/masked/Ensembl98.primary.to_mask.bed"))
HLA_genes.masked <- gtf.primary[!(gtf.primary$gene_name %in% hla_genes_to_mask$V4), ]


rtracklayer::export(HLA_genes.masked, 
                    paste0(DIR_single_cell,
                           "/02_1.personalized_mapping/pers_refs/masked/GRCh38.primary.Ensembl98.HLA_masked.sorted.gtf") )
