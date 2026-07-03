library(tidyverse)
library(data.table)

DIR_tensorQTL="/path/eQTLresults"
DIR_MAIN="/path/coloc/"

####################################################
## 1. load data: bim_and_maf
####################################################
source(paste0(DIR_MAIN,
              "/scripts/00_functions/load_data.R"))
print(dim(bim_and_maf))

###########################################################################
## 2. combine all cis_nominal results from OPTIM_pcs
###########################################################################

all_cis_eqtl_conditional <- list.files(path=DIR_tensorQTL,
                                       pattern ="Cis_eqtls_independent.tsv",
                                       #pattern ="eqtls",
                                       recursive = T,
                                       full.names = T)
all_cis_eqtl_conditional <- all_cis_eqtl_conditional[grep("OPTIM_pcs", all_cis_eqtl_conditional)]
all_cis_eqtl_conditional <- all_cis_eqtl_conditional[grep("results", all_cis_eqtl_conditional)]

all_cis_eqtl_conditional2 <- data.frame(cell_type=gsub(".*dMean\\_\\_(.*)\\_all\\/OPT.*", "\\1", all_cis_eqtl_conditional),
                                        file=all_cis_eqtl_conditional)

mainDir=paste0(DIR_MAIN, 
               "/02_2.load_eqtls_data/")
for(i in 1:nrow(all_cis_eqtl_conditional2) ) {
  
  cell_type <- all_cis_eqtl_conditional2$cell_type[i]
  print(cell_type)
  dir.create(file.path(mainDir, cell_type))
  #file.copy(cis_out, paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
  
  f <- fread(all_cis_eqtl_conditional2$file[i]) %>% as.data.frame()
  gene_ids <- unique(f$phenotype_id)
  length(gene_ids)
  
  dir_cell_type = paste0(DIR_tensorQTL, "/", cell_type, "/results")
  
  get_e_genes <- do.call(rbind, lapply(list.files(dir_cell_type, 
                                                  pattern="cis_nominal1.cis_qtl_pairs.*", 
                                                  recursive = T, 
                                                  full.names = T), 
                                       function(file.name) {
                                         fread(paste0(file.name)) %>%
                                           as.data.frame() %>%
                                           dplyr::filter(phenotype_id %in% gene_ids) %>%
                                           dplyr::mutate(chr_pos = gsub("(chr.*\\_\\d+)\\_.*\\_.*", "\\1", variant_id))
                                       })) %>%
    mutate(., tmp_n=1:nrow(.)) %>%
    #merge(., bim_and_maf, by.x="variant_id", by.y="snp") %>%
    merge(., bim_and_maf, by="chr_pos") %>%
    #merge(., bim_and_maf, by.x="variant_id", by.y="chr_pos_variants") %>%
    arrange(., tmp_n)
  
  dim(get_e_genes)
  
  OUT_FILE=gzfile( paste0( mainDir, "/", 
                           cell_type, 
                           "/cis_nominal_", cell_type, ".txt.gz")  )
  write.table(get_e_genes, OUT_FILE,
              quote=F, row.names=F, sep="\t")
}


###########################################################################
## 3. get conditional eQTL data
###########################################################################
# th= 1e-04
th= 0.05

for(i in 2:nrow(all_cis_eqtl_conditional2)) {
  print(i)
  
  cell_type <- all_cis_eqtl_conditional2$cell_type[i]
  print(cell_type)
  
  dir.create(file.path(mainDir, cell_type))
    
  # 1. copy Cis_eqtls_qval.tsv
  file.copy(gsub("Cis_eqtls_independent.tsv", "Cis_eqtls_qval.tsv", 
                 all_cis_eqtl_conditional2$file[i]), 
            paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
  
  # 2. copy Cis_eqtls.tsv
  file.copy(gsub("Cis_eqtls_independent.tsv", "Cis_eqtls.tsv", 
                 all_cis_eqtl_conditional2$file[i]), 
            paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
  
  f <- fread(gsub("Cis_eqtls_independent.tsv", "Cis_eqtls_qval.tsv", 
                  all_cis_eqtl_conditional2$file[i])) %>% 
    as.data.frame() %>%
    dplyr::filter(., qval < th) # pval_perm or 
  gene_ids <- unique(f$phenotype_id)
  length(gene_ids)
  dim(f)
  
  tmp <- f %>% dplyr::select(., phenotype_id, variant_id, pval_nominal, pval_perm, qval ) %>%
    merge(., bim_and_maf, by.x="variant_id", by.y="chr_pos_variants") %>%
    mutate(cell_type=cell_type)
  dim(tmp)
  # need to check this
  nrow(f) == nrow(tmp)
  
  # a <- tmp %>% group_by(variant_id)  %>% filter(n() > 1)
  
  
  OUT_FILE=gzfile( paste0(mainDir, "/",
                          cell_type, "/Cis_eqtls_qval_", th, "_POS_added.txt.gz")) 
  write.table(tmp, OUT_FILE,
              quote=F, row.names=F, sep="\t")
  
  # 3. copy Cis_eqtls_independent.tsv
  file.copy(all_cis_eqtl_conditional2$file[i],
            paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )

  f <- fread(all_cis_eqtl_conditional2$file[i]) %>% as.data.frame()
  gene_ids <- unique(f$phenotype_id)
  length(gene_ids)
  dim(f)
  tmp <- f %>% dplyr::select(., phenotype_id, variant_id, rank ) %>%
    merge(., bim_and_maf, by.x="variant_id", by.y="chr_pos_variants") %>%
    mutate(cell_type=cell_type)
  dim(tmp)

  # need to check this
  nrow(f) == nrow(tmp)

  OUT_FILE=gzfile( paste0(mainDir, "/",
                          cell_type, "/Cis_eqtls_independent_POS_added.txt.gz"))
  write.table(tmp, OUT_FILE,
              quote=F, row.names=F, sep="\t")

}

for(i in 2:nrow(all_cis_eqtl_conditional2)){
  print(i)
  cell_type <- all_cis_eqtl_conditional2$cell_type[i]
  print(cell_type)
  norm <- paste0(DIR_tensorQTL, "/", 
                 cell_type, "/results/norm_data/dMean__", 
                 cell_type, "_all/normalised_phenotype.tsv")
  norm2 <- paste0( gsub("normalised_phenotype.tsv", "", norm),
                   "dMean__", cell_type, "___phenotype_file.tsv" )

  
  dir.create(file.path(mainDir, cell_type))
  file.copy(norm, 
            paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
  file.copy(norm2, 
            paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
}




