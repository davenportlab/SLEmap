library(tidyverse)
library(data.table)

DIR_tensorQTL="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs"
DIR_tensorQTL_all="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells"

mainDir=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Colocalisation/load_eqtls_data")

#DIR_MAIN="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/wl2/colocalization/"

####################################################
## 1. load data: bim_and_maf
####################################################

##To get maf and N_CHR info for coloc run the following
#plink2 --vcf /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/eb30/results/oneK1K_imputation/Imputation_TOPMed_run3/plink_conversion_QC2/onek1k_imputed_allchr_afterQC2_new_filter.vcf.gz --freq --update-sex /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/eb30/results/oneK1K_imputation/Imputation_TOPMed_run3/plink_conversion_QC2/sex.txt --out /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Genotypes/freq/onek1k_imputed_allchr_afterQC2_new_filter
freq <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Genotypes/freq/onek1k_imputed_allchr_afterQC2_new_filter.afreq")
colnames(freq) <- c("CHROM","variant_id","REF","ALT","ALT_FREQS","OBS_CT")
freq <- freq %>%
  mutate(
    ALT_FREQS = as.numeric(ALT_FREQS),
    REF_FREQ  = 1 - ALT_FREQS,
    minor_allele = if_else(ALT_FREQS <= REF_FREQ, ALT, REF),
    minor_allele_frq = pmin(ALT_FREQS, REF_FREQ)
  )

###########################################################################
## 2. combine all cis_nominal results from OPTIM_pcs
###########################################################################

all_cis_eqtl_conditional <- list.files(path=DIR_tensorQTL,
                                       pattern ="Cis_eqtls_independent.tsv",
                                       #pattern ="eqtls",
                                       recursive = T,
                                       full.names = T)
all_cis_eqtl_conditional <- all_cis_eqtl_conditional[!grepl("/work/", all_cis_eqtl_conditional)]
all_cis_eqtl_conditional <- append(all_cis_eqtl_conditional,"/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv")

all_cis_eqtl_conditional2 <- data.frame(cell_type=gsub(".*dMean\\_\\_(.*)\\_all\\/OPT.*", "\\1",all_cis_eqtl_conditional),file=all_cis_eqtl_conditional)

for(i in 1:nrow(all_cis_eqtl_conditional2) ) {
  
  cell_type <- all_cis_eqtl_conditional2$cell_type[i]
  print(cell_type)
  dir.create(file.path(mainDir, cell_type))
  #file.copy(cis_out, paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
  
  f <- fread(all_cis_eqtl_conditional2$file[i]) %>% as.data.frame()
  gene_ids <- unique(f$phenotype_id)
  length(gene_ids)
  
  if(cell_type == "All"){
    dir_cell_type = paste0(DIR_tensorQTL_all, "/results")
  }else{
    dir_cell_type = paste0(DIR_tensorQTL, "/", cell_type, "/results")
  }
  
  get_e_genes <- do.call(rbind, lapply(list.files(dir_cell_type, 
                                                  pattern="cis_nominal1.cis_qtl_pairs.*", 
                                                  recursive = T, 
                                                  full.names = T), 
                                       function(file.name) {
                                         fread(paste0(file.name)) %>%
                                           as.data.frame() %>%
                                           dplyr::filter(phenotype_id %in% gene_ids) %>%
                                           dplyr::mutate(chr_pos = gsub("(chr.*\\_\\d+)\\_.*\\_.*", "\\1", variant_id))%>%
                                           separate_wider_delim(variant_id,delim = "_",names = c("CHROM", "POS", "GTF_ALT", "GTF_REF"),cols_remove = FALSE) #REF ALT is flipped on purpose due to legacy code. This is flipped back and checked with the gtf file before interpretation. 
                                       })) %>%
    mutate(., tmp_n=1:nrow(.)) %>%
    merge(., freq[,c("variant_id","OBS_CT","minor_allele","minor_allele_frq")], by="variant_id") %>%
    #merge(., bim_and_maf, by.x="variant_id", by.y="snp") %>%
    #merge(., bim_and_maf, by="chr_pos") %>%
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

for(i in 1:nrow(all_cis_eqtl_conditional2)) {
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
    #merge(., bim_and_maf, by.x="variant_id", by.y="chr_pos_variants") %>%
    separate_wider_delim(variant_id,delim = "_",names = c("CHROM", "POS", "GTF_ALT", "GTF_REF"),cols_remove = FALSE) %>% #REF ALT is flipped on purpose due to legacy code. This is flipped back before interpretation. 
    merge(., freq[,c("variant_id","OBS_CT","minor_allele","minor_allele_frq")], by="variant_id") %>%
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
    separate_wider_delim(variant_id,delim = "_",names = c("CHROM", "POS", "GTF_ALT", "GTF_REF"),cols_remove = FALSE) %>% #REF ALT is flipped on purpose due to legacy code. This is flipped back before interpretation. 
    #merge(., bim_and_maf, by.x="variant_id", by.y="chr_pos_variants") %>%
    mutate(cell_type=cell_type)
  dim(tmp)
  
  # need to check this
  nrow(f) == nrow(tmp)
  
  OUT_FILE=gzfile( paste0(mainDir, "/",
                          cell_type, "/Cis_eqtls_independent_POS_added.txt.gz"))
  write.table(tmp, OUT_FILE,
              quote=F, row.names=F, sep="\t")
  
  if(cell_type == "All"){
    norm <- paste0(DIR_tensorQTL_all, "/results/norm_data/dMean__", 
                   cell_type, "_all/normalised_phenotype.tsv")
  }else{
    norm <- paste0(DIR_tensorQTL, "/", 
                   cell_type, "/results/norm_data/dMean__", 
                   cell_type, "_all/normalised_phenotype.tsv")
  }
  
  norm2 <- paste0( gsub("normalised_phenotype.tsv", "", norm),
                   "dMean__", cell_type, "___phenotype_file.tsv" )
  
  file.copy(norm, 
            paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
  file.copy(norm2, 
            paste0(mainDir, "/", cell_type, "/"), overwrite = TRUE )
}


