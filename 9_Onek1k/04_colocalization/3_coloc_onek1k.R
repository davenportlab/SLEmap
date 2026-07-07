##onek1k coloc script with locus breaker windows

library(tidyverse)
library(data.table)
library(coloc)

GWAS_input= "/path/colocalization/inputs/gwas/"
DIR_output="/path/onek1k_locus_breaker_coloc"
eQTL_nominal_p_source="/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs/"
eQTL_nominal_p_source_allcells="/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/"

GWAS_threshold=1e-5
minSNPs = 1

######colocalisation functions
get_window_data <- function(data, GWAShit_chr, GWAShit_pos, gene_id, locusStart, locusEnd) {
  subset_data <- data %>% as.data.frame(.) %>% 
    dplyr::filter(., phenotype_id==gene_id )
  
  
  windowed_data <- dplyr::filter(subset_data, 
                                 CHROM==GWAShit_chr & 
                                   POS >= locusStart & 
                                   POS <= locusEnd ) %>%
    dplyr::arrange(., POS) 
  return(windowed_data)
}

get_window_data_v2 <- function(data, GWAShit_chr, GWAShit_pos, locusStart, locusEnd) {
  # get chromosome and position columns
  data$chr <- as.character(data$chr)
  windowed_data <- data[data$chr==GWAShit_chr & 
                          data$pos >= locusStart &  
                          data$pos <= locusEnd,]
  return(windowed_data)
}

prepare_data_eQTL <- function(data, type) {
  data <- data %>% 
    dplyr::mutate(., N=OBS_CT/2) %>%
    ungroup() %>%
    dplyr::select(variant_id, slope, slope_se, pval_nominal, minor_allele_frq, N)
  
  colnames(data) <- c("SNP", "beta", "se", "pval", "MAF", "N")
  data$se <- as.numeric(data$se)
  data$beta <- as.numeric(data$beta)
  data$pval <- as.numeric(data$pval)
  
  ## add MAF
  coloc_format <- list(beta = data$beta,varbeta = data$se^2,snp = data$SNP,type = type,N = data$N,MAF = data$MAF)
  return(coloc_format)
}

prepare_data_GWAS <- function(data, type, N, s) {
  data <- data %>% 
    dplyr::mutate(., N=N) %>%
    dplyr::select(variant_id, beta, standard_error, p_value, effect_allele_frequency)
  
  colnames(data) <- c("SNP", "beta", "se", "pval", "effect_allele_frequency")
  data$se <- as.numeric(data$se)
  data$beta <- as.numeric(data$beta)
  data$pval <- as.numeric(data$pval)
  
  ## add MAF
  coloc_format <- list(pvalues = as.numeric(data$pval),beta = as.numeric(data$beta),varbeta = data$se^2,snp = data$SNP,type = "cc",N = N,s = s)
  return(coloc_format)
}

coloc_between_eqtl_and_gwas <- function(eQTL_INPUT, GWAS_INPUT, GWAS_ID, N, s, ld_eqtl,
                                        cell_type, locusStart, locusEnd, gene_name, gene_id, 
                                        lead_snp, lead_snp_chr, lead_snp_pos,
                                        scaled_expression, mean_expression,
                                        do_plot="N",GWAShit_chr,GWAShit_pos,nominal_pval_thres){
  
  # scaled_expression, mean_expression : these are for plotting
  input_title1=paste0("H1:", cell_type)
  input_title2=paste0("H2:", GWAS_ID) 
  print(input_title1)
  print(input_title2)
  
  
  ## 1) filter variants
  
  # Subset data around a specific lead SNP
  eqtl_window_1 <- get_window_data(data=eQTL_INPUT, 
                                   GWAShit_chr=GWAShit_chr,
                                   GWAShit_pos=GWAShit_pos, 
                                   gene_id=gene_id, 
                                   locusStart=locusStart,
                                   locusEnd=locusEnd) #contains duplicates
  dim(eqtl_window_1)
  in_eqtl_window <- eqtl_window_1 %>% 
    dplyr::arrange(., pval_nominal)
  
  gwas_window_1 <- get_window_data_v2(GWAS_INPUT, 
                                      GWAShit_chr=GWAShit_chr, 
                                      GWAShit_pos=GWAShit_pos, 
                                      locusStart=locusStart,
                                      locusEnd=locusEnd) %>%
    dplyr::rename(., POS=pos)
  dim(gwas_window_1)
  in_gwas_window <- gwas_window_1 %>% 
    dplyr::arrange(., p_value)
  
  gwas_window_1 <- gwas_window_1[!duplicated(gwas_window_1[c("POS")]), ] 
  
  ## 2)extract dup in eqtl and check alleles with gwas
  eqtl_window_1_notdup <- eqtl_window_1 %>%group_by(across(POS)) %>%dplyr::filter(n() == 1) %>% ungroup()
  eqtl_window_1_dup <- eqtl_window_1 %>%group_by(across(POS)) %>%dplyr::filter(n() != 1)%>% ungroup()
  eqtl_window_1_dup <- eqtl_window_1_dup[!duplicated(eqtl_window_1_dup[c("variant_id_v2")]), ]
  eqtl_window_1_dup <- eqtl_window_1_dup[eqtl_window_1_dup$variant_id_v2 %in% c(paste0(gwas_window_1$variant_id,"_",gwas_window_1$other_allele,"_",gwas_window_1$effect_allele), paste0(gwas_window_1$variant_id,"_",gwas_window_1$effect_allele,"_",gwas_window_1$other_allele)),] #for eQTLs with multiple minor alleles, use result that matches the gwas allele
  eqtl_window_1 <- rbind(eqtl_window_1_notdup,eqtl_window_1_dup)
  
  dim(eqtl_window_1)
  dim(gwas_window_1)
  
  
  ## 3) get common variants 
  common_snps <- intersect(eqtl_window_1$variant_id, gwas_window_1$variant_id)
  print ( paste0("common SNPs between datasets: ", length(common_snps) ) )
  
  if (length(common_snps) < minSNPs) {
    print(paste0("Too small number of common SNPs between datasets:", length(common_snps)))
    return("Too small number of common SNPs between datasets")
  }
  
  in_eqtl_window$common <- ifelse(in_eqtl_window$variant_id %in% common_snps, "Y", "N")
  in_gwas_window$common <- ifelse(in_gwas_window$variant_id %in% common_snps, "Y", "N")
  
  
  ## 4) number of variants for coloc
  cat("Running coloc analysis with", length(common_snps), "common SNPs\n")
  
  eqtl_window_1 <- eqtl_window_1[eqtl_window_1$variant_id %in% common_snps,] %>% 
    dplyr::arrange(POS)
  gwas_window_1 <- gwas_window_1[gwas_window_1$variant_id %in% common_snps,] %>% 
    dplyr::arrange(POS)
  dim(eqtl_window_1)
  dim(gwas_window_1)
  
  all(eqtl_window_1$variant_id==gwas_window_1$variant_id)
  all(eqtl_window_1$POS==gwas_window_1$POS)
  all(eqtl_window_1$variant_id_v2==paste0(gwas_window_1$variant_id,"_",gwas_window_1$effect_allele,"_",gwas_window_1$other_allele))
  
  if (all(eqtl_window_1$POS==gwas_window_1$POS)==FALSE) {
    print("SNP positions are different between datasets")
    return("SNP positions are different between datasets")
  }
  
  if( all(eqtl_window_1$variant_id_v2==paste0(gwas_window_1$variant_id,"_",gwas_window_1$effect_allele,"_",gwas_window_1$other_allele)) == FALSE){
    print("effect allele different between datasets")
    return("effect allele different between datasets")
  }
  
  if( min(gwas_window_1$p_value) > GWAS_threshold ){
    print("minimum GWAS p value among common SNPs is bigger than 1e-5") # this can still happen when there's no variants < 1e-5 in the common variants
    return("minimum GWAS p value among common SNPs is bigger than 1e-5")
  }
  
  if( min(eqtl_window_1$pval_nominal) > nominal_pval_thres ){
    print("minimum eQTL p value among common SNPs is bigger than nominal pval threshold") 
    return("minimum eQTL p value among common SNPs is bigger than nominal pval threshold")
  }
  
  ## 5) prepare input data for coloc
  coloc_input_1 <- prepare_data_eQTL(data=eqtl_window_1, type="quant")
  coloc_input_2 <- prepare_data_GWAS(data=gwas_window_1, type="cc", 
                                     N=N, s=s)
  
  ## 6) Perform colocalization
  coloc_results <- coloc.abf(
    dataset1 = coloc_input_1,
    dataset2 = coloc_input_2
  )
  
  ## 7) summerize results
  print(coloc_results$summary)
  
  output <- data.frame(cell_type=cell_type,
                       gene_id=gene_id,
                       lead_snp=lead_snp,
                       lead_snp_chr=lead_snp_chr,
                       lead_snp_pos=lead_snp_pos,
                       gwas_hit=paste0(GWAShit_chr,"_",GWAShit_pos),
                       t(coloc_results$summary)) 
  dashboard <- "no_fig"
  
  return(all_output = list(coloc_all=coloc_results,
                           coloc_summary=output,
                           eqtl_window=eqtl_window_1,
                           gwas_window=gwas_window_1,
                           out_fig=dashboard,
                           input_eqtl=in_eqtl_window,
                           input_gwas=in_gwas_window,
                           variants_for_coloc=common_snps))
  
}

## specify GWAS
GWAS_ID="GCST90270940"
N = 5877 + 188588 + 14355 + 505956 + 1393 + 2327
s = (5877 + 14355 + 1393) / (188588 + 505956 + 2327)

## make output directory
DIR_OUT = paste0(DIR_output, 
                 "/", GWAS_ID)
print(DIR_OUT)
if (!dir.exists(DIR_OUT)) {
  dir.create(DIR_OUT)
}
if (!dir.exists(paste0(DIR_OUT, "/all_checksigeQTL_checkallele"))) {
  dir.create(paste0(DIR_OUT, "/all_checksigeQTL_checkallele"))
}
if (!dir.exists(paste0(DIR_OUT, "/sig_checksigeQTL_checkallele"))) {
  dir.create(paste0(DIR_OUT, "/sig_checksigeQTL_checkallele"))
}

#### load GWAS data- and remove HLA region (ref 26-24Mb: https://www.nature.com/articles/ng.3434#MOESM28 - increased it to start at 25Mb to cover whole region)
GWAS_INPUT <- fread (paste0(GWAS_input,GWAS_ID,"/",GWAS_ID,"_for_eqtl.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::rename(., variant_id_v2=variant_id) %>%
  dplyr::mutate(variant_id = paste0(chr, "_", pos))# 7071163

GWAS_INPUT <- GWAS_INPUT[!((GWAS_INPUT$chr == "chr6") & (GWAS_INPUT$pos > 25000000)&(GWAS_INPUT$pos < 34000000)),] 
dim(GWAS_INPUT)
# 7018195

GWAS_INPUT <- GWAS_INPUT[is.na(GWAS_INPUT$variant_id) == F,]

locusbreaker <- read.csv("/path/locus_breaker_coloc/gwashit_locusbreaker_forinput.csv")

#### run coloc ####
cell_types <- c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","Nonclassical_Monocytes","MAIT_and_GammaDelta_T_cells","All")

for(cell_type in cell_types){
  DIR_eqtl_data=paste0("/path/Onek1k/Colocalisation/load_eqtls_data/")
  
  ## 1) conditional output
  cis_eqtl_f <- paste0(DIR_eqtl_data,
                       cell_type,
                       "/Cis_eqtls_independent_POS_added.txt.gz") #significant results
  cis_eqtl <- fread(cis_eqtl_f) %>%
    as.data.frame() %>%
    dplyr::rename(., variant_id_v2=variant_id) %>%
    dplyr::mutate(variant_id = paste0(CHROM, "_", POS))
  dim(cis_eqtl)
  
  ## add nominal pval threshold
  if(cell_type == "All"){
    nominal_pval_thres <- read.table(paste0(eQTL_nominal_p_source_allcells,"/results/TensorQTL_eQTLS/dMean__",cell_type,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"),fill=T,header=T)
  }else{
    nominal_pval_thres <- read.table(paste0(eQTL_nominal_p_source,cell_type,"/results/TensorQTL_eQTLS/dMean__",cell_type,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"),fill=T,header=T)
  }
  cis_eqtl$nominal_pval_thres <- nominal_pval_thres$pval_nominal_threshold[match(cis_eqtl$phenotype_id,nominal_pval_thres$phenotype_id)]
  
  ## 1-1) conditional output : only 1 eSNP per eGene
  cis_eqtl_1 <- cis_eqtl %>%
    dplyr::group_by(phenotype_id) %>%
    dplyr::filter(n() == 1) %>%
    ungroup() %>% 
    dplyr::mutate(eQTL_group="A") %>%
    as.data.frame()
  
  ## 1-2) conditional output : multiple eSNPs per eGene
  cis_eqtl_2 <- cis_eqtl %>%
    dplyr::group_by(phenotype_id) %>%
    dplyr::filter(n() > 1) %>%
    ungroup() %>% 
    dplyr::mutate(eQTL_group="B") %>%
    as.data.frame()
  length(unique(cis_eqtl_2$phenotype_id))
  
  cis_eqtl <- rbind(cis_eqtl_1, cis_eqtl_2)
  
  ## 3-1) load eqtl nominal p-values
  ## This is for eGenes that are associated with only one eSNP.
  eQTL_INPUT_1 <- fread( paste0(DIR_eqtl_data, 
                                cell_type,  "/cis_nominal_", cell_type, ".txt.gz") ) %>%
    as.data.frame() %>%
    dplyr::rename(., variant_id_v2=variant_id) %>%
    dplyr::rename(variant_id = chr_pos)
  dim(eQTL_INPUT_1)
  
  
  ## 3-2) load eqtl nominal p-values
  ## This is for eGenes that are associated with multiple eSNP.
  
  eQTL_INPUT_2 <- fread( paste0(DIR_eqtl_data, 
                                cell_type,
                                "/cis_nominal_", cell_type, ".multiple_eSNPs.txt.gz") ) %>%
    as.data.frame()
  
  dim(eQTL_INPUT_2)
  head(eQTL_INPUT_2)
  
  
  ## 4) scaled_expression
  scaled_expression = read.table(paste0(DIR_eqtl_data,
                                        cell_type, "/normalised_phenotype.tsv"),
                                 sep="\t")
  dim(scaled_expression)
  
  ####################################################
  ## Starting coloc
  ####################################################
  
  ## define a name of output file
  
  #OUT_FILE <- paste0(DIR_OUT, "/all/", cell_type, ".txt")
  OUT_FILE <- paste0(DIR_OUT, "/all_checksigeQTL_checkallele/", cell_type, ".txt")
  print(OUT_FILE)
  out_to_write <- data.frame(cell_type = character(),
                             gene_id =  character(),
                             lead_snp =  character(),
                             lead_snp_chr =  character(),
                             lead_snp_pos =  integer(),
                             gwas_hit = character(),
                             nsnps =  numeric(),
                             PP.H0.abf =  numeric(),
                             PP.H1.abf =  numeric(),
                             PP.H2.abf =  numeric(),
                             PP.H3.abf =  numeric(),
                             PP.H4.abf =  numeric(),
                             lead_H4_variant =  character(),
                             lead_H4_SNP.PP.H4 =  numeric(),
                             n_variants_for_coloc =  integer(),
                             eQTL_group=character(),
                             eQTL_rank=numeric(),
                             gwas_inwindow=numeric(),
                             status= character())
  write.table(out_to_write, 
              OUT_FILE,
              sep="\t", quote=F, row.names = F, col.names = T,append = F)
  
  for(i in 1:nrow(cis_eqtl)) {
    
    print(i)
    
    gene_name = cis_eqtl$phenotype_id[i]
    gene_id = cis_eqtl$phenotype_id[i]
    lead_snp = cis_eqtl$variant_id[i]
    lead_snp_chr = cis_eqtl$CHROM[i]
    lead_snp_pos = cis_eqtl$POS[i]
    eQTL_rank = cis_eqtl$rank[i]
    eQTL_INPUT <- eQTL_INPUT_1
    nominal_pval_thres <- cis_eqtl$nominal_pval_thres[i]
    
    if(cis_eqtl$eQTL_group[i]=="B"){
      eQTL_INPUT <- eQTL_INPUT_2 %>% 
        dplyr::filter(., phenotype_id==gene_id & eSNP_position==lead_snp_pos)
    }
    
    GWAShit <- locusbreaker[locusbreaker$chr == lead_snp_chr & locusbreaker$locusStart <= lead_snp_pos & locusbreaker$locusEnd >= lead_snp_pos, ]
    
    if(nrow(GWAShit) > 0){
      GWAShit_chr <- GWAShit$chr[1]
      GWAShit_pos <- GWAShit$pos[1]
      locusStart <- GWAShit$locusStart[1]
      locusEnd <- GWAShit$locusEnd[1]
      
      eQTL_INPUT$POS <- as.numeric(eQTL_INPUT$POS)
      GWAS_INPUT$pos <- as.numeric(GWAS_INPUT$pos)
      
      out <- coloc_between_eqtl_and_gwas (eQTL_INPUT, GWAS_INPUT, GWAS_ID, N, s, ld_eqtl,
                                          cell_type, locusStart, locusEnd, gene_name, gene_id, 
                                          lead_snp, lead_snp_chr, lead_snp_pos,
                                          scaled_expression, mean_expression,
                                          do_plot="N",GWAShit_chr,GWAShit_pos,nominal_pval_thres)
      
      if(length(out)>1){
        coloc_all <- out$coloc_all$results %>% 
          dplyr::slice_max(SNP.PP.H4, n = 1)
        
        out_to_write <- data.frame(out[[2]],
                                   lead_H4_variant =  coloc_all$snp,
                                   lead_H4_SNP.PP.H4 =  coloc_all$SNP.PP.H4,
                                   n_variants_for_coloc=length(out$variants_for_coloc),
                                   eQTL_group=cis_eqtl$eQTL_group[i],
                                   eQTL_rank=cis_eqtl$rank[i],
                                   gwas_inwindow=1,
                                   status="tested for coloc")

        print(out_to_write)
        
        write.table(out_to_write, 
                    OUT_FILE,
                    sep="\t", quote=F,row.names = F, col.names = F, append = T)
      }
      if(length(out)==1){
        ## this happens when < 100 snps overlapping, SNP positions are different between GWAS and eQTL, when there's no common variants with p-value < 1e-5, when there's no common variants with eQTL p-value < nominal p-val threshold, alleles between gwas and eQTL not match
        out_to_write <- data.frame(cell_type = cell_type,
                                   gene_id =  gene_id,
                                   lead_snp =  lead_snp,
                                   lead_snp_chr = lead_snp_chr,
                                   lead_snp_pos =  lead_snp_pos,
                                   gwas_hit=paste0(GWAShit_chr,"_",GWAShit_pos),
                                   nsnps =  NA,
                                   PP.H0.abf =  NA,
                                   PP.H1.abf =  NA,
                                   PP.H2.abf =  NA,
                                   PP.H3.abf =  NA,
                                   PP.H4.abf =  NA,
                                   lead_H4_variant = NA,
                                   lead_H4_SNP.PP.H4 = NA,
                                   n_variants_for_coloc =  NA,
                                   eQTL_group=cis_eqtl$eQTL_group[i],
                                   eQTL_rank=eQTL_rank,
                                   gwas_inwindow=1,
                                   status=out)
        write.table(out_to_write, 
                    OUT_FILE,
                    sep="\t", quote=F,row.names = F, col.names = F, append = T)
      }
    }
  }
  
}

