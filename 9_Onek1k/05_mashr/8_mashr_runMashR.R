###prep for mashr 
# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(flashier)
library(argparse)
library(rhdf5)
library(zoo)
set.seed(123)

###allcells added to onek1k qval file
qval_onek1k <- read.table("/path/onek1k_locus_breaker_coloc/mashr/input/all_qval_onek1k.tsv",header=T,fill=T,sep="\t")
colnames(qval_onek1k)[colnames(qval_onek1k) == "CD56Bright_NK_cells"] <- "cell_type"
all_cells <- fread("/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv",header=T,fill=T,sep="\t")
all_cells$cell_type <- "All"
qval_onek1k <- rbind(qval_onek1k,all_cells)
write.table(qval_onek1k,"/path/onek1k_locus_breaker_coloc/mashr/input/all_qval_onek1k.tsv",row.names = F,sep="\t")


### mash
## we are doing a pairwise comparison, hence if snp is missing in one cell type, we are dropping it and not filling it in.
h5_mashr_per_gene <- function(f, gene, max.missing=Inf, replace_beta, replace_beta_se, verbose=T) {
  temp = h5read(f, paste0("/", gene))
  rownames <- temp$rownames
  # add gene
  rownames = paste0(gene, "_", rownames)
  colnames <- temp$colnames
  
  betas <- t(temp$beta)
  rownames(betas) <- rownames
  colnames(betas) <- colnames
  
  error <- t(temp$se)
  row.names(error) <- rownames
  colnames(error) <- colnames
  
  
  if(any(!complete.cases(betas))){
    nrow_start <- nrow(betas)
    keep <- rowSums(is.na(error)) <= max.missing & rowSums(is.na(betas)) <= max.missing
    error <- error[keep, ]
    betas <- betas[keep, ]
    if(verbose){
      message(match(gene,genes),"-Dropping ", nrow_start-nrow(betas), " because # NAs > ", 
              max.missing, ". ", nrow(betas), " tests remaining.")
      
      message("Have ", signif(100*sum(is.na(betas))/length(betas),2), "% missing values. Filling...")
    }
    if(replace_beta  == "median"){
      betas_median <- median(betas, na.rm=TRUE)
      if(verbose){
        message("Replacing remaining beta NAs with median betas: ", betas_median)
      }
      betas[is.na(betas)] <- betas_median
    } else{
      betas[is.na(betas)] <- as.numeric(replace_beta)
    }
    
    if(replace_beta_se  == "median"){
      error_median <- median(error, na.rm=TRUE)
      if(verbose){
        message("Replacing remaining beta NAs with median error: ", error_median)
      }
      error[is.na(error)] <- error_median
    } else {
      error[is.na(error)] <- as.numeric(replace_beta_se)
    }
    
  }
  
  mashData <- mash_set_data(as.matrix(betas), 
                            as.matrix(error))
  return(mashData)
}

args_commandline <- commandArgs(trailingOnly = TRUE)
celltype <- args_commandline[1]

args=list(input_h5 = paste0("/path/onek1k_locus_breaker_coloc/mashr/input/fastqtl_to_mash_output/",celltype,"/",celltype,".h5"), 
          max_missing="0", ###if missing in one or the other cell type
          replace_beta = "0", 
          replace_se = "1", 
          qval_slemap = "/path/onek1k_locus_breaker_coloc/mashr/input/all_qval_slemapEUR.tsv",
          qval_onek1k = "/path/onek1k_locus_breaker_coloc/mashr/input/all_qval_onek1k.tsv",
          lfsr_thres=0.05,
          factor=0.5,
          use_complete_random = "FALSE", 
          nrand = "200000", 
          outdir = paste0("/path/onek1k_locus_breaker_coloc/mashr/output/",celltype), 
          matrices = "pca,canonical,flash", 
          reference="None", 
          dd_matrix_version="mashr")
          
if(file.exists(paste0(args$outdir)) == F){
  dir.create(args$outdir)
}

######### Prepping data ###########
# load in the complete merged dataframe matrix
setwd("/path/onek1k_locus_breaker_coloc/mashr")
print(paste0("Reading from: ", args$input_h5))
h5file <- H5Fopen(args$input_h5)
overview = h5ls(h5file)
genes = overview$name[grep("ENSG", overview$name)]
mash_per_gene = lapply(genes, function(x){
  h5_mashr_per_gene(f=args$input_h5, g=x, max.missing=args$max_missing, replace_beta=args$replace_beta, replace_beta_se=args$replace_se)
})

# Combine into a single object
beta_all = lapply(mash_per_gene, function(x){
  return(x$Bhat)
})
beta = do.call(rbind, beta_all)
se_all = lapply(mash_per_gene, function(x){
  return(x$Shat)
})
se = do.call(rbind, se_all)

######## For the prediction of effects across conditions, we want to find the top gene-snp pair per condition ########
# Get all qvalues
all_qval_slemap = read.table(args$qval_slemap,header=T,fill=T,sep="\t")
colnames(all_qval_slemap)[ncol(all_qval_slemap)] <- "cell_type"
all_qval_slemap$gene_variant = paste0(all_qval_slemap$phenotype_id, "_", all_qval_slemap$variant_id)

all_qval_onek1k = read.table(args$qval_onek1k,header=T,fill=T,sep="\t")
all_qval_onek1k$gene_variant = paste0(all_qval_onek1k$phenotype_id, "_", all_qval_onek1k$variant_id)


# *******BUG FIX: 'chr' from variant_ids are lost from h5 prep. Means we fail to test ******** #
if(grepl("_chr", all_qval_slemap$gene_variant[1])){
  all_qval_slemap$gene_variant = gsub("_chr", "_", all_qval_slemap$gene_variant)
}
if(grepl("_chr", all_qval_onek1k$gene_variant[1])){
  all_qval_onek1k$gene_variant = gsub("_chr", "_", all_qval_onek1k$gene_variant)
}

#### make strong set with significant eQTL in either slemap or onek1k - warning: if the eSNP is missing in the other cell type, it will not be picked up.
all_qval <- rbind(all_qval_onek1k,all_qval_slemap)
all_qval <- all_qval[all_qval$cell_type %in% c(celltype),]
all_qval_strong <- all_qval[all_qval$qval < 0.05,] ## this can include duplicates in variants
strong_variant <- all_qval_strong$gene_variant[!duplicated(all_qval_strong$gene_variant)]
beta_strong = beta[rownames(beta) %in% strong_variant,]
se_strong = se[rownames(se) %in% strong_variant,]
strong.subset = mash_set_data(beta_strong, se_strong)

# Intersect beta/se, for the variants that were top per condition - in the qval file: this is the set of variants x QTLs that we want to predict for - snp with lowest pval for all genes
beta_pred_celltype = beta[rownames(beta) %in% all_qval_slemap$gene_variant[all_qval_slemap$cell_type == celltype],]
se_pred_celltype = se[rownames(se) %in% all_qval_slemap$gene_variant[all_qval_slemap$cell_type == celltype],]

######## Get random subset from the entire dataset ########
random.subset = sample(1:nrow(beta), args$nrand)
temp.random = mash_set_data(beta[random.subset,], se[random.subset,])

######### Correlation structure ###########
# Estimating correlation structure from the random tests
Vhat = estimate_null_correlation_simple(temp.random)

######### Set up proper strong/random subsets using the correlation structure ########
data.random = mash_set_data(beta[random.subset,], se[random.subset,], V=Vhat)
data.strong=mash_set_data(beta_strong, se_strong,V=Vhat)

# Update reference if desired
if( args$reference != "None"){
  data.random <- mash_update_data(data.random, ref = args$reference)
  data.strong <- mash_update_data(data.strong, ref = args$reference)
}

######### Derive data-driven covariance matrices ########
if(args$dd_matrix_version == "mashr"){
  want_matrices = unlist(strsplit(args$matrices, ","))
  dd_matrices = NULL
  if("pca" %in% want_matrices){
    U.pca = cov_pca(data.strong,2)
    dd_matrices = c(U.pca)
  }
  if("flash" %in% want_matrices){
    U.f = cov_flash(data.strong, factors="nonneg", tag="non_neg")
    if(length(dd_matrices) == 0){
      dd_matrices[[1]] = U.f
    } else {
      dd_matrices = c(dd_matrices, U.f)
    }
  }
  # Apply extreme deconvolution 
  U.ed = cov_ed(data.strong, dd_matrices)
}

########## Now fit the mash model (data-driven + canonical) to random tests ##########
if("canonical" %in% want_matrices){
  U.c = cov_canonical(data.random)
  Ulist = c(dd_matrices, U.c)
} else {
  Ulist = c(dd_matrices)
}
m = mash(data.random, Ulist = Ulist, outputlevel = 1)

# Save
saveRDS(m, paste0(args$outdir, "/model.rds"))

########## Get posterior summaries #######
# Predict for the predicted set - now incorporating correlations
mpred_celltype = mash_set_data(beta_pred_celltype, se_pred_celltype, V=Vhat)
m2_celltype= mash(mpred_celltype, g=get_fitted_g(m), fixg=TRUE)

########## Make files to summarise the model #######
if(file.exists(paste0(args$outdir, "/model_summary")) == F){
  dir.create(paste0(args$outdir, "/model_summary"))
  dir.create(paste0(args$outdir, "/model_summary/predicted_results"))
}
plotout = paste0(args$outdir, "/model_summary")

get_shared_signals <- function(mash_model, lfsr_thres = args$lfsr_thres, factor = args$factor) {
  pm   <- get_pm(mash_model)
  lfsr <- get_lfsr(mash_model)
  
  idx <- which(
    lfsr[, 1] < lfsr_thres &
    lfsr[, 2] < lfsr_thres &
    sign(pm[, 1]) == sign(pm[, 2]) &
    abs(pm[, 1] - pm[, 2]) < factor * (abs(pm[, 1]) + abs(pm[, 2])) / 2
  )
  rownames(pm)[idx]
}

lfsr <- get_lfsr(m2_celltype)
pm <- get_pm(m2_celltype)
psd <- get_psd(m2_celltype)
groups = colnames(lfsr)
for(i in seq_along(groups)){
  df = data.frame("posterior_means" = pm[,i], "posterior_sd" = psd[,i], "lfsr" = lfsr[,i])
  fname = gsub(".tsv", "", groups[i])
  fname = gsub("merged_", "", fname)
  write.table(df, paste0(plotout, "/predicted_results/",celltype,"_",fname,"_mashr.txt"))
}

pm <- data.frame(pm)
lfsr <- data.frame(lfsr)
colnames(pm) <- c("onek1k_pp","slemapEUR_pp")
colnames(lfsr) <- c("onek1k_lfsr","slemapEUR_lfsr")

pm$onek1k_lfsr <- lfsr$onek1k_lfsr[match(rownames(pm),rownames(lfsr))]
pm$slemapEUR_lfsr <- lfsr$slemapEUR_lfsr[match(rownames(pm),rownames(lfsr))]
pm$eQTL <- NA
pm$eQTL[(pm$onek1k_lfsr > 0.05) & (pm$slemapEUR_lfsr > 0.05)] <- "Not significant in either"
pm$eQTL[(pm$onek1k_lfsr < 0.05) & (pm$slemapEUR_lfsr > 0.05)] <- "Only significant in Onek1k"
pm$eQTL[(pm$onek1k_lfsr > 0.05) & (pm$slemapEUR_lfsr < 0.05)] <- "Only significant in SLEmap"
pm$eQTL[(pm$onek1k_lfsr < 0.05) & (pm$slemapEUR_lfsr < 0.05) &(pm$slemapEUR_pp/pm$onek1k_pp > 2)] <- "Significant in both, larger effect in SLEmap"
pm$eQTL[(pm$onek1k_lfsr < 0.05) & (pm$slemapEUR_lfsr < 0.05) &(pm$slemapEUR_pp/pm$onek1k_pp < 0.5)] <- "Significant in both, smaller effect in SLEmap"
pm$eQTL[(pm$onek1k_lfsr < 0.05) & (pm$slemapEUR_lfsr < 0.05) & (pm$onek1k_pp < 0)& (pm$slemapEUR_pp > 0)] <- "Opposite effect"
pm$eQTL[(pm$onek1k_lfsr < 0.05) & (pm$slemapEUR_lfsr < 0.05) & (pm$onek1k_pp > 0)& (pm$slemapEUR_pp < 0)] <- "Opposite effect"
pm$eQTL[(pm$onek1k_lfsr < 0.05) & (pm$slemapEUR_lfsr < 0.05) & (pm$slemapEUR_pp/pm$onek1k_pp >= 0.5)& (pm$slemapEUR_pp/pm$onek1k_pp <= 2)] <- "Shared effect"
table(pm$eQTL,useNA="always")
ggplot(pm, aes(x=slemapEUR_pp, y=onek1k_pp, color=eQTL)) + geom_point(size=1) +xlab("SLEmap_EUR posterior effect size")+ylab("Onek1k posterior effect size")+theme_classic()+geom_hline(yintercept=0)+geom_vline(xintercept=0)+ggtitle(paste0(celltype,"- n SLEenhanced: ",nrow(pm[pm$eQTL == "Significant in both, larger effect in SLEmap" & is.na(pm$eQTL) == F,])))
ggsave(paste0(plotout, "/predicted_results/",celltype,"_groups.pdf"),width=9,height=5)


pm$gene <- sub("_.*", "", rownames(pm))
gene_symbol <- read.csv("/path/ensemblID_to_genesymbol.csv")
pm$gene_symbol <- gene_symbol$gene_symbols[match(pm$gene,gene_symbol$X)]

write.csv(pm,paste0(plotout, "/predicted_results/",celltype,"_grouped.csv"))

