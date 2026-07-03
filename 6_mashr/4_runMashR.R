###mashr 
##run using: while IFS=$'\t' read -r n1 n2; do bsub -o logs/mashr_${n1}_${n2}.o -e logs/mashr_${n1}_${n2}.e -R"select[mem>60000] rusage[mem=60000]" -M60000 Rscript /path/mash/4_runMashR.R $n1 $n2 done < /path/mashresults/pairs.tsv

# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(flashier)
library(argparse)
library(rhdf5)
library(zoo)
set.seed(123)

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
celltype1 <- args_commandline[1]
celltype2 <- args_commandline[2]

args=list(input_h5 = paste0("/path/mashresults/input/fastqtl_to_mash_output/",celltype1,"_",celltype2,"/merged_test_conditions_",celltype1,"_",celltype2,".h5"), 
          max_missing="0", ###if missing in one or the other cell type
          replace_beta = "0", 
          replace_se = "1", 
          qval_f = "/path/mashresults/input/all_qval.tsv", 
          lfsr_thres=0.05,
          factor=0.5,
          use_complete_random = "FALSE", 
          nrand = "200000", 
          outdir = paste0("/path/mashresults/output/",celltype1,"_",celltype2), 
          matrices = "pca,canonical,flash", 
          reference="None", 
          dd_matrix_version="mashr")
          
if(file.exists(paste0(args$outdir)) == F){
  dir.create(args$outdir)
}

######### Prepping data ###########
# load in the complete merged dataframe matrix
setwd("/path/mashresults")
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
all_qval = read.table(args$qval_f,header=T)
all_qval$gene_variant = paste0(all_qval$phenotype_id, "_", all_qval$variant_id)

# *******BUG FIX: 'chr' from variant_ids are lost from h5 prep. Means we fail to test ******** #
if(grepl("_chr", all_qval$gene_variant[1])){
  all_qval$gene_variant = gsub("_chr", "_", all_qval$gene_variant)
}

#### make strong set with significant eQTL in either celltype - warning: if the eSNP is missing in the other cell type, it will not be picked up.
all_qval <- all_qval[all_qval$cell_type %in% c(celltype1,celltype2),]
all_qval_strong <- all_qval[all_qval$qval < 0.05,] ## this can include duplicates in variants
strong_variant <- all_qval_strong$gene_variant[!duplicated(all_qval_strong$gene_variant)]
beta_strong = beta[rownames(beta) %in% strong_variant,]
se_strong = se[rownames(se) %in% strong_variant,]
strong.subset = mash_set_data(beta_strong, se_strong)
writeLines(paste0("Out of ", length(strong_variant), " eSNP significant in either ",celltype1, " or ", celltype2, ", missing info for SNP in other cell type in ",sum(!strong_variant %in% rownames(beta))," (total gene count = ",length(genes),")"), paste0(args$outdir, "/missingleadSNPs.txt")) ####this seems larger because it includes genes not tested in both cell types. if run next time, count after genes not tested in both cell types are removed.

# Intersect beta/se, for the variants that were top per condition - in the qval file: this is the set of variants x QTLs that we want to predict for - doing separately for cell type 1 and cell type 2 - snp with lowest pval for all genes
beta_pred_celltype1 = beta[rownames(beta) %in% all_qval$gene_variant[all_qval$cell_type == celltype1],]
se_pred_celltype1 = se[rownames(se) %in% all_qval$gene_variant[all_qval$cell_type == celltype1],]

beta_pred_celltype2 = beta[rownames(beta) %in% all_qval$gene_variant[all_qval$cell_type == celltype2],]
se_pred_celltype2 = se[rownames(se) %in% all_qval$gene_variant[all_qval$cell_type == celltype2],]

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
} else {
  # Apply ultimate deconvolution
  # TO DO
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
mpred_celltype1 = mash_set_data(beta_pred_celltype1, se_pred_celltype1, V=Vhat)
m2_celltype1 = mash(mpred_celltype1, g=get_fitted_g(m), fixg=TRUE)

mpred_celltype2 = mash_set_data(beta_pred_celltype2, se_pred_celltype2, V=Vhat)
m2_celltype2 = mash(mpred_celltype2, g=get_fitted_g(m), fixg=TRUE)

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

##from celltype1 to celltype 2
celltype1_sig <- all_qval_strong$gene_variant[all_qval_strong$cell_type == celltype1]
celltype1_sig_tested <- celltype1_sig[celltype1_sig %in% rownames(m2_celltype1$result$PosteriorMean)]
celltype1_shared <- get_shared_signals(m2_celltype1)
prop_celltype1_to2 <- length(intersect(celltype1_shared,celltype1_sig_tested))/length(celltype1_sig_tested)

lfsr <- get_lfsr(m2_celltype1)
pm <- get_pm(m2_celltype1)
psd <- get_psd(m2_celltype1)
groups = colnames(lfsr)
for(i in seq_along(groups)){
  df = data.frame("posterior_means" = pm[,i], "posterior_sd" = psd[,i], "lfsr" = lfsr[,i])
  fname = gsub(".tsv", "", groups[i])
  fname = gsub("merged_", "", fname)
  write.table(df, paste0(plotout, "/predicted_results/",celltype1,"_to_", fname, "_mashr.txt"))
}

##from celltype2 to celltype1
celltype2_sig <- all_qval_strong$gene_variant[all_qval_strong$cell_type == celltype2]
celltype2_sig_tested <- celltype2_sig[celltype2_sig %in% rownames(m2_celltype2$result$PosteriorMean)]
celltype2_shared <- get_shared_signals(m2_celltype2)
prop_celltype2_to1 <- length(intersect(celltype2_shared,celltype2_sig_tested))/length(celltype2_sig_tested)

lfsr <- get_lfsr(m2_celltype2)
pm <- get_pm(m2_celltype2)
psd <- get_psd(m2_celltype2)
groups = colnames(lfsr)
for(i in seq_along(groups)){
  df = data.frame("posterior_means" = pm[,i], "posterior_sd" = psd[,i], "lfsr" = lfsr[,i])
  fname = gsub(".tsv", "", groups[i])
  fname = gsub("merged_", "", fname)
  write.table(df, paste0(plotout, "/predicted_results/",celltype2,"_to_", fname, "_mashr.txt"))
}

###save prop sharing
df <- data.frame(from = c(celltype1,celltype2),to=c(celltype2,celltype1),prop = c(prop_celltype1_to2,prop_celltype2_to1))
write.csv(df,paste0(args$outdir,"/",celltype1,"_",celltype2,"_prop_sharing.csv"),row.names=F)
