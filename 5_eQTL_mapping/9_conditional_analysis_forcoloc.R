### get conditionally independent eQTL mapping results from the sig results outputted from tensorQTL 
library(tidyr)
library(dplyr)
library(data.table)

celltype <- commandArgs(trailingOnly = TRUE)
tensordir <- "/path/eQTLresults/"

##get num optim exp PCs
optimPCs <- read.csv("/path/optim_pcs.csv")
optimPCs <- optimPCs$PCs_use[optimPCs$Celltype == celltype]

## get genes with more than 1 independent eQTLs
sig_indep_all <- read.table(paste0(tensordir,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
sig_indep_gene <- sig_indep_all[sig_indep_all$rank ==2,]$phenotype_id

skips <- c()
for(k in 1:length(sig_indep_gene)){
  gene <- sig_indep_gene[k]
  print(paste0(k,"-",gene))
  
  indep_snp <- sig_indep_all[sig_indep_all$phenotype_id == gene,]$variant_id
  snps_to_test <- read.table(paste0(tensordir,"/2_indep_coloc/",celltype,"/genotypes_tidy/",gene,"_genotypes_tidy.txt"),sep="\t",header=T,row.names=1)
  exp_covariates <- read.table(paste0(tensordir,"/2_indep_coloc/",celltype,"/model/",gene,"_model.txt"),sep="\t",header=T,row.names=1)
  
  for(indep_snp_use in indep_snp){
    print(paste0(k,"-",gene,"-",indep_snp_use))
    regressed_snps <- indep_snp[!indep_snp %in% indep_snp_use] #this is all the independent snps except for the one we are testing for = snps we are regressing out
    exp_covariates_use <- exp_covariates[,!(names(exp_covariates) == indep_snp_use)] #remove column that includes the independent variant
    
    output <- data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("celltype","phenotype_id","independent_variant","regressed_variants","variant_id", "slope", "slope_se","tval","pval_nominal"))))
    
    for(snp in colnames(snps_to_test)){
      if(!(snp %in% regressed_snps)){
          exp_covariates_use_tmp <- transform(merge(exp_covariates_use,snps_to_test[snp],by=0), row.names=Row.names, Row.names=NULL) #add snp to test
          null_model <- lm(exp_covariates_use_tmp)
          if(snp %in% row.names(summary(null_model)$coefficients)){
            output[nrow(output) + 1,] = c(celltype,gene,indep_snp_use,paste(indep_snp[!indep_snp %in% indep_snp_use], collapse = "__"),snp,summary(null_model)$coefficients[snp, ])
          }else{
            skips <- append(skips,paste0(k,"-",gene,"-",indep_snp_use,"-",snp))
          }
        }
      }
    }
    write.csv(output,paste0(tensordir,"/2_indep_coloc/",celltype,"/nominal_p/",gene,"_regress_",paste0(substr(regressed_snps, 1, nchar(regressed_snps) - 4), collapse="_"),".csv"))
}


