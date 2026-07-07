### get conditionally independent eQTL mapping results from the sig results outputted from tensorQTL - for onek1k
## before running this, chop up the vcf file by chromosome, make directories for outputs
#run in parallel through cell types

library(tidyr)
library(dplyr)
library(data.table)

# start here
celltype <- commandArgs(trailingOnly = TRUE)
if(celltype == "All"){
  tensordir <- "/path/onek1k_eQTLresults/ManualPCs/"
}else{
  tensordir <- paste0("/path/onek1k_eQTLresults/ManualPCs/",celltype)
}
resultdir <- paste0("/path/onek1k/Colocalisation/indep_coloc/",celltype)

###get genes with more than one indep signal
sig_indep_all <- read.csv("/path/onek1k_eQTLresults/Cis_eqtls_independent_allcelltypes_wbulklike.csv",row.names=1)
sig_indep_all <- sig_indep_all[sig_indep_all$celltype == celltype,]
sig_indep <- sig_indep_all[sig_indep_all$rank ==2,] #get genes with at more than one independent result
sig_indep <- sig_indep %>%
  separate(variant_id, into = c("chr", NA), sep = "_", extra = "drop", remove = FALSE)
sig_indep$chr <- substr(sig_indep$chr, 4, nchar(sig_indep$chr))

## get mean expression
norm_data_dir <- paste0(tensordir,"/results/norm_data/dMean__",celltype,"_all/")
mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)

##input covariates
covariates <- read.table(paste0(tensordir,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)

sampleIDs <- as.character(covariates[1,])
row.names(covariates) <- paste0(covariates$V1,covariates$V2)
covariates <- covariates[-1,-c(1:2)]    
colnames(covariates) <- sampleIDs[1:ncol(covariates)]
covariates <- mutate_all(as.data.frame(t(covariates)), function(x) as.numeric(as.character(x)))
covariates_names <- colnames(covariates)


## get snps tested by chromosome from tensorQTL results 
skips <- c()
for (i in unique(sig_indep$chr)){
  sig_indep_chr <- sig_indep[sig_indep$chr == i,]
  sig_indep_chr_genes <- unique(sig_indep_chr$phenotype_id)
  sig_indep_chr_genes_snps <- fread(paste0(tensordir,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",i,".tsv"), select = c("phenotype_id", "variant_id"))[phenotype_id %in% sig_indep_chr_genes]
  
  ## for each gene in indep results, 
  for (j in sig_indep_chr_genes){
    #save list of snps for that gene
    sig_indep_chr_genes_snps_for1gene <- sig_indep_chr_genes_snps$variant_id[sig_indep_chr_genes_snps$phenotype_id == j]
    write.table(sig_indep_chr_genes_snps_for1gene,row.names=F,col.names=F,quote=F, file=paste0(resultdir,"/snps/",j,"_snps.txt"))
    
    #get genotypes
    system2("bash", args = c("/path/genotype_multisnps_onek1k.sh",i,celltype,j))
    indiv <- read.table(paste0(resultdir,"/genotypes/",j,"_genotypes.012.indv"),sep="\t")
    genotype <- read.table(paste0(resultdir,"/genotypes/",j,"_genotypes.012.gz"),sep="\t",row.names = 1)
    row.names(genotype) <- indiv$V1
    colnames(genotype) <- sig_indep_chr_genes_snps_for1gene
    
    #leave only genotype from individuals included in the eQTL mapping
    genotype <- genotype[row.names(genotype) %in% row.names(covariates),]
    
    ##tensorQTL imputes genotypes for missing values: they use the mean genotype (not integer)
    for(k in sig_indep_chr_genes_snps_for1gene){genotype[[k]][genotype[[k]] == -1] <- mean(genotype[[k]][genotype[[k]] != -1],na.rm=T)}
    write.table(genotype,paste0(resultdir,"/genotypes_tidy/",j,"_genotypes_tidy.txt"),sep="\t", col.names=NA)
    
    #get all independent eQTLs for the gene
    indep_snp <- sig_indep_all[sig_indep_all$phenotype_id == j,]$variant_id
    
    #get gene expression
    mean_use <- as.data.frame(t(mean[j,]))
    mean_use$sampleID <- sub(".*_(\\d+_\\d+)$", "\\1", row.names(mean_use))
    mean_use <- merge(mean_use,genotype[colnames(genotype) %in% indep_snp],by.x="sampleID",by.y=0)
    
    #get covariates
    mean_use <- merge(mean_use,covariates,by.x="sampleID",by.y=0)
    
    ## write file with all components of model
    write.table(mean_use,paste0(resultdir,"/model/",j,"_model.txt"),sep="\t",row.names = F)
  }
  
}

