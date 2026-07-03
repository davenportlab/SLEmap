### get conditionally independent eQTL mapping results from the sig results outputted from tensorQTL 

## before running this make index file for vcf file, chop up the vcf file by chromosome, make directories for outputs 
# run this file with each cell type in parallel

library(tidyr)
library(dplyr)
library(data.table)

celltype <- commandArgs(trailingOnly = TRUE)
tensordir <- "/path/eQTLresults/"

## get genes with more than 1 indep results
sig_indep_all <- read.table(paste0(tensordir,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
sig_indep <- sig_indep_all[sig_indep_all$rank ==2,] #get genes with at more than one independent result
sig_indep <- sig_indep %>%
  separate(variant_id, into = c("chr", NA), sep = "_", extra = "drop", remove = FALSE)
sig_indep$chr <- substr(sig_indep$chr, 4, nchar(sig_indep$chr))

## get mean expression
norm_data_dir <- paste0(tensordir,celltype,"/results/norm_data/dMean__",celltype,"_all/")
mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)

##input covariates
covariates <- read.table(paste0(tensordir,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)

sampleIDs <- as.character(covariates[1,])
row.names(covariates) <- paste0(covariates$V1,covariates$V2)
covariates <- covariates[-1,-c(1:2)]    
colnames(covariates) <- sampleIDs[1:ncol(covariates)]
covariates <- mutate_all(as.data.frame(t(covariates)), function(x) as.numeric(as.character(x)))
covariates_names <- colnames(covariates)

dir.create(paste0(tensordir,"/2_indep_coloc/",celltype,"/snps"), showWarnings = FALSE)
dir.create(paste0(tensordir,"/2_indep_coloc/",celltype,"/genotypes"), showWarnings = FALSE)
dir.create(paste0(tensordir,"/2_indep_coloc/",celltype,"/genotypes_tidy"), showWarnings = FALSE)
dir.create(paste0(tensordir,"/2_indep_coloc/",celltype,"/model"), showWarnings = FALSE)
dir.create(paste0(tensordir,"/2_indep_coloc/",celltype,"/nominal_p"), showWarnings = FALSE)

## get snps tested by chromosome from tensorQTL results 
skips <- c()
for (i in unique(sig_indep$chr)){
  sig_indep_chr <- sig_indep[sig_indep$chr == i,]
  sig_indep_chr_genes <- unique(sig_indep_chr$phenotype_id)
  sig_indep_chr_genes_snps <- fread(paste0(tensordir,celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.",i,".tsv"), select = c("phenotype_id", "variant_id"))[phenotype_id %in% sig_indep_chr_genes]
  
  ## for each gene in indep results, 
  for (j in sig_indep_chr_genes){
    #save list of snps for that gene
    sig_indep_chr_genes_snps_for1gene <- sig_indep_chr_genes_snps$variant_id[sig_indep_chr_genes_snps$phenotype_id == j]
    write.table(sig_indep_chr_genes_snps_for1gene,row.names=F,col.names=F,quote=F, file=paste0(tensordir,"/2_indep_coloc/",celltype,"/snps/",j,"_snps.txt"))
    
    #get genotypes
    system2("bash", args = c("/path/genotype_multisnps.sh",i,celltype,j))
    indiv <- read.table(paste0(tensordir,"/2_indep_coloc/",celltype,"/genotypes/",j,"_genotypes.012.indv"),sep="\t")
    genotype <- read.table(paste0(tensordir,"/2_indep_coloc/",celltype,"/genotypes/",j,"_genotypes.012.gz"),sep="\t",row.names = 1)
    row.names(genotype) <- indiv$V1
    colnames(genotype) <- sig_indep_chr_genes_snps_for1gene
    
    #leave only genotype from individuals included in the eQTL mapping
    genotype <- genotype[row.names(genotype) %in% row.names(covariates),]
    
    ##tensorQTL imputes genotypes for missing values: they use the mean genotype (not integer)
    for(k in sig_indep_chr_genes_snps_for1gene){genotype[[k]][genotype[[k]] == -1] <- mean(genotype[[k]][genotype[[k]] != -1],na.rm=T)}
    write.table(genotype,paste0(tensordir,"/2_indep_coloc/",celltype,"/genotypes_tidy/",j,"_genotypes_tidy.txt"),sep="\t", col.names=NA)
    
    #get all independent eQTLs for the gene
    indep_snp <- sig_indep_all[sig_indep_all$phenotype_id == j,]$variant_id
    
    #get gene expression
    mean_use <- as.data.frame(t(mean[j,]))
    mean_use$sampleID <- paste0("SLE_",row.names(mean_use) %>% strsplit( ".",fixed=T) %>%  sapply( "[", 3 ) %>% strsplit( "_",fixed=T) %>%  sapply( "[", 3 ))
    mean_use <- merge(mean_use,genotype[colnames(genotype) %in% indep_snp],by.x="sampleID",by.y=0)
    
    #get covariates
    mean_use <- merge(mean_use,covariates,by.x="sampleID",by.y=0)
    
    ## write file with all components of model
    write.table(mean_use,paste0(tensordir,"/2_indep_coloc/",celltype,"/model/",j,"_model.txt"),sep="\t",row.names = F)
  }
  
}




