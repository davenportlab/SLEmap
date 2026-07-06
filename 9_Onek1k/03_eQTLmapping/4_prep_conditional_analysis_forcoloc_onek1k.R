### get conditionally independent eQTL mapping results from the sig results outputted from tensorQTL - for onek1k

## before running this, /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Genotypes/byCHR/split_chr.sh has been run to chop up the vcf file by chromosome - done
## also make directories for outputs by running /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Colocalisation/indep_coloc/makedirectories.sh
# to run this file with each cell type in parallel: bash /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/12_Onek1k/4_prep_conditional_analysis_forcoloc_onek1k.sh

library(tidyr)
library(dplyr)
library(data.table)

## get indep results all for onek1k
# allindep <- data.frame()
# for(celltype in c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","Nonclassical_Monocytes","MAIT_and_GammaDelta_T_cells")){
#   indepresults <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
#   indepresults$celltype <- celltype
#   allindep <- rbind(allindep,indepresults)
# }
# 
# ##add all cells
# indepresults <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",header=T)
# indepresults$celltype <- "All"
# allindep <- rbind(allindep,indepresults)
# write.csv(allindep, "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Cis_eqtls_independent_allcelltypes_wbulklike.csv")

#for testing
# celltype <- "CD56Bright_NK_cells"
# j <- "ENSG00000160213"
# i <- 21

# start here
celltype <- commandArgs(trailingOnly = TRUE)
if(celltype == "All"){
  tensordir <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/"
}else{
  tensordir <- paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/eQTLmapping/eQTL_mapping_manualPCs/",celltype)
}
resultdir <- paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Colocalisation/indep_coloc/",celltype)

###get genes with more than one indep signal
sig_indep_all <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Cis_eqtls_independent_allcelltypes_wbulklike.csv",row.names=1)
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
    system2("bash", args = c("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/12_Onek1k/SNPs/genotype_multisnps_onek1k.sh",i,celltype,j))
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



###check that all files were made
# sig_indep_all <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Cis_eqtls_independent_allcelltypes_wbulklike.csv",row.names=1)
# sig_indep <- sig_indep_all[sig_indep_all$rank ==2,]
# 
# numfiles_made <- data.frame()
# for(celltype in c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","Nonclassical_Monocytes","MAIT_and_GammaDelta_T_cells","All")){
#   resultdir <- paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/12_Onek1k/Colocalisation/indep_coloc/",celltype)
#   resultsdir_files <- list.files(paste0(resultdir,"/model/"))
#   sig_indep_use <- sig_indep[sig_indep$celltype == celltype,]
#   numfiles_made <- rbind(numfiles_made,c(celltype,length(resultsdir_files),nrow(sig_indep_use)))
# }
# colnames(numfiles_made) <- c("celltype","numfiles","numgenes_withindep")
# all(numfiles_made$numfiles == numfiles_made$numgenes_withindep)
