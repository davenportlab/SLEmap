### identify examples of shared/unshared signals and plot with adjusted expression values for paper fig2- must use module with vcftools (ex.HGI/softpack/users/hj10/SLEmap_HJ/5)
# eQTL sharing framework

setwd("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/kb21/flanders")
# .libPaths("~/R/x86_64-pc-linux-gnu-library/4.1")
#library(UpSetR)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
# library(gridGraphics)
library(tibble)
library(pheatmap)
library(igraph)

#library(data.table)
#library(dplyr)

# read in files
loci <- list.files("flanders_updated_outputs/results/coloc_info_tables/", full.names = T)
finemapped.loci <- lapply(loci, fread)
finemapped.loci <- rbindlist(finemapped.loci)

coloc.results <- read.delim("flanders_updated_outputs/results/coloc/coloc_run_colocalization.table.all.tsv")

coloc.info.tables <- list.files("flanders_updated_outputs/results/coloc_info_tables/", full.names = T)
coloc.info.tables <- lapply(coloc.info.tables, fread)
coloc.info <- rbindlist(coloc.info.tables)

load("results_signal_sharing.RData")
load("results_membership_df.RData")

# list of cell types
all.cell.types <- c("Naive_CD4_T_cells","CM_CD4_T_cells","EM_CD4_T_cells",
                    "Cytotoxic_CD4_T_cells","Regulatory_CD4_T_cells","Naive_CD8_T_cells",
                    "CM_CD8_T_cells","EM_CD8_T_cells","TEMRA","DN_T_cells",
                    "Naive_B_cells","Memory_B_cells","Classical_Monocytes",
                    "CD56Bright_NK_cells","CD56Dim_NK_cells", "all_cells")

# make expression matrix i.e. is each gene expressed in each cell type
expressed.genes <- list()
gene.expression <- read.table("../../hj10/Results/5_eQTL_SLEmap_findPCs_X_1.6_fixsort_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Expression_Data.sorted.bed")
expressed.genes[["all_cells"]] <- gene.expression$V4
all_genes <- gene.expression$V4

for(i in all.cell.types){
  if(i == "all_cells"){
    print(i)
  } else {
    print(i)
    gene.expression <- read.table(paste0("../../hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",
                                         i, "/results/TensorQTL_eQTLS/dMean__", i, "_all/OPTIM_pcs/base_output/base/Expression_Data.sorted.bed"))
    gene.expression <- gene.expression$V4
    expressed.genes[[i]] <- gene.expression
  }
}

# list of all expressed genes
all_genes <- unique(unlist(expressed.genes))

# matrix of cell type expression by cell type
gene_matrix <- sapply(expressed.genes, function(genes) all_genes %in% genes)
rownames(gene_matrix) <- all_genes

# signal presence matrix: is each eQTL signal detected in each cell type (from previous script)
# Combine list of shared signals
all_gene_signal_sharing_named <- lapply(all_gene_signal_sharing, function(df) {
  df$rownames <- rownames(df)   # Save rownames as a column
  return(df)
})

# Bind rows, filling missing columns with NA
combined <- bind_rows(all_gene_signal_sharing_named)

# restore rownames
rownames(combined) <- combined$rownames
combined$rownames <- NULL  # Drop the temporary column

# Replace NA with FALSE
signal_matrix <- combined
signal_matrix[is.na(combined)] <- FALSE

# signal info (long form)
all_genes_membership_df <- rbindlist(all_membership_df)
# add short signal name key to membership table
all_genes_membership_df$gene <- unlist(strsplit(all_genes_membership_df$node, "::"))[seq(2, by=3, length.out=nrow(all_genes_membership_df))]
all_genes_membership_df$signal <- paste0(all_genes_membership_df$gene, "_", all_genes_membership_df$group)

# key matching signal name to gene
signal_to_gene <- data.frame("Signal"=rownames(signal_matrix),
                             "Gene"=substr(rownames(signal_matrix), 1, nchar(rownames(signal_matrix))-2))

#########find good examples - discordant #######
gene_symbols <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/slemap/singlecell/ensemblID_to_genesymbol.csv")

discordant <- list.files("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/kb21/flanders/discordant/convincing")
discordant <- gsub(".pdf","",discordant)
discordant <- gene_symbols$gene_symbols[match(discordant,gene_symbols$X)]
discordant

signal_matrix_onlysc <- signal_matrix[signal_matrix$all_cells == F,]
signal_matrix_onlysc$gene <- signal_to_gene$Gene[match(row.names(signal_matrix_onlysc),signal_to_gene$Signal)]
signal_matrix_onlysc$gene_symbols <- gene_symbols$gene_symbols[match(signal_matrix_onlysc$gene,gene_symbols$X)]
signal_matrix_onlysc[signal_matrix_onlysc$gene_symbols == "LY96",] # opposite direction between classical mono and naive CD8 - SLE related (part of TLR4), different signal from all cells (overlap in credible set but does not coloc) 
all_genes_membership_df[all_genes_membership_df$gene == "ENSG00000154589",] ##ol.cs.group is chr8:73989162:A:G - probably not the lead coloc snp, not sure how to get it and what to plot

#######find good examples - concordant #####
gene_symbols <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/slemap/singlecell/ensemblID_to_genesymbol.csv")

signal_matrix_B <- signal_matrix[signal_matrix$all_cells == F & signal_matrix$Naive_B_cells == T & signal_matrix$Memory_B_cells == T,]
signal_matrix_B$sig_count <- rowSums(signal_matrix_B == TRUE)
signal_matrix_B <- signal_matrix_B[signal_matrix_B$sig_count == 2,]

signal_matrix_B$gene <- signal_to_gene$Gene[match(row.names(signal_matrix_B),signal_to_gene$Signal)]
signal_matrix_B$gene_symbols <- gene_symbols$gene_symbols[match(signal_matrix_B$gene,gene_symbols$X)]
signal_matrix_B$gene_symbols

signal_matrix_B[signal_matrix_B$gene_symbols == "CD96",] # eQTL in naive and memory B cells but not in any other including all cells - interesting because CD96 is an immune check point receptor in T cells and NK cells mainly, with more expression in these cell types, but eQTL is only found in B cells
all_genes_membership_df[all_genes_membership_df$gene == "ENSG00000153283",] 

##########plot adjusted expression ############
#library(biomaRt)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#ensembl <- useEnsembl("snp",dataset = "hsapiens_snp")

## set working directory somewhere with write permission - this will be where the genotype files are made
setwd("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/3_eQTL_prep/SNPs")


getadjexp <- function(celltype,Gene,SNP,color="purple",minmax=NA){
  
  if(celltype == "All"){
    eQTLrun <- "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells"
  }else{
    eQTLrun <- paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype)
  }
  
  
  ## set this to where the QTLite run results are. There should be a norm_data directory with each cell type
  norm_data_dir <- paste0(eQTLrun,"/results/norm_data/dMean__",celltype,"_all/")
  
  ## get mean expression
  mean <- read.table(paste0(norm_data_dir,"normalised_phenotype.tsv"),sep="\t",header=T,row.names = 1)
  
  # get SNP position
  SNPs <- setNames(data.frame(unlist(strsplit(SNP, "_"))[1], unlist(strsplit(SNP, "_"))[2], unlist(strsplit(SNP, "_"))[2]), c('chr_name', 'chrom_start','chrom_end'))
  
  # ## get SNP position fom rsid 
  # SNPs <- getBM(attributes=c("refsnp_id","chr_name","chrom_start","chrom_end"),filters ="snp_filter", values =SNP, mart = ensembl, uniqueRows=TRUE)
  # SNPs$chr_name <- as.numeric(SNPs$chr_name)#sometimes there's two SNPs with the same rsID, the wrond one seems to be on a unknown chromosome name
  # SNPs <- SNPs[(SNPs$chr_name < 23) & (is.na(SNPs$chr_name) == F) ,] 
  # system2("bash", args = c("/lustre/scratch126/opentargets/opentargets/OTAR2064/working/users/hj10/Scripts/3_eQTL_prep/SNPs/genotype.sh", paste0("chr",SNPs$chr_name), SNPs$chrom_start,SNPs$chrom_end,SNP))
  
  ## get genotype
  #system2("bash", args = c("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/3_eQTL_prep/SNPs/genotype.sh", SNPs$chr_name, SNPs$chrom_start,SNPs$chrom_end,SNP))
  indiv <- read.table(paste0(SNP,".012.indv"),sep="\t")
  genotype <- read.table(paste0(SNP,".012.gz"),sep="\t")
  #system2("rm",args= c(paste0(SNP,"*")))
  indiv$geno <- genotype$V2
  
  ## plot (mean) gene exp
  mean_use <- as.data.frame(t(mean[Gene,]))
  mean_use$sampleID <- paste0("SLE_",row.names(mean_use) %>% strsplit( ".",fixed=T) %>%  sapply( "[", 3 ) %>% strsplit( "_",fixed=T) %>%  sapply( "[", 3 ))
  mean_use$geno <- indiv$geno[match(mean_use$sampleID, indiv$V1)]
  colnames(mean_use) <- c("Mean_exp","sampleID","Genotype")
  
  ##input covariates
  covariates <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_input/base_output__base/Covariates.tsv"),fill=T)
  
  sampleIDs <- as.character(covariates[1,])
  row.names(covariates) <- paste0(covariates$V1,covariates$V2)
  covariates <- covariates[-1,-c(1:2)]    
  colnames(covariates) <- sampleIDs[1:ncol(covariates)]
  covariates <- mutate_all(as.data.frame(t(covariates)), function(x) as.numeric(as.character(x)))
  
  mean_use <- merge(mean_use,covariates,by.x="sampleID",by.y=0)
  
  # Get optim phenotype PCs and 5 genotypePCs
  optimPCs <- read.table(paste0(eQTLrun,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output__base/optim_pcs.txt"))
  optimPCs <- as.numeric(optimPCs[1,1])
  
  genotypePCs <- paste0("GenotypePC",1:5)
  phenotypePCs <- paste0("PhenotypePC",1:optimPCs)
  
  ##tensorQTL imputes genotypes for missing values: they use the mean genotype (not integer)
  mean_use$Genotype_imputed <- mean_use$Genotype
  mean_use$Genotype_imputed[mean_use$Genotype_imputed == -1] <- mean(mean_use$Genotype_imputed[mean_use$Genotype_imputed != -1])
  
  # get full model: the coefficient for genotypes should match the slope from tensorQTL
  null_model <- lm(mean_use[c("Mean_exp", "Genotype_imputed", genotypePCs, phenotypePCs)]) 
  
  coefs <- null_model[["coefficients"]][append(genotypePCs,phenotypePCs)]
  
  mean_use$adjusted_expr <- mean_use$Mean_exp - rowSums(sweep(mean_use[append(genotypePCs,phenotypePCs)], 2, coefs, FUN = "*"))
  beta <- as.numeric(null_model[["coefficients"]][["Genotype_imputed"]])
  
  if(color=="purple"){
    colorpalette <- c("#756BB1","#756BB1","#756BB1")
  }else if(color=="blue"){
    colorpalette <- c("#5C6BC0","#5C6BC0","#5C6BC0")
  }
  
  if(all(is.na(minmax))==T){
      # Plot adjusted expression by genotype
      p2 <- ggplot(mean_use[mean_use$Genotype != -1,], aes(y = adjusted_expr, x = Genotype, group=Genotype,color=factor(Genotype))) + scale_x_continuous(breaks = c(0, 1, 2))+scale_color_manual(values=colorpalette) + geom_jitter(size=1,width=0.25,alpha=0.7)+ geom_boxplot(outlier.shape = NA,alpha=0)+ggtitle(paste0(celltype,Gene,"-",SNP,", n=",nrow(mean_use[mean_use$Genotype != -1,])))+theme_classic()+annotate("text", x = 2, y = min(mean_use[mean_use$Genotype != -1,]$adjusted_expr), label = paste0("b=",round(beta,3)))+ylab("")+xlab("")+ theme(legend.position = "none")

  }else{
    min <- minmax[[1]]
    max <- minmax[[2]]
    p2 <- ggplot(mean_use[mean_use$Genotype != -1,], aes(y = adjusted_expr, x = Genotype, group=Genotype,color=factor(Genotype))) + scale_x_continuous(breaks = c(0, 1, 2))+scale_color_manual(values=colorpalette) + geom_jitter(size=1,width=0.25,alpha=0.7)+ geom_boxplot(outlier.shape = NA,alpha=0)+ggtitle(paste0(celltype,Gene,"-",SNP,", n=",nrow(mean_use[mean_use$Genotype != -1,])))+theme_classic()+annotate("text", x = 2, y = min, label = paste0("b=",round(beta,3)))+ylab("")+xlab("")+ theme(legend.position = "none")+ylim(min,max)
    
  }

  
  return(list(mean_use,p2))
  
}

###ENSG00000154589 - chr8_73989162_G (ref)_A (alt)
LY96_Classical_Monocytes <- getadjexp("Classical_Monocytes","ENSG00000154589","chr8_73989162_G_A","purple",minmax=c(0,0.6))
LY96_Naive_CD8_T_cells <- getadjexp("Naive_CD8_T_cells","ENSG00000154589","chr8_73989162_G_A","purple")
LY96_Allcells <- getadjexp("All","ENSG00000154589","chr8_73989162_G_A","purple",minmax=c(0,0.2))
ggarrange(LY96_Classical_Monocytes[[2]],LY96_Naive_CD8_T_cells[[2]],LY96_Allcells[[2]],nrow=1,ncol=3,common.legend = T)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/0_plots/LY96_example.pdf",width=6.2,height=2.7)

# min <- min(min(LY96_Classical_Monocytes[[1]]$adjusted_expr),min(LY96_Naive_CD8_T_cells[[1]]$adjusted_expr),min(LY96_Allcells[[1]]$adjusted_expr))
# max <- max(max(LY96_Classical_Monocytes[[1]]$adjusted_expr),max(LY96_Naive_CD8_T_cells[[1]]$adjusted_expr),max(LY96_Allcells[[1]]$adjusted_expr))
# minmax <- c(min,max)
# 
# LY96_Classical_Monocytes <- getadjexp("Classical_Monocytes","ENSG00000154589","chr8_73989162_G_A","purple",minmax)
# LY96_Naive_CD8_T_cells <- getadjexp("Naive_CD8_T_cells","ENSG00000154589","chr8_73989162_G_A","purple",minmax)
# LY96_Allcells <- getadjexp("All","ENSG00000154589","chr8_73989162_G_A","purple",minmax)
# ggarrange(LY96_Classical_Monocytes[[2]],LY96_Naive_CD8_T_cells[[2]],LY96_Allcells[[2]],nrow=1,ncol=3,common.legend = T)


ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/0_plots/LY96_example.pdf",width=6.3,height=3)

## check that eQTL is not sig in all cells
## check nominal p-val thresholdand whether theres a sig eQTL 
allcells_nominal_qval <- fread("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv")
allcells_nominal_qval[allcells_nominal_qval$phenotype_id == "ENSG00000154589",] #qval = 0.2620144 -- how did this input to flanders??? 

allcells_nominal <- fread(paste0("grep ENSG00000154589 ","/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.8.tsv"))
colnames(allcells_nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
allcells_nominal[allcells_nominal$variant_id == "chr8_73989162_G_A",]


###ENSG00000153283 - chr3_111535372_A (ref) _C (alt)
CD96_Naive_B_cells <- getadjexp("Naive_B_cells","ENSG00000153283","chr3_111535372_A_C","blue")
CD96_Memory_B_cells <- getadjexp("Memory_B_cells","ENSG00000153283","chr3_111535372_A_C","blue")
CD96_Allcells <- getadjexp("All","ENSG00000153283","chr3_111535372_A_C","blue",minmax=c(0,1.15))
ggarrange(CD96_Naive_B_cells[[2]],CD96_Memory_B_cells[[2]],CD96_Allcells[[2]],nrow=1,ncol=3,common.legend = T)
ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/0_plots/CD96_example.pdf",width=6.2,height=2.7)

min <- min(min(CD96_Naive_B_cells[[1]]$adjusted_expr),min(CD96_Memory_B_cells[[1]]$adjusted_expr),min(CD96_Allcells[[1]]$adjusted_expr))

CD96_Naive_B_cells <- getadjexp("Naive_B_cells","ENSG00000153283","chr3_111535372_A_C","blue",min)
CD96_Memory_B_cells <- getadjexp("Memory_B_cells","ENSG00000153283","chr3_111535372_A_C","blue",min)
CD96_Allcells <- getadjexp("All","ENSG00000153283","chr3_111535372_A_C","blue",min)
ggarrange(CD96_Naive_B_cells[[2]],CD96_Memory_B_cells[[2]],CD96_Allcells[[2]],nrow=1,ncol=3,common.legend = T)

ggsave("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/0_plots/CD96_example.pdf",width=6.3,height=3)

