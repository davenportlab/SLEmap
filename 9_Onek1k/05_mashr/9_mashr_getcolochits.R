###get mash results for colocalised signals only 

# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(flashier)
library(argparse)
library(rhdf5)
library(zoo)
set.seed(123)
library(ggrepel)

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

coloc <- read.csv(paste0("/path/onek1k_locus_breaker_coloc/","1_csvfiles/SLEmap_Onek1kgroup.csv"),row.names=1)
coloc$gene_variant <- paste0(coloc$gene_id,"_",gsub("chr","",coloc$lead_H4_variant),"_",coloc$lead_H4_variant_GTF_REF,"_",coloc$lead_H4_variant_GTF_ALT)

coloc_mash <- data.frame()
for(celltype in c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","Memory_B_cells","Naive_B_cells","All")){
  print(celltype)
  outdir = paste0("/path/onek1k_locus_breaker_coloc/mashr/output/",celltype)
  input_h5 = paste0("/path/onek1k_locus_breaker_coloc/mashr/input/fastqtl_to_mash_output/",celltype,"/",celltype,".h5")
  
  ###get beta and se 
  h5file <- H5Fopen(input_h5)
  overview = h5ls(h5file)
  genes = overview$name[grep("ENSG", overview$name)]
  
  mash_per_gene = lapply(genes, function(x){
    h5_mashr_per_gene(f=input_h5, g=x, max.missing="0", replace_beta = "0", replace_beta_se = "1",verbose=F) ##nothings actually being replaced, if snp is missing in either dataset, it is dropped
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
  
  ####get the mash model
  m <- readRDS(paste0(outdir, "/model.rds"))
  
  ###get random hits
  random.subset = sample(1:nrow(beta), "200000")
  temp.random = mash_set_data(beta[random.subset,], se[random.subset,])
  Vhat = estimate_null_correlation_simple(temp.random)
  
  if(all(is.na(beta[rownames(beta) %in% coloc$gene_variant[coloc$cell_type == celltype],])) == F){
    ###get coloc hits - to predict
    beta_pred_celltype = beta[rownames(beta) %in% coloc$gene_variant[coloc$cell_type == celltype],]
    se_pred_celltype = se[rownames(se) %in% coloc$gene_variant[coloc$cell_type == celltype],]
    
    if(length(beta_pred_celltype) == 2){
      ##theres a bug where it doesnt work when predicting for one snp (length 2) - so duplicating the snp and taking just the first result.
      beta_pred_celltype <- t(as.matrix(beta_pred_celltype))
      se_pred_celltype   <- t(as.matrix(se_pred_celltype))
      rownames(beta_pred_celltype) <- rownames(beta)[rownames(beta) %in% coloc$gene_variant[coloc$cell_type == celltype]]
      rownames(se_pred_celltype) <- rownames(se)[rownames(se) %in% coloc$gene_variant[coloc$cell_type == celltype]]
      beta_pred_celltype <- rbind(beta_pred_celltype, beta_pred_celltype)
      se_pred_celltype   <- rbind(se_pred_celltype,   se_pred_celltype)
      
      ###predict
      mpred_celltype = mash_set_data(beta_pred_celltype, se_pred_celltype, V=Vhat)
      m2_celltype= mash(mpred_celltype, g=get_fitted_g(m), fixg=TRUE)
      pm  <- data.frame(get_pm(m2_celltype))[1, , drop = FALSE]
      lfsr <- data.frame(get_lfsr(m2_celltype))[1, , drop = FALSE]
      
    }else{
      beta_pred_celltype <- as.matrix(beta_pred_celltype)
      se_pred_celltype   <- as.matrix(se_pred_celltype)
      ###predict
      mpred_celltype = mash_set_data(beta_pred_celltype, se_pred_celltype, V=Vhat)
      m2_celltype= mash(mpred_celltype, g=get_fitted_g(m), fixg=TRUE)
      lfsr <- data.frame(get_lfsr(m2_celltype))
      pm <- data.frame(get_pm(m2_celltype))
    }

    
    colnames(pm) <- c("onek1k_pp","slemapEUR_pp")
    colnames(lfsr) <- c("onek1k_lfsr","slemapEUR_lfsr")
    pm$onek1k_lfsr <- lfsr$onek1k_lfsr[match(rownames(pm),rownames(lfsr))]
    pm$slemapEUR_lfsr <- lfsr$slemapEUR_lfsr[match(rownames(pm),rownames(lfsr))]
    pm$celltype <- celltype
    pm$chr_pos <- row.names(pm)
    coloc_mash <- rbind(coloc_mash,pm)
  }
}

write.csv(coloc_mash,"/path/onek1k_locus_breaker_coloc/mashr/output_coloc/mashr_coloc.csv")  ##will not have hits where the lead pph4 snp is not in onek1k dataset

###look at results
coloc_mash <- read.csv("/path/onek1k_locus_breaker_coloc/mashr/output_coloc/mashr_coloc.csv",row.names=1)
coloc <- read.csv(paste0("/path/onek1k_locus_breaker_coloc/","1_csvfiles/SLEmap_Onek1kgroup.csv"),row.names=1)
coloc$gene_variant <- paste0(coloc$gene_id,"_",gsub("chr","",coloc$lead_H4_variant),"_",coloc$lead_H4_variant_GTF_REF,"_",coloc$lead_H4_variant_GTF_ALT)
coloc_testedinonek1k <- coloc[coloc$group_onek1k != "Not tested for eQTL",] 
coloc_testedinonek1k[!paste0(coloc_testedinonek1k$cell_type,coloc_testedinonek1k$gene_variant) %in% paste0(coloc_mash$celltype, coloc_mash$chr_pos),] ## these lead coloc snps were not tested in onek1k - might be maf difference - only two hits

ggplot(coloc_mash, aes(x=slemapEUR_pp, y=onek1k_pp)) +theme_classic()+ geom_hline(yintercept=0)+ geom_vline(xintercept=0) + geom_point()+ geom_abline(intercept = 0, slope = 1,linetype="dashed")
ggplot(coloc_mash, aes(x=onek1k_pp, y=slemapEUR_pp)) +theme_classic()+ geom_hline(yintercept=0)+ geom_vline(xintercept=0) + geom_point()+ geom_abline(intercept = 0, slope = 1,linetype="dashed")

coloc_mash$eQTL <- NA
coloc_mash$eQTL[(coloc_mash$slemapEUR_pp/coloc_mash$onek1k_pp > 2)] <- "larger effect in SLEmap"
coloc_mash$eQTL[(coloc_mash$slemapEUR_pp/coloc_mash$onek1k_pp < 0.5)] <- "smaller effect in SLEmap"
coloc_mash$eQTL[(coloc_mash$slemapEUR_pp/coloc_mash$onek1k_pp >= 0.5)& (coloc_mash$slemapEUR_pp/coloc_mash$onek1k_pp <= 2)] <- "Shared effect"
coloc_mash$eQTL[(coloc_mash$onek1k_pp < 0)& (coloc_mash$slemapEUR_pp > 0)] <- "Opposite effect"
coloc_mash$eQTL[(coloc_mash$onek1k_pp > 0)& (coloc_mash$slemapEUR_pp < 0)] <- "Opposite effect"
coloc_mash$eQTL[(coloc_mash$slemapEUR_pp/coloc_mash$onek1k_pp >= 0.5)& (coloc_mash$slemapEUR_pp/coloc_mash$onek1k_pp <= 2)] <- "Shared effect"

table(coloc_mash$eQTL,useNA="always")
ggplot(coloc_mash, aes(x=onek1k_pp, y=slemapEUR_pp, color=eQTL)) + geom_point(size=1) +xlab("Onek1k posterior effect size")+ylab("SLEmap EUR posterior effect size")+theme_classic()+geom_hline(yintercept=0)+geom_vline(xintercept=0)+ geom_abline(intercept = 0, slope = 0.5,linetype="dashed",alpha=0.3)+ geom_abline(intercept = 0, slope = 2,linetype="dashed",alpha=0.3)
ggsave(paste0("/path/onek1k_locus_breaker_coloc/mashr/output_coloc/coloc_onek1k_slemapEUR.pdf"),width=7,height=5) 

ggplot(coloc_mash, aes(x=onek1k_pp, y=slemapEUR_pp, color=gsub("_"," ",celltype))) + geom_point(size=1) +xlab("OneK1K eQTL posterior effect size")+ylab("SLEmap (European ancestry only) eQTL posterior effect size")+theme_classic()+geom_hline(yintercept=0)+geom_vline(xintercept=0)+ geom_abline(intercept = 0, slope = 0.5,linetype="dashed",alpha=0.3)+ geom_abline(intercept = 0, slope = 2,linetype="dashed",alpha=0.3)
ggsave(paste0("/path/onek1k_locus_breaker_coloc/mashr/output_coloc/coloc_onek1k_slemapEUR_celltype.pdf"),width=7,height=5) 

coloc_mash[coloc_mash$eQTL != "Shared effect",]

coloc_mash$gene <-sub("_.*", "", coloc_mash$chr_pos)

gene_symbol <- read.csv("/path/ensemblID_to_genesymbol.csv")
coloc_mash$gene_symbol <- gene_symbol$gene_symbols[match(coloc_mash$gene,gene_symbol$X)]
coloc_mash[coloc_mash$eQTL == "larger effect in SLEmap",]
coloc_mash$SLEenhanced[coloc_mash$eQTL == "larger effect in SLEmap"] <- coloc_mash$gene_symbol[coloc_mash$eQTL == "larger effect in SLEmap"]

coloc_mash$group <- coloc$group_onek1k_forplot[match(paste0(coloc_mash$celltype, coloc_mash$chr_pos),paste0(coloc$cell_type,coloc$gene_variant))]

coloc_mash$group <- factor(coloc_mash$group, levels = c("Colocalized in both SLEmap and Onek1k", "No eQTL detected in Onek1k", "Colocalised only in SLEmap\n(eQTL in both)"))
coloc_mash <- coloc_mash[order(coloc_mash$group), ]
ggplot(coloc_mash, aes(x=onek1k_pp, y=slemapEUR_pp, color=group))+geom_hline(yintercept=0)+geom_vline(xintercept=0) + geom_point(size=1.4,alpha=0.75) +xlab("Onek1k posterior effect size")+ylab("SLEmap_EUR posterior effect size")+theme_classic()+ geom_abline(intercept = 0, slope = 0.5,linetype="dashed",alpha=0.3)+ geom_abline(intercept = 0, slope = 2,linetype="dashed",alpha=0.3)+scale_color_manual(values = c("Colocalized in both SLEmap and Onek1k" = "#00468B", "No eQTL detected in Onek1k" = "#D55E00","Colocalised only in SLEmap\n(eQTL in both)" = "#00C9A7"))
ggsave(paste0("/path/onek1k_locus_breaker_coloc/mashr/output_coloc/coloc_onek1k_slemapEUR_group_forpaper.pdf"),width=7.5,height=4.5) 