## plot all coloced genes across all cell types - works with module load HGI/softpack/users/hj10/SLEmap_HJ/8
library(ggplot2)
library(ggpubr)
library(rlist)
library(data.table)
library(dplyr)
DIR_MAIN="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs"

## FOR GWAS
GWAS_input= "/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/14_colocalization/inputs/gwas/"
#GWAS_ID="GCST003156"
GWAS_ID="GCST90270940"
GWAS_INPUT <- fread (paste0(GWAS_input,GWAS_ID,"/",GWAS_ID,"_for_eqtl.txt.gz")) %>%
  as.data.frame() %>%
  dplyr::rename(., variant_id_v2=variant_id) %>%
  dplyr::mutate(variant_id = paste0(chr, "_", pos))
GWAS_INPUT <- GWAS_INPUT[!((GWAS_INPUT$chr == "chr6") & (GWAS_INPUT$pos > 25000000)&(GWAS_INPUT$pos < 34000000)),] 

GWAS_INPUT$pos <- as.numeric(GWAS_INPUT$pos)
GWAS_INPUT$p_value <- as.numeric(GWAS_INPUT$p_value)

locusbreaker <- read.csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/gwashit_locusbreaker_forinput.csv")
coloc_group <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_genegroup_",GWAS_ID,"_checksigeQTL.csv"),nrow=F)

##GET coloc results
coloc <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_sigresults_",GWAS_ID,"_checksigeQTL.csv"))

##plotting
plotlist_list <- list()
colocresults_inothercelltypes <- data.frame(matrix(ncol=10))
colnames(colocresults_inothercelltypes) <- c("phenotype_id","variant_id","pval_nominal","slope","slope_se","POS","chr_pos","coloc_leadsnp","celltype","independent_variant")

for(i in 1:nrow(coloc)){
  
  colocvariant <- coloc$lead_H4_variant[i]
  coloc_genesymbol <- coloc$gene_symbol[i]
  print(coloc_genesymbol)
  
  ##get GWAS
  GWAShit = coloc$gwas_hit[i]
  GWAShit_chr = sapply(strsplit(GWAShit, "_"), `[`, 1)
  GWAShit_pos = as.numeric(sapply(strsplit(GWAShit, "_"), `[`, 2))
  pos_start = locusbreaker$locusStart[locusbreaker$variant_id == GWAShit]
  pos_end = locusbreaker$locusEnd[locusbreaker$variant_id == GWAShit]
  
  GWAS_input_plot <- GWAS_INPUT[GWAS_INPUT$chr == GWAShit_chr,]
  GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos < pos_end,]
  GWAS_input_plot <- GWAS_input_plot[GWAS_input_plot$pos > pos_start,]
  
  GWAS_input_plot$color <- "notlead_colocsnp"
  GWAS_input_plot$color[GWAS_input_plot$variant_id ==colocvariant] <- "lead_colocsnp"
  GWAS_input_plot <- GWAS_input_plot %>% arrange(desc(color))
  
  plotlist <- list()
  ## get eQTL
  for(celltype in c("CD56Bright_NK_cells","CD56Dim_NK_cells","Classical_Monocytes","CM_CD4_T_cells","Cytotoxic_CD4_T_cells","EM_CD4_T_cells","Naive_CD4_T_cells","Regulatory_CD4_T_cells","CM_CD8_T_cells","EM_CD8_T_cells", "Naive_CD8_T_cells","TEMRA","DN_T_cells","Memory_B_cells","Naive_B_cells","All")){
    print(celltype)
    gene=coloc$gene_id[i]
    gene_symbol=coloc$gene_symbol[i]
    gene_chr=coloc$lead_snp_chr[i]
    gene_chr <- gsub("chr", "",gene_chr)
    eQTL_group <- coloc_group$group[(coloc_group$gene == gene)&(coloc_group$celltype == celltype)]
    
    if(!eQTL_group %in% c("Not tested for eQTL")){
      if(celltype == "All"){
        EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
      }else{
        EQTL_input=paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/")
        }
      
      independent=read.table(paste0(EQTL_input,"Cis_eqtls_independent.tsv"),header=T)
      independent_n <- nrow(independent[independent$phenotype_id == gene,])
      
      plot <- list()
      if(independent_n > 1){
        if(celltype == "All"){
          indeplist <- list.files(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/2_indep_coloc/",celltype,"/nominal_p/"),gene,full.names = TRUE)
        }else{
          indeplist <- list.files(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/2_indep_coloc/",celltype,"/nominal_p/"),gene,full.names = TRUE)
        }
        
        for(m in 1:length(indeplist)){
          indepresults <- read.csv(indeplist[m],header=T)
          indepresults$POS <- as.numeric(sapply(strsplit(indepresults$variant_id, "_"), `[`, 2))
          indepresults$chr_pos <- sapply(strsplit(indepresults$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
          
          overlapping_variants <- intersect(indepresults$chr_pos,GWAS_input_plot$variant_id)
          indepresults$overlap <- "not_overlap"
          indepresults$overlap[indepresults$chr_pos %in% overlapping_variants] <- "overlap"
          indepresults$coloc_leadsnp <- "notlead_colocsnp"
          indepresults$coloc_leadsnp[indepresults$chr_pos == colocvariant] <- "lead_colocsnp"
          indepresults <- indepresults %>% arrange(desc(coloc_leadsnp))
          
          plot <- list.append(plot,ggplot(indepresults, aes(x = POS, y = -log10(pval_nominal),color=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol,"-",eQTL_group))+theme_bw()) # to show all snps in region while having alpha higher in overlapping ones
          
          colocresults_inothercelltypes <- rbind(colocresults_inothercelltypes,indepresults[indepresults$chr_pos == colocvariant,c("phenotype_id","variant_id","pval_nominal","slope","slope_se","POS","chr_pos","coloc_leadsnp","celltype","independent_variant")])
        }
        
        GWAS_input_plot$overlap <- "not_overlap"
        GWAS_input_plot$overlap[GWAS_input_plot$variant_id %in% overlapping_variants] <- "overlap"
        plot <- list.append(plot,ggplot(GWAS_input_plot, aes(x = pos, y = -log10(p_value),color=color)) + geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw())# to show all snps in region while having alpha higher in overlapping ones
        
        
      }else{
        
        nominal <- fread(paste0("grep ",gene," ",EQTL_input,"cis_nominal1.cis_qtl_pairs.",gene_chr,".tsv"))
        colnames(nominal) <- c("phenotype_id","variant_id","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se")
        nominal <- nominal[nominal$phenotype_id == gene,]
        nominal$POS <- as.numeric(sapply(strsplit(nominal$variant_id, "_"), `[`, 2))
        nominal$chr_pos <- sapply(strsplit(nominal$variant_id, "_"), function(x) paste(x[1], x[2], sep = "_"))
        
        overlapping_variants <- intersect(nominal$chr_pos,GWAS_input_plot$variant_id)
        nominal$overlap <- "not_overlap"
        nominal$overlap[nominal$chr_pos %in% overlapping_variants] <- "overlap"
        nominal$coloc_leadsnp <- "notlead_colocsnp"
        nominal$coloc_leadsnp[nominal$chr_pos == colocvariant] <- "lead_colocsnp"
        nominal <- nominal %>% arrange(desc(coloc_leadsnp))
        
        plot <- list.append(plot,ggplot(nominal, aes(x = POS, y = -log10(pval_nominal),color=coloc_leadsnp)) +geom_point(aes(size=overlap, alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(celltype," ",coloc_genesymbol,"-",eQTL_group))+theme_bw()) # to show all snps in region while having alpha higher in overlapping ones
        
        GWAS_input_plot$overlap <- "not_overlap"
        GWAS_input_plot$overlap[GWAS_input_plot$variant_id %in% overlapping_variants] <- "overlap"
        plot <- list.append(plot,ggplot(GWAS_input_plot, aes(x = pos, y = -log10(p_value),color=color)) + geom_point(aes(size=overlap,alpha = overlap))+scale_alpha_manual(values = c("not_overlap" = 0.2, "overlap" = 0.7))+scale_size_manual(values = c("not_overlap" = 1, "overlap" = 1.5)) +scale_x_continuous(labels = scales::comma, limits=c(pos_start,pos_end))+ggtitle(paste0(GWAS_ID,": GWAS hit:",GWAShit))+theme_bw())# to show all snps in region while having alpha higher in overlapping ones
        
        independent_variant=NA
        colocresults_inothercelltypes <- rbind(colocresults_inothercelltypes,cbind(nominal[nominal$chr_pos == colocvariant,c("phenotype_id","variant_id","pval_nominal","slope","slope_se","POS","chr_pos","coloc_leadsnp")],celltype,independent_variant))
      }
      plotlist[[celltype]] <- do.call(ggarrange, c(plot, ncol = 1, align = "v",common.legend=T))
    }
  }
  plotlist_list[[i]] <- plotlist
}


names(plotlist_list) <- coloc$gene_symbol

if (!dir.exists(paste0(DIR_MAIN,"/0_plots/coloc_regional_pairs_",GWAS_ID))) {
  dir.create(paste0(DIR_MAIN,"/0_plots/coloc_regional_pairs_",GWAS_ID))
}
for (i in 1:length(plotlist_list)) {
  pdf(paste0(DIR_MAIN,"/0_plots/coloc_regional_pairs_",GWAS_ID,"/coloc_regionalplots_",names(plotlist_list)[i],".pdf"))
  print(plotlist_list[[i]])
  dev.off()
}

colocresults_inothercelltypes <- colocresults_inothercelltypes[-1,]
write.csv(colocresults_inothercelltypes,paste0(DIR_MAIN,"/1_csvfiles/colocresults_inothercelltypes_",GWAS_ID,"_checksigeQTL.csv"),row.names=F)


