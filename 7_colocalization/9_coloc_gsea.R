####pathway analysis on coloc genes

# Load libraries
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

GWAS_ID="GCST90270940"
DIR_MAIN="/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/23_locus_breaker_coloc/outputs"

coloc <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_sigresults_",GWAS_ID,"_checksigeQTL.csv"))

coloc_count <- as.data.frame(table(coloc$cell_type))
celltypes <- as.character(coloc_count$Var1[coloc_count$Freq > 9]) ###taking cell types with at least 10 colocs

#gene_group <- read.csv(paste0(DIR_MAIN,"/1_csvfiles/coloc_genegroup_",GWAS_ID,"_checksigeQTL.csv"),row.names=1)

gsea <- function(){
  
  allresults <- data.frame()
  
  for (celltype in celltypes){
    
    # list of colocalized genes
    genes <- coloc$gene_id[coloc$cell_type == celltype]
    
    for(background in c("testedforcoloc","sigeQTL","allgenestested")){
      
      if(background == "testedforcoloc"){
        use <- read.table(paste0(DIR_MAIN,"/",GWAS_ID,"/all_checksigeQTL_checkallele/",celltype,".txt"),fill=T,header=T,sep="\t")
        background_genes <- use$gene_id[use$status == "tested for coloc"]
        
      }else if(background == "sigeQTL"){
        if(celltype == "All"){
          use <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",header=T)
          background_genes <- use$phenotype_id
        }else{
          use <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
          background_genes <- use$phenotype_id
        }
        
      }else if(background == "allgenestested"){
        if(celltype == "All"){
          use <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv",header=T,fill=T)
          background_genes <- use$phenotype_id
        }else{
          use <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T,fill=T)
          background_genes <- use$phenotype_id
        }
      }
      
      # Convert both foreground and background to Entrez IDs
      entrez_foreground <- bitr(genes,
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
      
      entrez_background <- bitr(background_genes,
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
      
      # Run ORA with custom background
      reactome_results <- enrichPathway(gene         = entrez_foreground$ENTREZID,
                                        universe     = entrez_background$ENTREZID,
                                        organism     = "human",
                                        pvalueCutoff = 1,  # keep all for filtering
                                        readable     = TRUE,
                                        minGSSize = 0)
      
      # Convert IDs back to gene symbols for labels
      reactome_results <- setReadable(reactome_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      reactome_results <- reactome_results@result
      reactome_results$celltype <- celltype
      reactome_results$background <- background
      allresults <- rbind(allresults,reactome_results)
      
    }
    
  }
  return(allresults)
}

gsea <- gsea()

write.csv(gsea,paste0(DIR_MAIN,"/1_csvfiles/coloc_gsea_",GWAS_ID,".csv"))

tested_for_coloc <- gsea[gsea$background == "testedforcoloc",]
sigeQTL <- gsea[gsea$background == "sigeQTL",]
sigeQTL <- sigeQTL[sigeQTL$Count >2,]

###run by cell lineage
coloc_count <- as.data.frame(table(coloc$cellgroup))
cellgroups <- as.character(coloc_count$Var1[coloc_count$Freq > 5]) ###dropping other T

gsea_bylineage <- function(){
  
  allresults <- data.frame()
  
  for (cellgroup in cellgroups){
    
    # list of colocalized genes
    genes <- coloc$gene_id[coloc$cellgroup == cellgroup]
    
    #get celltypes in cellgroup
    celltypes <- unique(coloc$cell_type[coloc$cellgroup == cellgroup])
    
    for(background in c("testedforcoloc","sigeQTL")){
      
    background_genes <- c()  
      if(background == "testedforcoloc"){
        for(celltype in celltypes){
          use <- read.table(paste0(DIR_MAIN,"/",GWAS_ID,"/all_checksigeQTL_checkallele/",celltype,".txt"),fill=T,header=T,sep="\t")
          background_genes <- append(background_genes,use$gene_id[use$status == "tested for coloc"])
        }
        
      }else if(background == "sigeQTL"){
        if(cellgroup == "All"){
          use <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",header=T)
          background_genes <- use$phenotype_id
        }else{
          for(celltype in celltypes){
            use <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
            background_genes <- append(background_genes,use$phenotype_id)
          }
          
        }
        
      }
    background_genes <- background_genes[!duplicated(background_genes)]
    genes <- genes[!duplicated(genes)]
      # Convert both foreground and background to Entrez IDs
      entrez_foreground <- bitr(genes,
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
      
      entrez_background <- bitr(background_genes,
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
      
      # Run ORA with custom background
      reactome_results <- enrichPathway(gene         = entrez_foreground$ENTREZID,
                                        universe     = entrez_background$ENTREZID,
                                        organism     = "human",
                                        pvalueCutoff = 1,  # keep all for filtering
                                        readable     = TRUE,minGSSize = 0)
      
      # Convert IDs back to gene symbols for labels
      reactome_results <- setReadable(reactome_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      reactome_results <- reactome_results@result
      reactome_results$cellgroup <- cellgroup
      reactome_results$background <- background
      allresults <- rbind(allresults,reactome_results)
      
    }
    
  }
  return(allresults)
}
gsea_bylineage <- gsea_bylineage()

reactome_results@result <- reactome_results@result[reactome_results@result$Count > 2,]
cnetplot(reactome_results)

write.csv(gsea_bylineage,paste0(DIR_MAIN,"/1_csvfiles/coloc_gsea_bylineage_",GWAS_ID,".csv"))
tested_for_coloc <- gsea_bylineage[gsea_bylineage$background == "testedforcoloc",]
tested_for_coloc <- tested_for_coloc[tested_for_coloc$Count >2,]
sigeQTL <- gsea_bylineage[gsea_bylineage$background == "sigeQTL",]
sigeQTL <- sigeQTL[sigeQTL$Count >2,]

### run with GO
gsea_GO_bylineage <- function(){
  
  allresults <- data.frame()
  
  for (cellgroup in cellgroups){
    
    # list of colocalized genes
    genes <- coloc$gene_id[coloc$cellgroup == cellgroup]
    
    #get celltypes in cellgroup
    celltypes <- unique(coloc$cell_type[coloc$cellgroup == cellgroup])
    
    for(background in c("testedforcoloc","sigeQTL")){
      
      background_genes <- c()  
      if(background == "testedforcoloc"){
        for(celltype in celltypes){
          use <- read.table(paste0(DIR_MAIN,"/",GWAS_ID,"/all_checksigeQTL_checkallele/",celltype,".txt"),fill=T,header=T,sep="\t")
          background_genes <- append(background_genes,use$gene_id[use$status == "tested for coloc"])
        }
        
      }else if(background == "sigeQTL"){
        if(cellgroup == "All"){
          use <- read.table("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsoft_allSNP_allcells/results/TensorQTL_eQTLS/dMean__All_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv",header=T)
          background_genes <- use$phenotype_id
        }else{
          for(celltype in celltypes){
            use <- read.table(paste0("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/5_eQTL_SLEmap_manualPCs_X_1.6_fixsort_allSNP/",celltype,"/results/TensorQTL_eQTLS/dMean__",celltype,"_all/OPTIM_pcs/base_output/base/Cis_eqtls_independent.tsv"),header=T)
            background_genes <- append(background_genes,use$phenotype_id)
          }
          
        }
        
      }
      background_genes <- background_genes[!duplicated(background_genes)]
      genes <- genes[!duplicated(genes)]
      # Convert both foreground and background to Entrez IDs
      entrez_foreground <- bitr(genes,
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
      
      entrez_background <- bitr(background_genes,
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
      
      # Run ORA with custom background
      ego <- enrichGO(gene          = entrez_foreground$ENTREZID,
                      universe      = entrez_background$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      pvalueCutoff = 1,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      readable      = TRUE,qvalueCutoff=1,minGSSize = 0,maxGSSize = 100)
      
      # Convert IDs back to gene symbols for labels
      ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      ego <- ego@result
      ego$cellgroup <- cellgroup
      ego$background <- background
      allresults <- rbind(allresults,ego)
      
    }
    
  }
  return(allresults)
}
gsea_GO_bylineage <- gsea_GO_bylineage()
write.csv(gsea_GO_bylineage,paste0(DIR_MAIN,"/1_csvfiles/coloc_gsea_GO_bylineage_",GWAS_ID,".csv"))

gsea_GO_bylineage <- gsea_GO_bylineage[gsea_GO_bylineage$Count > 2,]
