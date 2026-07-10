## This script is to recalculate 3MAD thresholds (the YASCP pipeline has done 3MAD on both above and below for total_counts and n_genes_by_counts but we only want below.) 3MAD is calculated after the hard filter.  
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(UpSetR)
library(ggpubr)
library(schard)

## get obs object from totalVI results
obs = schard::h5ad2data.frame('/path/YASCPresults.h5ad','obs')
obs$pool <- do.call(rbind, strsplit(obs$experiment_id, "_"))[,3]

## 0. make empty df to store newly calculated 3MAD
QCforfiltering <- data.frame(matrix(ncol = 0, nrow = nrow(obs)))
QCforfiltering$index <- obs$X_index
QCforfiltering$experiment_id <- obs$experiment_id
QCforfiltering$pool <- obs$pool

## 1. hardfilter- see distribution of cells not passing hard filter
obs$hard_filter_manual <- "PASShardfilter"
obs$hard_filter_manual[(obs$n_genes_by_counts < 500) | (obs$log10_ngenes_by_count < 0.8) | (obs$pct_counts_gene_group__mito_transcript > 20)] <- "NOTPASShardfilter"
QCforfiltering$hard_filter_manual <- obs$hard_filter_manual[match(QCforfiltering$index,obs$X_index)]
QCforfiltering$barcode <- paste0(sapply(strsplit(QCforfiltering$index, "-"), function(x) x[1]),"-1-",QCforfiltering$pool)
write.csv(QCforfiltering[,c("index","hard_filter_manual","barcode")],"/path/QC_hardfilters.csv")

obs <- obs[obs$hard_filter_manual == "PASShardfilter",]

## 2. make empty df to store newly calculated 3MAD
QCforfiltering <- data.frame(matrix(ncol = 0, nrow = nrow(obs)))
QCforfiltering$index <- obs$X_index
QCforfiltering$experiment_id <- obs$experiment_id
QCforfiltering$pool <- obs$pool

### 3-1. 3MAD-total_counts
QCforfiltering$total_counts <- obs$total_counts[match(QCforfiltering$index,obs$X_index)]
QCforfiltering$Low_total_counts <- NA
for(i in unique(QCforfiltering$pool)){
  mad <- as.numeric(mad(QCforfiltering$total_counts[QCforfiltering$pool == i],constant=1))
  low <- median(QCforfiltering$total_counts[QCforfiltering$pool == i]) - 3*mad
  QCforfiltering$Low_total_counts[QCforfiltering$pool == i] <- as.numeric(low)
}

### 3-2. 3MAD-n_genes_by_counts
QCforfiltering$n_genes_by_counts <- obs$n_genes_by_counts[match(QCforfiltering$index,obs$X_index)]
QCforfiltering$Low_n_genes_by_counts <- NA
for(i in unique(QCforfiltering$pool)){
  mad <- as.numeric(mad(QCforfiltering$n_genes_by_counts[QCforfiltering$pool == i],constant=1))
  low <- median(QCforfiltering$n_genes_by_counts[QCforfiltering$pool == i]) - 3*mad
  QCforfiltering$Low_n_genes_by_counts[QCforfiltering$pool == i] <- as.numeric(low)
}

### 3-3. 3MAD-pct_counts_gene_group__mito_transcript
QCforfiltering$pct_counts_gene_group__mito_transcript <- obs$pct_counts_gene_group__mito_transcript[match(QCforfiltering$index,obs$X_index)]
QCforfiltering$High_pct_counts_gene_group__mito_transcript <- NA
for(i in unique(QCforfiltering$pool)){
  mad <- as.numeric(mad(QCforfiltering$pct_counts_gene_group__mito_transcript[QCforfiltering$pool == i],constant=1))
  low <- median(QCforfiltering$pct_counts_gene_group__mito_transcript[QCforfiltering$pool == i]) + 3*mad
  QCforfiltering$High_pct_counts_gene_group__mito_transcript[QCforfiltering$pool == i] <- as.numeric(low)
}

### 3-4. Getting 3MAD together
QCforfiltering$PASS3MAD <- "NOTPASS"
QCforfiltering$PASS3MAD[(QCforfiltering$total_counts > QCforfiltering$Low_total_counts) & (QCforfiltering$n_genes_by_counts > QCforfiltering$Low_n_genes_by_counts) & (QCforfiltering$pct_counts_gene_group__mito_transcript < QCforfiltering$High_pct_counts_gene_group__mito_transcript)] <- "PASS"
write.csv(QCforfiltering,"/path/QC_3MADfilters_bypool.csv")

### 4. doublets
yascp_input <- read.table('/path/yascp_inputs.tsv',sep='\t',header=T)

doublet_list <- list()
for(i in 1:nrow(yascp_input)){
  pool <- yascp_input[i,1]
  data <- read.table(paste0('/path/doublets/',pool,'__doublet_results_combined.tsv'), sep = '\t', header = TRUE)
  data$index <- paste0(data$barcodes,"-",pool,"__donor")
  doublet_list[[i]] <- data
}

doublet_list <- bind_rows(doublet_list)
doublet_list <- doublet_list[doublet_list$index %in% obs$X_index,]

### 5. doublets in atleast 3 methods
doubletmethods = c("scDblFinder_DropletType","scds_DropletType","DoubletFinder_DropletType","Scrublet_DropletType")
SLEmap_doublet <- doublet_list[,colnames(doublet_list) %in% doubletmethods]
row.names(SLEmap_doublet) <- doublet_list$index
SLEmap_doublet[SLEmap_doublet == "singlet"] <- 0
SLEmap_doublet[SLEmap_doublet == "doublet"] <- 1
SLEmap_doublet <- mutate_all(SLEmap_doublet, function(x) as.numeric(as.character(x)))
#upset(SLEmap_doublet, order.by = "freq")

SLEmap_doublet$sum <- rowSums(SLEmap_doublet)
SLEmap_doublet$over3doublet <- "singlet"
SLEmap_doublet$over3doublet[SLEmap_doublet$sum >2] <- "doublet"
write.csv(SLEmap_doublet,"/path/QC_doublets.csv")

QCforfiltering$over3doublet <- SLEmap_doublet$over3doublet[match(QCforfiltering$index,row.names(SLEmap_doublet))]

### VDJ doublets
vdjdoub <- read.csv("/path/tcr_bcr_doublet_barcodes.csv")
vdjdoub$matchID <- paste0(vdjdoub$barcode,vdjdoub$pool)

vdj_doublets <- QCforfiltering[,c("index","pool")]
vdj_doublets$matchID <- paste0(sapply(strsplit(vdj_doublets$index, "-"), function(x) x[1]),"-1",sapply(strsplit(vdj_doublets$index, "_"), function(x) x[3]))
vdj_doublets$vdjfilter <- "NOTvdjdoublet"
vdj_doublets$vdjfilter[vdj_doublets$matchID %in% vdjdoub$matchID] <- "vdjdoublet"
QCforfiltering$vdjfilter <- vdj_doublets$vdjfilter[match(QCforfiltering$index,vdj_doublets$index)]

###save object for filtering
QCforfiltering$allpass <- "FILTER"
QCforfiltering$allpass[(QCforfiltering$PASS3MAD == "PASS") & (QCforfiltering$over3doublet == "singlet") & (QCforfiltering$vdjfilter == "NOTvdjdoublet")] <- "KEEP"

QCforfiltering$barcode <- paste0(sapply(strsplit(QCforfiltering$index, "-"), function(x) x[1]),"-1-",QCforfiltering$pool)
write.csv(QCforfiltering[,c("index","PASS3MAD","over3doublet","vdjfilter","barcode")],"/path/QC_allfilters.csv")
