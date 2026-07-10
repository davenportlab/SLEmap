
source("/path/single_cells/scripts/data_load.R")

DIR_personalizedHLA="/path/HLApm/personalisedHLAmapping/"
source( paste0(DIR_personalizedHLA, "/scripts/load_ref.R") )
source( paste0(DIR_personalizedHLA, "/scripts/align_and_adjust_annotation.R") )
source( paste0(DIR_personalizedHLA, "/scripts/make_personalized_HLA_ref.R") )

hla_to_mask <- read.table( paste0(DIR_personalizedHLA, 
                                  "/data/references/hg38/Ensembl98.primary.12_HLA_genes.bed") )



input_data <- paste0(DIR_single_cell, 
                     "/02_1.personalized_mapping/input_HLA_alleles.311_samples.txt")
input_alleles <- read.table(input_data, header=T)
input_alleles$HLA_allele <- paste0("HLA-", input_alleles$HLA_allele)
build_personalized_HLA_ref(input_alleles, 
                           output_directory=paste0(DIR_single_cell, 
                                                   "/02_1.personalized_mapping/pers_refs") )

