DIR_MAIN="/path"
DIR_single_cell=paste0(DIR_MAIN, "/single_cells")
DIR_wgs=paste0(DIR_MAIN, "/wgs")
DIR_pers_refs=paste0(DIR_single_cell, "/02_1.personalized_mapping/pers_refs")
DIR_pers_mapping=paste0(DIR_single_cell, "/02_1.personalized_mapping/pers_mapping")
DIR_pers_fastq=paste0(DIR_single_cell, "/02_1.personalized_mapping/pers_fastq")

cellranger_7="/path/cellranger-7.0.0/bin/cellranger"

Reference_Path="/path/downloaded_from_10X/refdata-gex-GRCh38-2020-A"
ref_genome=paste0(Reference_Path, "/fasta/genome.fa")
ref_genes=paste0(Reference_Path, "/genes/genes.gtf")

meta <- fread("/path/metadata/SLEmap_metadata.csv") %>%
  as.data.frame()

