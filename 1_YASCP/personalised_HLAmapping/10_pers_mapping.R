source("/path/data_load.R")

DIR_pers_refs=paste0(DIR_single_cell, "/02_1.personalized_mapping/pers_refs/Ensembl98")

DIR_demulti <- "/path/yascp/slemap"
donor_wgs_f <- paste0(DIR_demulti, "/gtmatch/assignments_all_pools.tsv")

donor_wgs <- fread(donor_wgs_f) %>% as.data.frame() %>%
  mutate(pool_donor=paste0(pool, ".", donor_query)) %>% 
  select(., donor_gt, final_panel, pool_donor, pool) %>%
  filter(., donor_gt!="NONE")
dim(donor_wgs)

run_cellranger="/path/cellranger-7.0.0/bin/cellranger"
ncpus=20

DIR_demultiplex <- "/path/yascp/slemap/handover/Donor_Quantification"

output <- paste0(DIR_single_cell, "/scripts/02_1.personalized_mapping/tmp.03_2.sh")
write.table("#", output, append = F,
            quote=F, row.names = F, col.names = F)

for(i in 1:nrow(donor_wgs)){

  print(i)
  sample_id <- paste0( donor_wgs$donor_gt[i], "_", donor_wgs$pool[i])
  fastq_path=paste0(DIR_pers_fastq, "/", sample_id)
  
  per_ref=paste0(DIR_pers_refs, "/", donor_wgs$donor_gt[i], "/", donor_wgs$donor_gt[i] )

  out_prefix=paste0(DIR_pers_mapping, "_Ensembl98/", sample_id)
  count <- paste0("./count.cellranger_v7.sh ", sample_id, " ", fastq_path, " ", per_ref, " ", out_prefix)
  print(count)
  cmd <- paste0("bsub -J count.", sample_id, " -q normal -R \"select[mem>80000] rusage[mem=80000]\" -M80000 ",
                "-n", ncpus, " -R \"span[hosts=1]\" ",
                "-o ../../logs/count.", sample_id, ".out.%J ",
                "-e ../../logs/count.", sample_id, ".error.%J '", count, "'" )
  
  write.table(cmd, output, append = T,
              quote=F, row.names = F, col.names = F)
}



