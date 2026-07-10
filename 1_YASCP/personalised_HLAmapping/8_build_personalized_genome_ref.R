library(Biostrings)
library(rtracklayer)

# cmd=paste0("R-4.1.0 --vanilla  < mkref_personalized.R --args ", sample_id, " ", HLA_alleles, " ", drb_hap)

source("/path/single_cells/scripts/data_load.R")


ref_primary <- paste0(DIR_pers_refs, "/masked/GRCh38.primary.HLA_masked.fa")
refSeq.primary.masked <- readBStringSet(ref_primary)

gtf_primary <- paste0(DIR_pers_refs, "/masked/GRCh38.primary.Ensembl98.HLA_masked.sorted.gtf")
refGTF.primary.masked <- rtracklayer::import(gtf_primary) 

samples <- list.files(path=paste0(DIR_pers_refs, "/out"), 
                      pattern = "*per.fa",
                      full.names = T, 
                      recursive = T)
samples

for(s in samples[1:100]){
  
  ind <- gsub(".*\\/out\\/(SLE.*)\\/SLE.*", "\\1", s)
  print(ind)
  
  p.fa <- s
  p.gtf <- gsub("fa", "gtf", s)
  print(p.fa)
  print(p.gtf)
  
  # reference sequences (final)
  
  refSeq.personalizedHLA <- readBStringSet(p.fa)
  refSeq.final <- c(refSeq.personalizedHLA, refSeq.primary.masked)
  Biostrings::writeXStringSet(refSeq.final, 
                              paste0(DIR_pers_refs, "/out/",
                                     ind, "/", 
                                     ind, ".primaryMasked_and_HLA.fa") )
  
  # gtf (final)
  
  refGTF.personalizedHLA <- rtracklayer::import(p.gtf) 
  refGTF.final <- c(refGTF.personalizedHLA, refGTF.primary.masked)
  rtracklayer::export(refGTF.final, 
                      paste0(DIR_pers_refs, "/out/",
                             ind, "/", 
                             ind, ".primaryMasked_and_HLA.Ensembl98.gtf"))
}
