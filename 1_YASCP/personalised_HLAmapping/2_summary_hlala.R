
library(data.table)
library(ggplot2)
library(dplyr)

DIR_MAIN="/path/wgs"
DIR_summary <- paste0(DIR_MAIN, "/02_1.hla_typing/summary.HLA_LA")
DIR_OUT= paste0(DIR_MAIN, "/02_1.hla_typing/HLA_LA")

input= fread( paste0(DIR_MAIN, 
                     "/../metadata/final.all_wgs/out/3.20240716.slemap_wgs.final.txt") ) %>%
  filter(., final_set=="Y") %>%
  as.data.frame()

##################################################################
### summarize results
##################################################################

aveCoverage <- c()
final_alleles <- c()

for(sn in input$SAMPLE){

  print(sn)
  out <- read.table( paste0(DIR_OUT, "/", sn, "/hla/R1_bestguess_G.txt"), header=T)
  
  aveCoverage <- rbind(aveCoverage,
                       data.frame(sample=sn,
                                  gene=out$Locus,
                                  AverageCoverage=out$AverageCoverage))
  final_alleles <- rbind(final_alleles,
                         data.frame(sample=sn,
                                    out[, 1:3]))
}

dim(final_alleles)

write.table(final_alleles, 
            paste0(DIR_summary, "/hlala.txt"), sep="\t", quote=F, row.names = F)
write.table(aveCoverage, 
            paste0(DIR_summary, "/hlala_aveCoverage.txt"), sep="\t", quote=F, row.names = F)


