library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(ggpubr)

coloc <- read.csv("/path/coloc/outputs/1_csvfiles/coloc_sigresults.csv")

for(i in 1:nrow(coloc)){
  gg_genetracks(locus(xrange=c(coloc$locusStart[i],coloc$locusEnd[i]), seqname=coloc$lead_snp_chr[i],chrom="chr",pos="pos",ens_db = "EnsDb.Hsapiens.v86"), filter_gene_biotype = 'protein_coding',highlight=coloc$gene_symbol[i],highlight_col = "#F9A825", gene_col = "#90A4AE",exon_col = "#90A4AE",exon_border = "#90A4AE",)
  ggsave(paste0("/path/coloc/outputs/0_plots/locuszoom/",coloc$gene_symbol[i],"_",coloc$gwas_hit[i],".pdf"),height=2,width=6)
}

for(i in 1:nrow(coloc)){
  gg_genetracks(locus(xrange=c(coloc$locusStart[i],coloc$locusEnd[i]), seqname=coloc$lead_snp_chr[i],chrom="chr",pos="pos",ens_db = "EnsDb.Hsapiens.v86"),highlight=coloc$gene_symbol[i],highlight_col = "#F9A825", gene_col = "#90A4AE",exon_col = "#90A4AE",exon_border = "#90A4AE",)
  ggsave(paste0("/path/coloc/outputs/0_plots/locuszoom_incnoncoding/",coloc$gene_symbol[i],"_",coloc$gwas_hit[i],".pdf"),height=2,width=6)
}
