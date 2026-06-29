#!/bin/bash

plink \
  --bfile 281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95 \
  --extract snplist.reliable_17491.txt \
  --make-bed \
  --out 281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_reliable_snp

cd King_out/

king -b King_reference/KGref.bed,281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_reliable_snp \
    --pca --projection --rplot

# Rscript king_ancestryplot_fixed_paper.R
