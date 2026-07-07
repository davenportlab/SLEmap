#!/bin/bash

repo=/path/onek1k_locus_breaker_coloc/mashr
outdir=${repo}/input
mkdir -p $outdir

for k in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells  Naive_CD8_T_cells TEMRA Memory_B_cells Naive_B_cells; do
    echo ${k}
    inputdir=/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs/${k}/results/TensorQTL_eQTLS/dMean__${k}_all/OPTIM_pcs/base_output/base 

    fnom1=${inputdir}/cis_nominal1.cis_qtl_pairs.1.tsv
    head -n 1 "$fnom1" > ${outdir}/merged_onek1k_${k}.tsv

    for c in $(seq 1 22); do
        echo $c
        fnom=${inputdir}/cis_nominal1.cis_qtl_pairs.${c}.tsv
        tail -n +2 "$fnom" >> ${outdir}/merged_onek1k_${k}.tsv
    done

    fnom=${inputdir}/cis_nominal1.cis_qtl_pairs.X.tsv
    tail -n +2 "$fnom" >> ${outdir}/merged_onek1k_${k}.tsv
    
    gzip ${outdir}/merged_onek1k_${k}.tsv
done

for k in All; do
    echo ${k}
    inputdir=/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs_allcells/results/TensorQTL_eQTLS/dMean__${k}_all/OPTIM_pcs/base_output/base 

    fnom1=${inputdir}/cis_nominal1.cis_qtl_pairs.1.tsv
    head -n 1 "$fnom1" > ${outdir}/merged_onek1k_${k}.tsv

    for c in $(seq 1 22); do
        echo $c
        fnom=${inputdir}/cis_nominal1.cis_qtl_pairs.${c}.tsv
        tail -n +2 "$fnom" >> ${outdir}/merged_onek1k_${k}.tsv
    done

    fnom=${inputdir}/cis_nominal1.cis_qtl_pairs.X.tsv
    tail -n +2 "$fnom" >> ${outdir}/merged_onek1k_${k}.tsv
    
    gzip ${outdir}/merged_onek1k_${k}.tsv
done