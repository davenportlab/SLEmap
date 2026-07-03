temp_file=$(mktemp)
output=/path/mashresults/input/all_qval.tsv
for i in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells Cytotoxic_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells  Naive_CD8_T_cells TEMRA DN_T_cells Memory_B_cells Naive_B_cells All; do
    f=/path/eQTLresults/${i}/results/TensorQTL_eQTLS/dMean__${i}_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv
    echo $f
    if [ ! -f ${output} ]; then
        echo "Getting header"
        head $f -n 1 >> $temp_file
        awk -v value="${i}" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$temp_file" >> ${output}
    fi
    awk -v value="${i}" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$f" | tail -n +2 >> ${output}
done