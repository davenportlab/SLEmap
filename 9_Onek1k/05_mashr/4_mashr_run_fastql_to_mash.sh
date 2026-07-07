## fastqtl_to_mash.ipynb script from mashr tutorial (from https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb and https://stephenslab.github.io/mashr/articles/eQTL_outline.html)

#make list of comparisons first
mkdir /path/onek1k_locus_breaker_coloc/mashr/input/merged_test_conditions && cd $_
for k in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells  Naive_CD8_T_cells TEMRA Memory_B_cells Naive_B_cells All; do
    printf "/path/onek1k_locus_breaker_coloc/mashr/input/merged_slemapEUR_%s.tsv.gz\n/path/onek1k_locus_breaker_coloc/mashr/input/merged_onek1k_%s.tsv.gz\n" "$k" "$k" > "${k}.list"
done


mkdir -p /path/onek1k_locus_breaker_coloc/mashr/input/fastqtl_to_mash_output

for k in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells  Naive_CD8_T_cells TEMRA Memory_B_cells Naive_B_cells All; do
    echo "Processing $k"

    outdir="/path/onek1k_locus_breaker_coloc/mashr/input/fastqtl_to_mash_output/${k}"
    mkdir -p "$outdir"

    # gene list paths
    gene_file="/path/onek1k_locus_breaker_coloc/mashr/overlap_genes/${k}.tsv"

    bsub -e logs/${k}.e -o logs/${k}.o -q normal -n 4 -M 100000 \
    -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" \
    singularity exec -B /path -B /software /path/hdf5tools.sif \
    sos run /path/24_onek1k_locus_breaker_coloc/fastqtl_to_mash.ipynb \
        --data-list /path/onek1k_locus_breaker_coloc/mashr/input/merged_test_conditions/${k}.list \
        --gene-list "$gene_file" \
        -j 8 \
        --best-per-gene 0 \
        --random-per-gene -1 \
        --random-snp-size 200000 \
        --cwd "$outdir"

done

#### get a file with all qvals from all celltypes from slemapEUR(things to test)
temp_file=$(mktemp)
output=/path/onek1k_locus_breaker_coloc/mashr/input/all_qval_slemapEUR.tsv
for i in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells Naive_CD8_T_cells TEMRA Memory_B_cells Naive_B_cells All; do
    f=/path/ancestry_coloc/EUR/ManualPCs/${i}/results/TensorQTL_eQTLS/dMean__${i}_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv
    echo $f
    if [ ! -f ${output} ]; then
        echo "Getting header"
        head $f -n 1 >> $temp_file
        awk -v value="${i}" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$temp_file" >> ${output}
    fi
    awk -v value="${i}" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$f" | tail -n +2 >> ${output}
done

#### get a file with all qvals from all celltypes from onek1k (all-pbmc will be added later)
temp_file=$(mktemp)
output=/path/onek1k_locus_breaker_coloc/mashr/input/all_qval_onek1k.tsv
for i in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells Naive_CD8_T_cells TEMRA Memory_B_cells Naive_B_cells; do
    f=/path/Onek1k/eQTLmapping/eQTL_mapping_manualPCs/${i}/results/TensorQTL_eQTLS/dMean__${i}_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv
    echo $f
    if [ ! -f ${output} ]; then
        echo "Getting header"
        head $f -n 1 >> $temp_file
        awk -v value="${i}" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$temp_file" >> ${output}
    fi
    awk -v value="${i}" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$f" | tail -n +2 >> ${output}
done
