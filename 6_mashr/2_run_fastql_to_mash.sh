#### run to get the random set (strong set will not be made from this output) using 2_run_fastql_to_mash.sh (scripts adapted from https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb and https://stephenslab.github.io/mashr/articles/eQTL_outline.html)

mkdir -p /path/mashresults/input/fastqtl_to_mash_output

while IFS=$'\t' read -r n1 n2; do
    echo "Processing $n1 and $n2"

    outdir="/path/mashresults/input/fastqtl_to_mash_output/${n1}_${n2}"
    mkdir -p "$outdir"

    # gene list paths
    gene_file_1="/path/mashresults/overlap_genes/${n1}_${n2}.tsv"
    gene_file_2="/path/mashresults/overlap_genes/${n2}_${n1}.tsv"

    # choose existing one
    if [[ -f "$gene_file_1" ]]; then
        gene_file="$gene_file_1"
    elif [[ -f "$gene_file_2" ]]; then
        gene_file="$gene_file_2"
    else
        echo "No gene list found for $n1 $n2, skipping"
        continue
    fi

    bsub -e logs/${n1}_${n2}.e -o logs/${n1}_${n2}.o -q normal -n 4 -M 100000 \
    -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" \
    singularity exec -B /path -B /software /path/hdf5tools.sif \
    sos run /path/mash/fastqtl_to_mash.ipynb \
        --data-list /path/mashresults/input/merged_test_conditions/merged_test_conditions_${n1}_${n2}.list \
        --gene-list "$gene_file" \
        -j 8 \
        --best-per-gene 0 \
        --random-per-gene -1 \
        --random-snp-size 200000 \
        --cwd "$outdir"

done < /path/mashresults/pairs.tsv
