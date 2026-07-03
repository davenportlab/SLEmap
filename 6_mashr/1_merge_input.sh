repo=/path/mashresults
outdir=${repo}/input
mkdir -p $outdir

for k in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells Cytotoxic_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells  Naive_CD8_T_cells TEMRA DN_T_cells Memory_B_cells Naive_B_cells All; do
    echo ${k}
    inputdir=/path/eQTLresults/${k}/results/TensorQTL_eQTLS/dMean__${k}_all/OPTIM_pcs/base_output/base 

    fnom1=${inputdir}/cis_nominal1.cis_qtl_pairs.1.tsv
    head -n 1 "$fnom1" > ${outdir}/merged_${k}.tsv

    for c in $(seq 1 22); do
        echo $c
        fnom=${inputdir}/cis_nominal1.cis_qtl_pairs.${c}.tsv
        tail -n +2 "$fnom" >> ${outdir}/merged_${k}.tsv
    done

    fnom=${inputdir}/cis_nominal1.cis_qtl_pairs.X.tsv
    tail -n +2 "$fnom" >> ${outdir}/merged_${k}.tsv
    
    gzip ${outdir}/merged_${k}.tsv
done

#### get list of merged result files 
cd /path/mashresults/input/
files=("$PWD"/*.tsv.gz); for ((i=0;i<${#files[@]};i++)); do for ((j=i+1;j<${#files[@]};j++)); do f1=${files[i]}; f2=${files[j]}; n1=$(basename "$f1" .tsv.gz | sed 's/^merged_//'); n2=$(basename "$f2" .tsv.gz | sed 's/^merged_//'); printf "%s\n%s\n" "$f1" "$f2" > "/path/mashresults/input/merged_test_conditions/merged_test_conditions_${n1}_${n2}.list"; done; done

###get list of all pairs
files=("$PWD"/*.tsv.gz); for ((i=0;i<${#files[@]};i++)); do for ((j=i+1;j<${#files[@]};j++)); do n1=$(basename "${files[i]}" .tsv.gz | sed 's/^merged_//'); n2=$(basename "${files[j]}" .tsv.gz | sed 's/^merged_//'); echo -e "$n1\t$n2"; done; done > /path/mashresults/pairs.tsv