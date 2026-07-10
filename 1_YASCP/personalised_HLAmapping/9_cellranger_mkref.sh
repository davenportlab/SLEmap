#!/bin/bash

run_cellranger=/path/cellranger-7.0.0/bin/cellranger

memory="25G"
ncpus=15


DIR_MAIN="/path"
DIR_pers_refs=${DIR_MAIN}"/single_cells/02_1.personalized_mapping/pers_refs/out"
DIR_pers_mapping=${DIR_MAIN}"/single_cells/02_1.personalized_mapping/pers_mapping"

input=${DIR_MAIN}"/metadata/final.all_wgs/out/slemap_wgs_samples.txt"

while IFS=$'\t' read -r EGAN SAMPLE batch;
do	
	per_fa=${DIR_pers_refs}"/"${SAMPLE}"/"${SAMPLE}".primaryMasked_and_HLA.fa"
	per_gtf=${DIR_pers_refs}"/"${SAMPLE}"/"${SAMPLE}".primaryMasked_and_HLA.Ensembl98.gtf"
  
  new_dir=${DIR_pers_refs}"/../Ensembl98/"${SAMPLE}

  mkref="sh ./mkref.cellranger_v7.sh "${SAMPLE}" "${per_fa}" "${per_gtf}" "${new_dir}

	bjob_1=$(echo "bsub -J mkref.$SAMPLE -q normal
    -R \"select[mem>$memory] rusage[mem=$memory]\" -M$memory
                -n$ncpus -R \"span[hosts=1]\"
                -o ../../../logs/mkref.$SAMPLE.out.log
                -e ../../../logs/mkref.$SAMPLE.err.log
                 '$mkref'");
  echo $bjob_1;

done < "$input"

