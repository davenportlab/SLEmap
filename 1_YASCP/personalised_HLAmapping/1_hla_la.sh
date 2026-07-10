#!/bin/bash

### for HLA-LA.pl: https://github.com/DiltheyLab/HLA-LA


### index graph PRG_MHC_GRCh38_withIMGT 
# HLA-LA --action prepareGraph  --PRG_graph_dir  /path/HLA-LA/PRG_MHC_GRCh38_withIMGT

# samtools index NA12878.mini.cram
# ./HLA-LA.pl --BAM NA12878.mini.cram --graph PRG_MHC_GRCh38_withIMGT --sampleID NA12878 --maxThreads 7


memory="55G"
ncpus=10

i=0

DIR_OUT="/path/wgs/02_1.hla_typing/HLA_LA"
DIR_cram="/path/wgs/crams"


DIR_meta="/path/metadata/final.all_wgs/out"
input=$(echo $DIR_meta/3.20240716.slemap_wgs.final.txt)

#while IFS="," read -r SLE_WGS batch; 
while IFS=$'\t' read -r SAMPLE EGAN batch version path final_set; 

do
	if [[ $final_set == "Y" ]] ; then
		
    base_name=$(echo $SAMPLE)
    filename=$(echo "$DIR_cram/$base_name.cram")
    echo ${batch}
    echo ${SAMPLE}
    
    run_hla_la=$(echo "HLA-LA.pl --BAM ${filename} \
	    --customGraphDir /path/HLA-LA/graphs \
	    --graph PRG_MHC_GRCh38_withIMGT \
	    --workingDir ${DIR_OUT} \
	    --sampleID ${base_name} \
	    --maxThreads $ncpus ")
    echo $run_hla_la;
    
    # submitting jobs
    bjob_1=$(echo "bsub -J hla_la.$base_name -q normal -G humgen-priority 
	    -R \"select[mem>$memory] rusage[mem=$memory]\" -M$memory 
	    -n$ncpus -R \"span[hosts=1]\" 
	    -o ../../logs/hla_la.$base_name.out.log  
	    -e ../../logs/hla_la.$base_name.err.log 
	    '$run_hla_la'");
	    
	    echo $bjob_1;
	fi
 
done < "$input"

