#!/bin/bash

mkdir -p $4
cd $4

/path/cellranger-7.0.0/bin/cellranger \
	count --id=$1 \
	--fastqs=$2 \
	--sample=$1 \
	--transcriptome=$3 \
	--chemistry=SC5P-R2 \
	--expect-cells=10000 \
	--localcores=20 \
	--localmem=80
cd -


