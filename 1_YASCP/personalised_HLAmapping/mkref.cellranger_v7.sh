#!/bin/bash

# $1 : genome - sample_name,  used to name output folder
# $2: fasta - primaryMasked_and_HLA.fa
# $3: gtf - primaryMasked_and_HLA.Ensembl98.gtf
# $4: out directory

mkdir -p $4
cd $4

/path/cellranger-7.0.0/bin/cellranger \
	mkref \
	--nthreads=15 \
	--genome=$1 \
	--fasta=$2 \
	--genes=$3
cd -

