#!/bin/bash

samples=$(awk '{print $1}' slemap_cellranger_to_pool.txt)
pools=$(awk '{print $2}' slemap_cellranger_to_pool.txt)

input_fasta=tcr/${this_pool}_filtered_contig.fasta
input_annotations=tcr/${this_pool}_filtered_contig_annotations.csv

outdir=IgBLAST_output_TCR

AssignGenes.py igblast -s ${input_fasta} -b igblast_databases_nov24/igblast \
   --organism human --loci tr --format blast --outdir ${outdir} --outname ${this_pool} --exec ncbi-igblast-1.22.0/bin/igblastn


MakeDb.py igblast -i ${outdir}/${this_pool}_igblast.fmt7 -s ${input_fasta} \
   -r igblast_databases_nov24/germlines/imgt/human/vdj --10x ${input_annotations} --extended --outdir ${outdir} --outname ${this_pool}

