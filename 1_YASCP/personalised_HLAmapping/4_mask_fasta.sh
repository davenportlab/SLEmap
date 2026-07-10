# primary fasta 
## use `bedtools maskfasta`

ref_genome="/path/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
DIR_OUT="/path/single_cells/02_1.personalized_mapping/pers_refs" 

bedtools maskfasta -fi ${ref_genome} \
	-bed ${DIR_OUT}"/masked/Ensembl98.primary.to_mask.bed" \
	-fo ${DIR_OUT}"/masked/GRCh38.primary.HLA_masked.fa"


ref_genome="/path/GRCh38.Ensembl113/GRCh38.primary_assembly.genome.fa.gz"
DIR_OUT="/path/single_cells/02_1.personalized_mapping/pers_refs" 

bedtools maskfasta -fi ${ref_genome} \
	-bed ${DIR_OUT}"/masked/Ensembl98.primary.to_mask.bed" \
	-fo ${DIR_OUT}"/masked/GRCh38.primary.HLA_masked.fa"
