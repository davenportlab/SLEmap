#!/bin/bash

module load HGI/softpack/users/hj10/SLEmap_HJ/15

#for i in CD56Bright_NK_cells CD56Dim_NK_cells Classical_Monocytes CM_CD4_T_cells EM_CD4_T_cells Naive_CD4_T_cells Regulatory_CD4_T_cells CM_CD8_T_cells EM_CD8_T_cells Naive_CD8_T_cells TEMRA Memory_B_cells Naive_B_cells Nonclassical_Monocytes MAIT_and_GammaDelta_T_cells
#do
#    bsub -o logs/cond_${i}.o -e logs/cond_${i}.e -R"select[mem>60000] rusage[mem=60000]" -M60000 Rscript /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/12_Onek1k/5_conditional_analysis_forcoloc_onek1k.R $i
#done

for i in All
do
    bsub -o logs/cond_${i}_645.o -e logs/cond_${i}_645.e -q yesterday -R"select[mem>60000] rusage[mem=60000]" -M60000 Rscript /lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Scripts/12_Onek1k/5_conditional_analysis_forcoloc_onek1k.R $i
done


