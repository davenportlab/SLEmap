// This file contains parameters used for running QC on SLEmap with the YASCP pipeline (version 1.7)

params{
    input_data_table="/path/yascp_inputs.tsv"
    input="cellranger"
    genotype_input {
        run_with_genotype_input=true
        vireo_with_gt=true
        posterior_assignment=true
        tsv_donor_panel_vcfs="/path/vcf_inputs.tsv"
        subset_genotypes=true
    }
}