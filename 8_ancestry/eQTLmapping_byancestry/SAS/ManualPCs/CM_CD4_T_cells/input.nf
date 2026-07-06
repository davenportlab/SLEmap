params{
    method='single_cell'
    input_vcf='/path/ancestry_coloc/SAS/inputs/281_samples.allChrs.hwe_0.000001_maf_0.05_miss_0.95_updateID_SAS_n62_ancestrymaf_0.05.vcf.gz'
    annotation_file = '/path/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    phenotype_file = '/path/ancestry_coloc/SAS/inputs/singlecell_input_n62_SAS.h5ad'
    aggregation_columns='Celltype_level1'
    gt_id_column='donor_id'
    sample_column='experiment_id'
    sample_covariates='' // does not do anything. there's a extra_covariates_file parameter below
    n_min_individ = '0'
    n_min_cells = '20' //more or equal to 
    hwe= 0.000001
    windowSize=1000000
    aggregation_method = 'dMean' // can be: dMean, dSum or both separated by comma
    inverse_normal_transform = 'FALSE' // Inverse normal trasnform data as part of normalisation
    numberOfPermutations=10000
    filter_method = 'None' // filterByExpr|HVG|None
    norm_method = 'NONE'  //'DESEQ|TMM|NONE'
    split_aggregation_adata=false
    percent_of_population_expressed=0.05 // whats the proportion of individuals that has to have the value !=0, Please do not go below 2% as this will result in a lot of low expressed genes which will cause issues in SAIGE and TensorQTL
    chunkSize=10
    cell_percentage_threshold=0.01
    outdir='results'
    maf=0.05
    dMean_norm_method = 'cp10k'
    plink2_filters = "--snps-only --update-sex /path/update_sex.txt --rm-dup exclude-all --merge-par --lax-chrx-import" //merge-par and lax-chrx-import added because the pseudoautosomal region can be treated like any other region of the X chr since all individuals are female. In the plink2 version we use (a.5.12), only using merge-par would produce an error 
    bcftools_filters = '--max-alleles 2 -m2 -M2 -v snps' // removed --known to run on all variants (with updatedIDs in vcf file, it should run the same with both with/without --known)
    copy_mode = 'copy'

    TensorQTL{
        run=true // Are we running TensorQTL?
        aggregation_subentry = 'CM_CD4_T_cells'
        chromosomes_to_test=''
        optimise_pcs = true // Whether to pick the most optimal PCs to use for the downstram analysis. 
        interaction_file='' // to run interaction, provide a TSV with genotype ID and interaction 
        interaction_maf = 0.1 // what interaction MAF to use in tensorqtl interactions test
        interaction_pc_cor_threshold = 0.25 // Drop PCs correlated with interaction above this threshold. Use 1 if you don not want to drop PCs
        interaction_gsea = false // Run GSEA on the interaction terms
        trans_by_cis=true // Run trans-by-cis analysis (all genes, limiting variants) following OPTIM PCs?
        trans_by_cis_variant_list='' // Provide a tsv with variant_id and condition_name OR leave empty to use all lead variants from eGenes
        trans_of_cis=true // Run trans-of-cis eQTL analysis (all variants, limiting genes), following OPTIM PCs?
        trans_by_cis_pval_threshold = 0.01 // when running trans by cis what pval threshold to use.
        alpha=0.05  ///What Alpha value to use in tensorqtl
        chrom_to_map_trans = '' // Option to run GWAS trans analysis - i.e if you specify 2 in here Tensorqtl will run all genes acros genome against all the SNPs on the chromosome 2.
        map_independent_qtls = true
    }
    
    genotypes{
        // subset_genotypes_to_available=false // if true, then the expression data will be first processed and then the samples availbale in all the expression data will be subset from genotype files
        // apply_bcftools_filters=true // if true and vcf file is provided then preprocessing will be done on the files.
        use_gt_dosage = true // whether to use dosage. This will convert the vcf/bed to pgen
    }
    use_gt_dosage = true

    covariates{
        nr_phenotype_pcs = '10' // this is used for Tensorqtl
        nr_genotype_pcs = 2
        extra_covariates_file = ''
    }

}

process{
    withName: NORMALISE_ANNDATA{
        memory = { 50.GB * task.attempt }
    }
    withName: H5AD_TO_SAIGE_FORMAT{
        memory = { 100.GB * task.attempt }
    }
    withName: PHENOTYPE_PCs{
        memory = { 500.GB * task.attempt }
        queue = 'yesterday'
    }
    withName: PREPERE_EXP_BED{
        cpus  = 1 
        memory = { 50.GB * task.attempt }
    }
    withName: GATHER_DATA{
        
        time   = 12.h
        queue = 'long'
        maxRetries = 1
    }
    withName: DETERMINE_TSS_AND_TEST_REGIONS{
        memory = { 8.GB * task.attempt }

    }
}
