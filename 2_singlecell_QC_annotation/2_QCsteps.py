# ### This is to apply QC steps on the (personal HLA) mapped SLEmap data on 298 donors

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import scvi
import re
import glob

pd.set_option('display.max_columns', None)
sc.set_figure_params(figsize=(5, 5),dpi=200)

SLEmap = sc.read('/path/YASCPresults.h5ad')

# ## Preparing anndata object for totalVI - filtering for hard filter, 3MAD, doublet, VDJdoublet

hardfilter = pd.read_csv("/path/QC_hardfilters.csv",index_col=0)
SLEmap.obs["index"] = SLEmap.obs.index
SLEmap.obs = SLEmap.obs.reset_index().merge(hardfilter, how = 'left',on="index",indicator=False).set_index('index')
newobs = pd.read_csv("/path/QC_allfilters.csv",index_col=0)
SLEmap.obs = SLEmap.obs.reset_index().merge(newobs, how = 'left',on="index",indicator=False).set_index('index')


all_cite_files = glob.glob("/path/yascp/citeseq/DSB/*/*/*.matrix.csv")
CITE = pd.DataFrame()
for f1 in all_cite_files:
    c1 = pd.read_csv(f1,index_col=0)
    CITE=pd.concat([CITE,c1])

SLEmap.obs = SLEmap.obs.merge(CITE, left_on=['barcode_y'],right_index=True,how='left', indicator=True)

# make citeseq data to obsm single cell data : required for totalVI
CITE_2 = SLEmap.obs[CITE.columns].copy()
SLEmap.obsm['protein_expression'] = CITE_2

# keep only cells passing QC
SLEmap = SLEmap[SLEmap.obs["hard_filter_manual"] == "PASShardfilter"]
SLEmap = SLEmap[SLEmap.obs["PASS3MAD"] == "PASS"]
SLEmap = SLEmap[SLEmap.obs["over3doublet"] ==  "singlet"]
SLEmap = SLEmap[SLEmap.obs["vdjfilter"] == "NOTvdjdoublet"]

## Preparing anndata object for totalVI: calculating new highly variable genes with the cells remaining: immune genes removed
SLEmap.obs['pool'] = SLEmap.obs['experiment_id'].str.split('_').str[2]
SLEmap.layers["counts"] = SLEmap.X.copy()

#remove immune receptor genes
v1 = 'IG[HKL][VDJ]|AC233755.*|IGH[GMDEA]|IGKC|IGLC|IGLL|TR[ABGD][CVDJ]'
regex_exclusions = set(SLEmap.var[SLEmap.var.gene_symbols.str.contains(v1)].index)
split_exclusions_gene_symbol = set(SLEmap.var[SLEmap.var['gene_symbols'].isin( v1.replace(',',';').split(';'))].index)
split_exclusions_gene_ensg = set(SLEmap.var[SLEmap.var.index.isin( v1.replace(',',';').split(';'))].index)
all_genes_to_exclude = regex_exclusions.union(split_exclusions_gene_symbol).union(split_exclusions_gene_ensg)
exclude_gene_list = pd.DataFrame(all_genes_to_exclude, columns=['ensembl_gene_id'])
SLEmap = SLEmap[:,list(set(SLEmap.var.index) - set(exclude_gene_list['ensembl_gene_id']))]

sc.pp.normalize_total(
        SLEmap,
        target_sum=1e4,
        exclude_highly_expressed=False,
        key_added='normalization_factor', 
        inplace=True
    )
sc.pp.log1p(SLEmap)

sc.pp.highly_variable_genes(
        SLEmap,
        flavor='seurat',
        n_top_genes=2000,  
        batch_key="pool",
        inplace=True
    )

SLEmap.write("/path/fortotalVI.h5ad")

## Running totalVI
SLEmap = SLEmap[:,SLEmap.var["highly_variable"]]
SLEmap.X = SLEmap.layers['counts'] 

SLEmap = SLEmap.to_memory().copy()
scvi.model.TOTALVI.setup_anndata(SLEmap, protein_expression_obsm_key="protein_expression",batch_key="pool")
model = scvi.model.TOTALVI(SLEmap,latent_distribution="normal",n_layers_decoder=2)
model.train()
model.save("/path/scvi_model",adata=SLEmap, overwrite=True)

# neighbourhood calculation, umap, clustering
SLEmap.obsm["X_totalVI"] = model.get_latent_representation()
sc.pp.neighbors(SLEmap, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap)
sc.tl.leiden(SLEmap, key_added="leiden_totalVI",resolution=2)
SLEmap.write("/path/afterfirsttotalVI.h5ad")

sc.pl.umap(
    SLEmap,
    color=["log1p_n_genes_by_counts","log1p_total_counts", 'log1p_total_counts_gene_group__mito_protein','log1p_total_counts_gene_group__ribo_protein'],size=0.7,ncols=2
)

sc.tl.rank_genes_groups(SLEmap, groupby="leiden_totalVI", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    SLEmap, groupby="leiden_totalVI", standard_scale="var", n_genes=2
)

# # Low Ribo clusters
sc.set_figure_params(figsize=(10, 5),dpi=150)
sc.pl.violin(SLEmap, keys=['log1p_total_counts_gene_group__ribo_protein'], groupby='leiden_totalVI', rotation=90)
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap,color=["leiden_totalVI",'log1p_total_counts_gene_group__ribo_protein'],groups=["25","35"],size=0.8,ncols=2)

# # CITE outlier clusters
protein_expression = SLEmap.obsm["protein_expression"]

protein_SLEmap = sc.AnnData(protein_expression.copy())
protein_SLEmap.raw = protein_SLEmap
protein_SLEmap.obsm["X_umap"] = SLEmap.obsm["X_umap"]
protein_SLEmap.obsm["X_totalVI"] = SLEmap.obsm["X_totalVI"]
protein_SLEmap.obs["leiden_totalVI"] = SLEmap.obs["leiden_totalVI"]

sc.pp.calculate_qc_metrics(
    protein_SLEmap, inplace=True, log1p=True,var_type='Protein',percent_top=None,expr_type='ADTlevels'
)

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(
    protein_SLEmap,
    color=['n_Protein_by_ADTlevels', 'log1p_n_Protein_by_ADTlevels', 'total_ADTlevels', 'log1p_total_ADTlevels'],
    size=0.7,ncols=2
)

sc.pp.log1p(protein_SLEmap)
sc.pl.umap(
    protein_SLEmap,
    color=['anti-human_CD4', 'anti-human_CD8', 'anti-human_CD19', 'anti-human_CD16'],
    size=0.7,ncols=2
)

sc.pl.umap(
    protein_SLEmap,
    color=['Mouse_IgG1_isotype_Ctrl',
 'Mouse_IgG2a_isotype_Ctrl',
 'Mouse_IgG2b_isotype_Ctrl',
 'Rat_IgG2b_Isotype_Ctrl',
 'Rat_IgG1_isotype_Ctrl',
 'Rat_IgG2a_Isotype_Ctrl',
 'Armenian_Hamster_IgG_Isotype_Ctrl'],
    size=0.7,ncols=2
)


sc.tl.rank_genes_groups(protein_SLEmap, groupby="leiden_totalVI", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(protein_SLEmap,groups=["20","24","37","39","27","30"], standard_scale="var", n_genes=5,save="group_dotplot_CITE.pdf")

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap,color=["leiden_totalVI"],groups=["20","24","37","39"],size=0.8,ncols=2)

## Erythrocytes and platelets
sc.pl.umap(
    SLEmap,
    color=["Azimuth:predicted.celltype.l2"],groups=["Eryth","Platelet"],
    size=0.5,ncols=2
)

# remove low ribo clusters
SLEmap = SLEmap[~SLEmap.obs["leiden_totalVI"].isin(["25","35"])]

# remove outlier CITE clusters
SLEmap = SLEmap[~SLEmap.obs["leiden_totalVI"].isin(["20","24","37","39"])]

# remove erthrocytes and platelets
SLEmap = SLEmap[~SLEmap.obs["Azimuth:predicted.celltype.l2"].isin(["Eryth","Platelet"])]

# Running totalVI
SLEmap = SLEmap.raw.to_adata()
SLEmap.raw = SLEmap.copy()
SLEmap = SLEmap[:,SLEmap.var["highly_variable"]] 
SLEmap = SLEmap.to_memory().copy()
scvi.settings.seed = 7
scvi.model.TOTALVI.setup_anndata(SLEmap, protein_expression_obsm_key="protein_expression",batch_key="poolID")
model = scvi.model.TOTALVI(SLEmap,latent_distribution="normal",n_layers_decoder=2)
model.train()
model.save("/path/scvi_model_second",adata=SLEmap, overwrite=True)

SLEmap.obsm["X_totalVI"] = model.get_latent_representation()
sc.pp.neighbors(SLEmap, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap)
sc.pl.umap(SLEmap,color="Azimuth:predicted.celltype.l2") 

sc.tl.leiden(SLEmap, key_added="leiden_totalVI",resolution=2)
sc.pl.umap(SLEmap,color="leiden_totalVI")

SLEmap.write("/path/aftersecondtotalVI.h5ad")
