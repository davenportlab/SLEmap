import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import scvi
import re
import seaborn

pd.set_option('display.max_columns', None)
sc.set_figure_params(figsize=(5, 5),dpi=200)

SLEmap=sc.read("/path/aftersecondtotalVI.h5ad")

## bring totalVI model back to get denoised cite exp
model = scvi.model.TOTALVI.load("/path/scvi_model_second",adata=SLEmap)

#dont know why but adding two random cite seq markers extracts all the protein data and none of the RNA which is perfect
rna, protein = model.get_normalized_expression(n_samples=25,return_mean=True,gene_list=['anti-human_CD270_(HVEM_TR2)', 'anti-human_CD155_(PVR)'])

SLEmap.obsm["denoised_protein"] = protein

sc.pl.umap(
    SLEmap,
    color=["Azimuth:predicted.celltype.l1","Azimuth:predicted.celltype.l2","VDJ","leiden_totalVI"],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=2,
)

## Get gene expression for all genes and normalizing RNA-seq data for annotation
SLEmap = SLEmap.raw.to_adata()
sc.pp.normalize_total(
        SLEmap,
        target_sum=1e4, #cp10k
        exclude_highly_expressed=False,
        key_added='normalization_factor',  
        inplace=True
    )
sc.pp.log1p(SLEmap)
# you can ignore "WARNING: adata.X seems to be already log-transformed." as this is triggered by SLEmap.uns['log1p'] but X is raw

## Annotating with CITE markers
protein_SLEmap = sc.AnnData(SLEmap.obsm["protein_expression"].copy(), obs=SLEmap.obs)
protein_SLEmap.raw = protein_SLEmap
protein_SLEmap.X = SLEmap.obsm["denoised_protein"]
protein_SLEmap.obsm["X_umap"] = SLEmap.obsm["X_umap"]
sc.pp.calculate_qc_metrics(
    protein_SLEmap, inplace=True, log1p=True,var_type='Protein',percent_top=None,expr_type='ADTlevels'
)

protein_SLEmap.write("/path/For_annotation_protein.h5ad")

sc.pl.umap(
    protein_SLEmap,
    color=['total_ADTlevels', 'log1p_total_ADTlevels'],
    size=2,ncols=2
)
sc.pp.log1p(protein_SLEmap) #https://discourse.scverse.org/t/totalvi-log-normalization-and-non-negativity/421/3
protein_SLEmap.obs["Azimuth:predicted.celltype.l1"].value_counts()

sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(protein_SLEmap, ['anti-human_CD3'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap, ['anti-human_CD4'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap, ['anti-human_CD8'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap, ['anti-human_CD19'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap, ['anti-human_CD14'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap, ['anti-human_CD16'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap, ['anti-human_CD11c'], groupby='leiden_totalVI', use_raw=False)

old_to_new = {
'0':'T',
'1':'B',
'2':'T',
'3':'T',
'4':'T',
'5':'Other',
'6':'T',
'7':'T',
'8':'T',
'9':'T',
'10':'T',
'11':'Other',
'12':'T',
'13':'T',
'14':'T',
'15':'T',
'16':'T',
'17':'B',
'18':'T',
'19':'T',
'20':'T',
'21':'B',
'22':'T',
'23':'T',
'24':'Other',
'25':'B',
'26':'T',
'27':'T',
'28':'T',
'29':'T',
'30':'Other',
'31':'T',
'32':'Other',
'33':'Other',
'34':'B',
'35':'Other',
}
SLEmap.obs["annotation_BTother"] = (
SLEmap.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)
protein_SLEmap.obs["annotation_BTother"] = (
protein_SLEmap.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

SLEmap.write("/path/For_annotation.h5ad")
protein_SLEmap.write("/path/For_annotation_protein.h5ad")
