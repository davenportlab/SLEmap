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

SLEmap = sc.read("path/For_annotation.h5ad")
protein_SLEmap = sc.read("/path/For_annotation_protein.h5ad")

## Subset to only the 'T' cells and re-cluster
SLEmap_T = SLEmap[SLEmap.obs["annotation_BTother"] == "T" ]
protein_SLEmap_T = protein_SLEmap[protein_SLEmap.obs["annotation_BTother"] == "T"]

sc.pp.neighbors(SLEmap_T, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap_T)
sc.tl.leiden(SLEmap_T, key_added="leiden_totalVI",resolution=2)

sc.pl.umap(
    SLEmap_T,
    color=["leiden_totalVI", "Azimuth:predicted.celltype.l2", "VDJ"],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=2,
)

def cluster_small_multiples(adata, clust_key, size=0.7, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

cluster_small_multiples(SLEmap_T, "Azimuth:predicted.celltype.l2")
cluster_small_multiples(SLEmap_T, "leiden_totalVI")

SLEmap_T.var_names = SLEmap_T.var["gene_symbols"]

sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(SLEmap_T, ['CCR7'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['SELL'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['ANXA1'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['FOXP3'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['CTLA4'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['TIGIT'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['RTKN2'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['PRF1'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['GZMB'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['GZMA'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['GNLY'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['NKG7'], groupby='leiden_totalVI', use_raw=False)

protein_SLEmap_T.obs["leiden_totalVI"] = SLEmap_T.obs["leiden_totalVI"] 
protein_SLEmap_T.obsm["X_umap"] = SLEmap_T.obsm["X_umap"] 
protein_SLEmap_T.uns["leiden_totalVI_colors"] = SLEmap_T.uns["leiden_totalVI_colors"] 

sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD3'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD4'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD45'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD8'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD19'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD45RO'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD45RA'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD62L'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD95_(Fas)'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD25'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD27'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_anti-human_CD127_(IL-7R_)'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD57_Recombinant'], groupby='leiden_totalVI', use_raw=False)


# %%
old_to_new = {
'0':'CD8',
'1':'CD4_naive',
'2':'CD8',
'3':'CD4_TCM',
'4':'CD4_naive',
'5':'CD4_naive',
'6':'CD8',
'7':'CD8',
'8':'CD4_Treg',
'9':'CD4_TCM',
'10':'CD4_naive',
'11':'CD4_TCM',
'12':'Further',
'13':'CD4_TEM',
'14':'Further',
'15':'CD8',
'16':'CD8',
'17':'CD4_CTL',
'18':'DN_T',
'19':'Further',
'20':'20',
'21':'CD4_naive',
'22':'CD8',
'23':'CD8',
'24':'CD8',
'25':'CD4_naive',
'26':'26',
'27':'CD4_TCM',
'28':'CD4_naive',
}

SLEmap_T.obs['annotation_T_l1'] = SLEmap_T.obs['leiden_totalVI'].copy()
SLEmap_T.obs['annotation_T_l1'] = (
SLEmap_T.obs['annotation_T_l1']
.map(old_to_new)
.astype('category')
)

protein_SLEmap_T.obs['annotation_T_l1'] = protein_SLEmap_T.obs['leiden_totalVI'].copy()
protein_SLEmap_T.obs['annotation_T_l1'] = (
protein_SLEmap_T.obs['annotation_T_l1']
.map(old_to_new)
.astype('category')
)

SLEmap_T.write("/path/1.SLEmap_T_resolution_2.h5ad")
protein_SLEmap_T.write("/path/1.protein_SLEmap_T_resolution_2.h5ad")
