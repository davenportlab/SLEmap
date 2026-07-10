import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scvi
import copy
import scanpy as sc
import re

pd.set_option('display.max_columns', None)
sc.set_figure_params(figsize=(5, 5),dpi=200)

## Read in data with B/T/Other annotation and check annotations
SLEmap = sc.read("/path/For_annotation.h5ad")
protein_SLEmap = sc.read("/path/For_annotation_protein.h5ad")

# ## Subset to only the 'Other' cells and recluster
SLEmap_other = SLEmap[SLEmap.obs["annotation_BTother"] == "Other"]
protein_SLEmap_other = protein_SLEmap[protein_SLEmap.obs["annotation_BTother"] == "Other"]

sc.pp.neighbors(SLEmap_other, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap_other)
sc.tl.leiden(SLEmap_other, key_added="leiden_totalVI",resolution=1)

# ## Annotation using RNA-seq gene markers

# to use gene symbolsl instead of ensembl ids 
SLEmap_other.var_names = SLEmap_other.var["gene_symbols"]

# level1 cell type markers
sc.settings.set_figure_params(figsize=(4,3))
marker_genes_dict = {
    'T-cell': ['CD3D'],
    'B-cell': ['CD79A', 'MS4A1','CD19'],
    'NK': ['GNLY', 'NKG7'],
    'Myeloid':['HLA-DRA','CST3'],
    'Monocytes': ['ITGAM'],
    'Dendritic': ['CD74','HLA-DPA1','FLT3'],
    'ILC':['IL7R','KIT'],
    'HSPC':['CD34']
}
sc.pl.dotplot(SLEmap_other, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

# level2 cell type markers
sc.settings.set_figure_params(figsize=(4,3))
marker_genes_dict = {
    'NK_CD56Bright': ['NCAM1','SELL','GZMK'],
    'NK_proliferating': ['MKI67'],
    'cDC': ['BST2'],
    'pDC':['IL3RA'],
    'Classical_Monocytes': ['CD14'],
    'NonClassical_Monocytes': ['FCGR3A'],
}
sc.pl.dotplot(SLEmap_other, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(SLEmap_other, ['MKI67'], groupby='leiden_totalVI', use_raw=False)

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(SLEmap_other, ['IL7R'], groupby='leiden_totalVI', use_raw=False)

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(SLEmap_other, ['KIT'], groupby='leiden_totalVI', use_raw=False)

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_other, use_raw=False,color=['KIT']) #Proliferating NK

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_other, use_raw=False,color=['NCAM1',"SELL","GZMK"]) #CD56Bright NK

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_other, use_raw=False,color=['HLA-DRA', 'CD74', 'HLA-DPA1', 'FLT3']) #DC

## Annotation using CITE-seq surface protein markers
protein_SLEmap_other.obs["leiden_totalVI"] = SLEmap_other.obs["leiden_totalVI"] 
protein_SLEmap_other.obsm["X_umap"] = SLEmap_other.obsm["X_umap"] 
protein_SLEmap_other.uns["leiden_totalVI_colors"] = SLEmap_other.uns["leiden_totalVI_colors"] 

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_other, use_raw=False,color=['anti-human_CD3','anti-human_CD4','anti-human_CD8','anti-human_CD19','anti-human_CD123','anti-human_CD45RA','anti-human_CD11c','anti-human_CD16','anti-human_CD14','anti-human_CD56','anti-human_CD161','anti-human_CD141_(Thrombomodulin)'])

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_other, use_raw=False,color=["leiden_totalVI",'anti-human_CD56','anti-human_CD62L']) #CD56Bright NK

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_other, use_raw=False,color=['anti-human_CD127_(IL-7R_)', 'anti-human_CD11c', 'anti-human_CD141_(Thrombomodulin)', 'anti-human_CD1c']) #DC 

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD11c'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD16'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD56'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD14'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_HLA-DR'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD123'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD11b'], groupby='leiden_totalVI', use_raw=False)

sc.settings.set_figure_params(figsize=(4,3))
marker_genes_dict = {
    'NK_CD56Bright': ['anti-human_CD56','anti-human_CD62L'],
    'Myeloid':['anti-human_HLA-DR'],
    'Dendritic':['anti-human_CD127_(IL-7R_)', 'anti-human_CD11c', 'anti-human_CD141_(Thrombomodulin)'],
    'cDC':['anti-human_CD1c'],
    'pDC':['anti-human_CD123','anti-human_CD45RA'],
    'Monocytes': ['anti-human_CD11b'],
    'Classical_Monocytes':['anti-human_CD14'],
    'Non_Classical_Monocytes':['anti-human_CD16'],
}
sc.pl.dotplot(protein_SLEmap_other, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)


### First annotation with those to further recluster
old_to_new = {
'0':'CD56Dim_NK',
'1':'CD56Dim_NK',
'2':'Classical_Mono',
'3':'Further',
'4':'CD56Dim_NK',
'5':'CD56Dim_NK',
'6':'Classical_Mono',
'7':'Classical_Mono',
'8':'Classical_Mono',
'9':'CD56Dim_NK',
'10':'Non_classical_Mono',
'11':'Further',
'12':'Further',
'13':'HSPC',
'14':'pDC',
'15':'Classical_Mono',
'16':'Classical_Mono',
'17':'Classical_Mono',
'18':'Classical_Mono',
'19':'Proliferating_NK',
'20':'CD56Dim_NK',
}
SLEmap_other.obs["leiden_totalVI_firstanot"] = (
SLEmap_other.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

sc.tl.leiden(SLEmap_other, key_added="leiden_totalVI_secondanot", restrict_to=("leiden_totalVI_firstanot",["Further"]),resolution=1)
sc.settings.set_figure_params(figsize=(5,5), dpi=150)
sc.pl.umap(
    SLEmap_other,
    color=["leiden_totalVI_secondanot"],
)

protein_SLEmap_other.obs["leiden_totalVI_secondanot"] = SLEmap_other.obs["leiden_totalVI_secondanot"]

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD56'], groupby='leiden_totalVI_secondanot', use_raw=False) #CD56BrightNK
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD62L'], groupby='leiden_totalVI_secondanot', use_raw=False)#CD56BrightNK

sc.pl.violin(SLEmap_other, ['NCAM1'], groupby='leiden_totalVI_secondanot', use_raw=False)#CD56BrightNK
sc.pl.violin(SLEmap_other, ['SELL'], groupby='leiden_totalVI_secondanot', use_raw=False)#CD56BrightNK
sc.pl.violin(SLEmap_other, ['GZMK'], groupby='leiden_totalVI_secondanot', use_raw=False)#CD56BrightNK

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD127_(IL-7R_)'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD11c'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD11b'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD141_(Thrombomodulin)'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC
sc.pl.violin(protein_SLEmap_other, ['anti-human_CD1c'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC

sc.pl.violin(SLEmap_other, ['CD74'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC
sc.pl.violin(SLEmap_other, ['HLA-DPA1'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC
sc.pl.violin(SLEmap_other, ['FLT3'], groupby='leiden_totalVI_secondanot', use_raw=False) #cDC


## Add new annotation to adata
old_to_new = {
'CD56Dim_NK':'CD56Dim_NK',
'Classical_Mono':'Classical_Mono',
'Non_classical_Mono':'Non_classical_Mono',
'HSPC':'HSPC',
'pDC':'pDC',
'Proliferating_NK':'Proliferating_NK',
'Further,0':'CD56Dim_NK',
'Further,1':'CD56Bright_NK',
'Further,2':'CD56Bright_NK',
'Further,3':'CD56Bright_NK',
'Further,4':'cDC',
'Further,5':'CD56Bright_NK',
'Further,6':'CD56Bright_NK',
'Further,7':'ILC',
'Further,8':'CD56Bright_NK',
'Further,9':'cDC',
'Further,10':'CD56Bright_NK',
'Further,11':'CD56Bright_NK',
'Further,12':'cDC'
}
SLEmap_other.obs['manual_annotation_l1'] = (
SLEmap_other.obs["leiden_totalVI_secondanot"]
.map(old_to_new)
.astype('category')
)

protein_SLEmap_other.obs['manual_annotation_l1'] = (
protein_SLEmap_other.obs["leiden_totalVI_secondanot"]
.map(old_to_new)
.astype('category')
)

old_to_new = {
'CD56Dim_NK':'CD56Dim_NK',
'Classical_Mono':'Classical_Mono',
'Non_classical_Mono':'Non_classical_Mono',
'HSPC':'HSPC',
'pDC':'pDC',
'Proliferating_NK':'Proliferating_NK',
'Further,0':'CD56Dim_NK',
'Further,1':'CD56Bright_NK',
'Further,2':'CD56Bright_NK',
'Further,3':'CD56Bright_NK',
'Further,4':'cDC',
'Further,5':'CD56Bright_NK',
'Further,6':'CD56Bright_NK',
'Further,7':'ILC',
'Further,8':'CD56Bright_NK',
'Further,9':'cDC',
'Further,10':'CD56Bright_NK',
'Further,11':'CD56Bright_NK',
'Further,12':'cDC'
}
SLEmap_other.obs['manual_annotation_l2'] = (
SLEmap_other.obs["leiden_totalVI_secondanot"]
.map(old_to_new)
.astype('category')
)

protein_SLEmap_other.obs['manual_annotation_l2'] = (
protein_SLEmap_other.obs["leiden_totalVI_secondanot"]
.map(old_to_new)
.astype('category')
)

# check obs for annotation
SLEmap_other.obs
sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_other, use_raw=False,color=['manual_annotation_l1'])

SLEmap_other.obs['manual_annotation_l1'].value_counts()
SLEmap_other.obs['manual_annotation_l2'].value_counts()

SLEmap_other.write("/path/manualannotation_Other.h5ad")
