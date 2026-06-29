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

SLEmap_T = sc.read("/path/1.SLEmap_T_resolution_2.h5ad")
protein_SLEmap_T = sc.read("/path/1.protein_SLEmap_T_resolution_2.h5ad")

sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD45RA'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_KLRG1_(MAFA)'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD57_Recombinant'], groupby='leiden_totalVI', use_raw=False)

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_T, use_raw=False,color=['NKG7',"GNLY","CCL4"]) 

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(SLEmap_T, ['NKG7'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['GNLY'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['CCL4'], groupby='leiden_totalVI', use_raw=False)

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_T, use_raw=False,color=['CCR7','CD44','CD69','SELL','IL7R'])

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(SLEmap_T, use_raw=False,color=['GZMK','CXCR3','CCL5'])

sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(SLEmap_T, ['GZMK'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['CXCR3'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(SLEmap_T, ['CCL5'], groupby='leiden_totalVI', use_raw=False)

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_T, use_raw=False,color=['anti-human_CD45RA', 'anti-human_CD45RO','anti-human_KLRG1_(MAFA)'])

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_T, use_raw=False,color=['anti-human_CD27'])

sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD27'], groupby='leiden_totalVI', use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CX3CR1'], groupby='leiden_totalVI', use_raw=False)

old_to_new = {
'0':'CD8_naive',
'1':'CD4_naive',
'2':'CD8_TEM',
'3':'CD4_TCM',
'4':'CD4_naive',
'5':'CD4_naive',
'6':'TEMRA',
'7':'CD8_naive',
'8':'CD4_Treg',
'9':'CD4_TCM',
'10':'CD4_naive',
'11':'CD4_TCM',
'12':'Further',
'13':'CD4_TEM',
'14':'Further',
'15':'CD8_TCM',
'16':'TEMRA',
'17':'CD4_CTL',
'18':'DN_T',
'19':'Further',
'20':'20',
'21':'CD4_naive',
'22':'22',
'23':'TEMRA',
'24':'CD8_naive',
'25':'CD4_naive',
'26':'26',
'27':'CD4_TCM',
'28':'CD4_naive',
}
SLEmap_T.obs["leiden_totalVI_firstanot"] = (
SLEmap_T.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)


# %%
sc.tl.leiden(SLEmap_T, key_added="leiden_totalVI_secondanot", restrict_to=("leiden_totalVI_firstanot",["Further"]),resolution=1)
sc.settings.set_figure_params(figsize=(5,5), dpi=150)
sc.pl.umap(
    SLEmap_T,
    color=["leiden_totalVI_firstanot","leiden_totalVI_secondanot"],
)

protein_SLEmap_T.obs["leiden_totalVI_secondanot"] = SLEmap_T.obs["leiden_totalVI_secondanot"]

# gdT and dnT
sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(protein_SLEmap_T, ['anti-human_TCR_'], groupby="leiden_totalVI_secondanot", use_raw=False)

## gdt
sc.pl.violin(protein_SLEmap_T, ['anti-human_TCR_V_2'], groupby="leiden_totalVI_secondanot", use_raw=False)

#mait
sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(SLEmap_T, ['SLC4A10'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.pl.violin(SLEmap_T, ['TRAV1-2'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.pl.violin(protein_SLEmap_T, [ 'anti-human_TCR_V_7.2'], groupby="leiden_totalVI_secondanot", use_raw=False)

# CD4 Treg
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD25'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD27'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD127_(IL-7R_)'], groupby="leiden_totalVI_secondanot", use_raw=False)

# CD4 TCM (CD62L+) TEM(CD62L-)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD62L'], groupby="leiden_totalVI_secondanot", use_raw=False)

#TEM
sc.pl.violin(SLEmap_T, ['GZMA'], groupby="leiden_totalVI_secondanot", use_raw=False)

# %%
old_to_new = {
'20':'20',
'22':'22',
'26':'26',
'CD4_CTL':'CD4_CTL',
'CD4_TCM':'CD4_TCM',
'CD4_TEM':'CD4_TEM',
'CD4_Treg':'CD4_Treg',
'CD4_naive':'CD4_naive',
'CD8_TCM':'CD8_TCM',
'CD8_TEM':'CD8_TEM',
'CD8_naive':'CD8_naive',
'DN_T':'DN_T',
'Further,0':'CD4_TCM',
'Further,1':'Further,1',
'Further,2':'CD4_TCM',
'Further,3':'CD4_TCM',
'Further,4':'GD_T',
'Further,5':'MAIT',
'Further,6':'Further,6',
'Further,7':'MAIT',
'Further,8':'CD4_TCM',
'Further,9':'CD4_TCM',
'Further,10':'CD4_TCM',
'TEMRA':'TEMRA',
}
SLEmap_T.obs['leiden_totalVI_thirdanot'] = (
SLEmap_T.obs['leiden_totalVI_secondanot']
.map(old_to_new)
.astype('category')
)

protein_SLEmap_T.obs['leiden_totalVI_thirdanot'] = (
protein_SLEmap_T.obs['leiden_totalVI_secondanot']
.map(old_to_new)
.astype('category')
)

sc.tl.rank_genes_groups(SLEmap_T, groupby='leiden_totalVI_thirdanot', method="wilcoxon", groups=['20','22','26','Further,1','Further,6'])
sc.tl.rank_genes_groups(protein_SLEmap_T, groupby='leiden_totalVI_thirdanot', method="wilcoxon", groups=['20','22','26','Further,1','Further,6'])
sc.get.rank_genes_groups_df(SLEmap_T, group="20").head(10) #CD4 Naive
sc.get.rank_genes_groups_df(protein_SLEmap_T, group="20").head(10)

sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD45RA'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD45RO'], groupby="leiden_totalVI_secondanot", use_raw=False)

sc.get.rank_genes_groups_df(SLEmap_T, group="22").head(10) # CD8_TEM
sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(SLEmap_T, ['CCL5'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.get.rank_genes_groups_df(protein_SLEmap_T, group="22").head(10)

sc.get.rank_genes_groups_df(SLEmap_T, group="26").head(10)
sc.get.rank_genes_groups_df(protein_SLEmap_T, group="26").head(10)
sc.settings.set_figure_params(figsize=(20,6))
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD161'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD194_(CCR4)'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD196_(CCR6)'], groupby="leiden_totalVI_secondanot", use_raw=False)

sc.pl.violin(SLEmap_T, ['RORC'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.get.rank_genes_groups_df(SLEmap_T, group="Further,1").head(10)
sc.get.rank_genes_groups_df(protein_SLEmap_T, group="Further,1").head(10)
sc.pl.violin(protein_SLEmap_T, ['anti-human_CD56'], groupby="leiden_totalVI_secondanot", use_raw=False)
sc.get.rank_genes_groups_df(SLEmap_T, group="Further,6").head(10)
sc.get.rank_genes_groups_df(protein_SLEmap_T, group="Further,6").head(10)

old_to_new = {
'20':'CD4_naive',
'22':'CD8_TEM',
'26':'CD4_TEM',
'CD4_CTL':'CD4_CTL',
'CD4_TCM':'CD4_TCM',
'CD4_TEM':'CD4_TEM',
'CD4_Treg':'CD4_Treg',
'CD4_naive':'CD4_naive',
'CD8_TCM':'CD8_TCM',
'CD8_TEM':'CD8_TEM',
'CD8_naive':'CD8_naive',
'DN_T':'DN_T',
'GD_T':'GD_T',
'MAIT':'MAIT',
'TEMRA':'TEMRA',
'Further,1':'CD4_TEM',
'Further,6':'CD4_TCM',
}
SLEmap_T.obs["manual_annotation_l1"] = (
SLEmap_T.obs['leiden_totalVI_thirdanot']
.map(old_to_new)
.astype('category')
)

protein_SLEmap_T.obs["manual_annotation_l1"] = (
protein_SLEmap_T.obs['leiden_totalVI_thirdanot']
.map(old_to_new)
.astype('category')
)

old_to_new = {
'20':'CD4_naive',
'22':'CD8_TEM',
'26':'CD4_TEM',
'CD4_CTL':'CD4_CTL',
'CD4_TCM':'CD4_TCM',
'CD4_TEM':'CD4_TEM',
'CD4_Treg':'CD4_Treg',
'CD4_naive':'CD4_naive',
'CD8_TCM':'CD8_TCM',
'CD8_TEM':'CD8_TEM',
'CD8_naive':'CD8_naive',
'DN_T':'DN_T',
'GD_T':'GD_T',
'MAIT':'MAIT',
'TEMRA':'TEMRA',
'Further,1':'CD4_TEM',
'Further,6':'CD4_TCM',
}
SLEmap_T.obs["manual_annotation_l2"] = (
SLEmap_T.obs['leiden_totalVI_thirdanot']
.map(old_to_new)
.astype('category')
)

protein_SLEmap_T.obs["manual_annotation_l2"] = (
protein_SLEmap_T.obs['leiden_totalVI_thirdanot']
.map(old_to_new)
.astype('category')
)

SLEmap_T.write("/path/T_annotation.h5ad")
protein_SLEmap_T.write("/path/protein_T_annotation.h5ad")
