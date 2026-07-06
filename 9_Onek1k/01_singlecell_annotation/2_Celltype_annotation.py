# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import scanpy as sc
import celltypist
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %%
import matplotlib

# %% [markdown]
# ## look into celltypist and azimuth annotations and leiden clustering

# %%
onek = sc.read("/path/Onek1k/Onek1k_female.h5ad")


# %%
def cluster_small_multiples(adata, clust_key, size=0.7, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

cluster_small_multiples(onek, "predicted.celltype.l2")

# %%
sc.pl.umap(onek, color = ['leiden_new'])

# %%
cluster_small_multiples(onek, "leiden_new")

# %% [markdown]
# ## Annotate T B Other

# %%
import seaborn as sns
ct = pd.crosstab(
    onek.obs["predicted.celltype.l2"],
    onek.obs["leiden_new"]
)
ct_frac = ct.div(ct.sum(axis=1), axis=0)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))

sns.scatterplot(
    data=df,
    x="leiden_new",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

plt.xticks(rotation=90)
plt.xlabel("leiden_new")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.grid()
plt.show()

# %%
import seaborn as sns
ct = pd.crosstab(
    onek.obs["predicted.celltype.l2"],
    onek.obs["leiden_new"]
)
ct_frac = ct.div(ct.sum(axis=0), axis=1)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))

sns.scatterplot(
    data=df,
    x="leiden_new",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

plt.xticks(rotation=90)
plt.xlabel("leiden_new")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.grid()
plt.show()

# %%
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek, ['ENSG00000160654'], groupby='leiden_new', use_raw=False,size=0) #CD3G

# %%
sc.pl.violin(onek, ['ENSG00000167286'], groupby='leiden_new', use_raw=False,size=0) #CD3D

# %%
sc.pl.violin(onek, ['ENSG00000198851'], groupby='leiden_new', use_raw=False,size=0)#CD3E

# %%
old_to_new = {
'0':'T',
'1':'T',
'2':'T',
'3':'T',
'4':'T',
'5':'T',
'6':'Other',
'7':'B',
'8':'T',
'9':'T',
'10':'B',
'11':'T',
'12':'T',
'13':'T',
'14':'Other',
'15':'T',
'16':'Other',
'17':'T',
'18':'T',
'19':'Other',
'20':'Other',
'21':'T',
'22':'T',
'23':'B',
'24':'Platelet',
'25':'Proliferating',
'26':'Other',
'27':'B',
'28':'B',
'29':'Other',
'30':'Other',
'31':'T',
'32':'Erythrocytes',
}
onek.obs["leiden_new_TBOther"] = (
onek.obs['leiden_new']
.map(old_to_new)
.astype('category')
)

# %%
sc.settings.set_figure_params(figsize=(8,8))
sc.pl.umap(onek, color = ['leiden_new_TBOther'],ncols=1)

# %%
onek.write("/path/Onek1k/Annotation/Onek1k_female_TBother.h5ad")

# %% [markdown]
# ## subset to TBOther and recluster - Other

# %%
onek_Other = onek[onek.obs["leiden_new_TBOther"] == "Other"]
sc.pp.neighbors(onek_Other, use_rep="X_harmony",n_neighbors=15)
sc.tl.umap(onek_Other)
sc.tl.leiden(onek_Other, key_added="leiden_new_Other",resolution=1)

# %%
sc.pl.umap(
    onek_Other,
    color=["leiden_new_Other","predicted.celltype.l2"],
    legend_loc='on data',
    legend_fontsize='x-small',
    ncols=2,
)

# %%
import seaborn as sns
ct = pd.crosstab(
    onek_Other.obs["predicted.celltype.l2"],
    onek_Other.obs["leiden_new_Other"]
)
ct_frac = ct.div(ct.sum(axis=0), axis=1)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new_Other",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))

sns.scatterplot(
    data=df,
    x="leiden_new_Other",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

#plt.xticks(rotation=90)
plt.xlabel("leiden_new_Other")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()

# %%
import seaborn as sns
ct = pd.crosstab(
    onek_Other.obs["predicted.celltype.l2"],
    onek_Other.obs["leiden_new_Other"]
)
ct_frac = ct.div(ct.sum(axis=1), axis=0)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new_Other",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))

sns.scatterplot(
    data=df,
    x="leiden_new_Other",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

#plt.xticks(rotation=90)
plt.xlabel("leiden_new_Other")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()

# %%
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_Other, ['ENSG00000149294'], groupby='leiden_new_Other', use_raw=False) #CD56

# %%
old_to_new = {
'0':'CD56Dim_NK_cells',
'1':'CD56Dim_NK_cells',
'2':'CD56Dim_NK_cells',
'3':'Classical_Monocytes',
'4':'CD56Bright_NK_cells',
'5':'Nonclassical_Monocytes',
'6':'Classical_Monocytes',
'7':'cDC',
'8':'pDC',
'9':'HSPC',
'10':'CD56Dim_NK_cells',
'11':'Classical_Monocytes',
'12':'Classical_Monocytes',
}
onek_Other.obs["leiden_new_Other_annotated"] = (
onek_Other.obs['leiden_new_Other']
.map(old_to_new)
.astype('category')
)

# %%
sc.settings.set_figure_params(figsize=(8,8))
sc.pl.umap(onek_Other, color = ["leiden_new_Other_annotated"],ncols=1)

# %%
onek_Other.write("/path/Onek1k/Annotation/Onek1k_female_Othercells.h5ad")

# %%
onek_Other = sc.read("/path/Onek1k/Annotation/Onek1k_female_Othercells.h5ad")
sc.settings.set_figure_params(figsize=(6,6))
sc.pl.umap(onek_Other, color = ["leiden_new_Other_annotated"],ncols=1)

# %% [markdown]
# ## subset to TBOther and recluster - T

# %%
onek_T = onek[onek.obs["leiden_new_TBOther"] == "T"]
sc.pp.neighbors(onek_T, use_rep="X_harmony",n_neighbors=15)


# %%
sc.tl.umap(onek_T)


# %%
sc.tl.leiden(onek_T, key_added="leiden_new_T",resolution=1)

# %%
onek_T = sc.read("/path/Onek1k/Annotation/Onek1k_female_Tcells.h5ad")

# %%
onek_T

# %%
sc.settings.set_figure_params(figsize=(5,5))
sc.pl.umap(
    onek_T,
    color=["leiden_new_T","predicted.celltype.l2",'cell_type'],
    legend_fontsize='x-small',
    ncols=2,
)


# %%
def cluster_small_multiples(adata, clust_key, size=0.7, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

cluster_small_multiples(onek_T, "leiden_new_T")

# %%
import seaborn as sns
ct = pd.crosstab(
    onek_T.obs["predicted.celltype.l2"],
    onek_T.obs["leiden_new_T"]
)
ct_frac = ct.div(ct.sum(axis=0), axis=1)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new_T",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))
plt.grid()
sns.scatterplot(
    data=df,
    x="leiden_new_T",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

plt.xlabel("leiden_new_T")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()

plt.show()

# %%
import seaborn as sns
ct = pd.crosstab(
    onek_T.obs["predicted.celltype.l2"],
    onek_T.obs["leiden_new_T"]
)
ct_frac = ct.div(ct.sum(axis=1), axis=0)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new_T",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))
plt.grid()
sns.scatterplot(
    data=df,
    x="leiden_new_T",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

#plt.xticks(rotation=90)
plt.xlabel("leiden_new_T")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.grid()
plt.show()

# %%
import seaborn as sns
ct = pd.crosstab(
    onek_T.obs["cell_type"],
    onek_T.obs["leiden_new_T"]
)
ct_frac = ct.div(ct.sum(axis=0), axis=1)
df = ct_frac.reset_index().melt(
    id_vars="cell_type",
    var_name="leiden_new_T",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))
plt.grid()
sns.scatterplot(
    data=df,
    x="leiden_new_T",
    y="cell_type",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

#plt.xticks(rotation=90)
plt.xlabel("leiden_new_T")
plt.ylabel("cell_type")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()

plt.show()

# %%
#CD4 CD8
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000010610"],ncols=1) #CD4
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000010610'], groupby='leiden_new_T', use_raw=False,size=0) #CD4

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000153563"],ncols=1) #CD8A
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000153563'], groupby='leiden_new_T', use_raw=False,size=0) #CD8A

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000172116"],ncols=1) #CD8B
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000172116'], groupby='leiden_new_T', use_raw=False,size=0) #CD8B


# %%
#naive CD4 T cells - CCR7, SELL
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000126353"],ncols=1) #CCR7
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000126353'], groupby='leiden_new_T', use_raw=False,size=0) #CCR7
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000188404"],ncols=1) #SELL
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000188404'], groupby='leiden_new_T', use_raw=False,size=0) #SELL
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000213203"],ncols=1) #GIMAP1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000213203'], groupby='leiden_new_T', use_raw=False,size=0) #GIMAP1


# %%
### CD4 TCM - CCR7, SELL, ANXA1
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000135046"],ncols=1) #ANXA1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000135046'], groupby='leiden_new_T', use_raw=False,size=0) #ANXA1


# %%
## Treg - FOXP3, CTLA4, TIGIT, RTKN2
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000049768"],ncols=1) #FOXP3
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000049768'], groupby='leiden_new_T', use_raw=False,size=0) #FOXP3

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000163599"],ncols=1) #CTLA4
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000163599'], groupby='leiden_new_T', use_raw=False,size=0) #CTLA4

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000181847"],ncols=1) #TIGIT
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000181847'], groupby='leiden_new_T', use_raw=False,size=0) #TIGIT

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000182010"],ncols=1) #RTKN2
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000182010'], groupby='leiden_new_T', use_raw=False,size=0) #RTKN2


# %%
## CD4 TEM - CCR7-, SELL-, ITGB1
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000150093"],ncols=1) #ITGB1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000150093'], groupby='leiden_new_T', use_raw=False,size=0) #ITGB1

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000271503"],ncols=1) #CCL5
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000271503'], groupby='leiden_new_T', use_raw=False,size=0) #CCL5



# %%
##CD4 CTL - PRF1, GZMB, GZMA, GZMH, GNLY, NKG7, CX3CR1, FGFBP2,  ITGAL, OASL
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000180644"],ncols=1) #PRF1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000180644'], groupby='leiden_new_T', use_raw=False,size=0) #PRF1

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000100453"],ncols=1) #GZMB
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000100453'], groupby='leiden_new_T', use_raw=False,size=0) #GZMB

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000145649"],ncols=1) #GZMA
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000145649'], groupby='leiden_new_T', use_raw=False,size=0) #GZMA 


# %%
##TEMRA - NKG7, GNLY, CCL4, CCR7-, ITGB1-
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000105374"],ncols=1) #NKG7
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000105374'], groupby='leiden_new_T', use_raw=False,size=0) #NKG7

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000115523"],ncols=1) #GNLY
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000115523'], groupby='leiden_new_T', use_raw=False,size=0) #GNLY

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000275302"],ncols=1) #CCL4
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000275302'], groupby='leiden_new_T', use_raw=False,size=0) #CCL4

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000126353"],ncols=1) #CCR7
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000126353'], groupby='leiden_new_T', use_raw=False,size=0) #CCR7

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000150093"],ncols=1) #ITGB1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000150093'], groupby='leiden_new_T', use_raw=False,size=0) #ITGB1


# %%
##MAIT - SLC4A10,TRAV1-2
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000144290"],ncols=1) #SLC4A10
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000144290'], groupby='leiden_new_T', use_raw=False,size=0) #SLC4A10


# %%
##gd T  (https://academic.oup.com/jleukbio/article/114/6/630/7223408?login=true)
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000109906"],ncols=1) #ZBTB16
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000109906'], groupby='leiden_new_T', use_raw=False,size=0) #ZBTB16

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000204475"],ncols=1) #NCR3
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000204475'], groupby='leiden_new_T', use_raw=False,size=0) #NCR3



# %%
##NK cells - CD3-, CD56 (for CD56Bright NK), GNLY, NKG7

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000160654"],ncols=1) #CD3G
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000160654'], groupby='leiden_new_T', use_raw=False,size=0) #CD3G

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ['ENSG00000167286'],ncols=1) #CD3E
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000167286'], groupby='leiden_new_T', use_raw=False,size=0) #CD3E

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000149294"],ncols=1) #CD56
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000149294'], groupby='leiden_new_T', use_raw=False,size=0) #CD56


sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000105374"],ncols=1) #NKG7
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000105374'], groupby='leiden_new_T', use_raw=False,size=0) #NKG7

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000115523"],ncols=1) #GNLY
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000115523'], groupby='leiden_new_T', use_raw=False,size=0) #GNLY



# %%
##Proliferating
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000148773"],ncols=1) #MKI67
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000148773'], groupby='leiden_new_T', use_raw=False,size=0) #MKI67



# %%
old_to_new = {
'0':'Naive_CD4_T_cells',
'1':'Naive_CD4_T_cells',
'2':'CM_CD4_T_cells',
'3':'EM_CD8_T_cells',
'4':'Further',
'5':'Further',
'6':'Naive_CD8_T_cells',
'7':'Further',
'8':'CM_CD8_T_cells',
'9':'Further',
'10':'Cytotoxic_CD4_T_cells',
'11':'Naive_CD4_T_cells',
'12':'DN_T_cells',

}
onek_T.obs["leiden_new_T_annotated"] = (
onek_T.obs['leiden_new_T']
.map(old_to_new)
.astype('category')
)

# %%

# %%
sc.settings.set_figure_params(figsize=(5,5))
sc.pl.umap(onek_T, color = ["leiden_new_T_annotated"],ncols=1)

# %%
sc.tl.leiden(onek_T, key_added="leiden_new_T_annotated_further", restrict_to=("leiden_new_T_annotated",["Further"]),resolution=0.5)

# %%
sc.settings.set_figure_params(figsize=(5,5))
sc.pl.umap(onek_T, color = ["leiden_new_T_annotated_further"],ncols=1)


# %%
def cluster_small_multiples(adata, clust_key, size=0.7, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

cluster_small_multiples(onek_T, "leiden_new_T_annotated_further")

# %%
import seaborn as sns
ct = pd.crosstab(
    onek_T.obs["predicted.celltype.l2"],
    onek_T.obs["leiden_new_T_annotated_further"]
)
ct_frac = ct.div(ct.sum(axis=0), axis=1)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new_T_annotated_further",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))
plt.grid()
sns.scatterplot(
    data=df,
    x="leiden_new_T_annotated_further",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

plt.xticks(rotation=90)
plt.xlabel("leiden_new_T_annotated_further")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()

plt.show()

# %%
#CD4 CD8
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000010610"],ncols=1) #CD4
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000010610'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CD4

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000153563"],ncols=1) #CD8A
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000153563'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CD8A

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000172116"],ncols=1) #CD8B
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000172116'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CD8B

# %%
#naive CD4 T cells - CCR7, SELL
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000126353"],ncols=1) #CCR7
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000126353'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CCR7
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000188404"],ncols=1) #SELL
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000188404'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #SELL
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000213203"],ncols=1) #GIMAP1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000213203'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #GIMAP1

# %%
### CD4 TCM - CCR7, SELL, ANXA1
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000135046"],ncols=1) #ANXA1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000135046'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #ANXA1

# %%
## CD4 TEM - CCR7-, SELL-, ITGB1
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000150093"],ncols=1) #ITGB1
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000150093'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #ITGB1

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000271503"],ncols=1) #CCL5
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000271503'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CCL5



# %%
## Treg - FOXP3, CTLA4, TIGIT, RTKN2
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000049768"],ncols=1) #FOXP3
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000049768'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #FOXP3

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000163599"],ncols=1) #CTLA4
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000163599'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CTLA4

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000181847"],ncols=1) #TIGIT
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000181847'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #TIGIT

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000182010"],ncols=1) #RTKN2
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000182010'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #RTKN2

# %%
##NK cells - CD3-, CD56 (for CD56Bright NK), GNLY, NKG7

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000160654"],ncols=1) #CD3G
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000160654'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CD3G

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ['ENSG00000167286'],ncols=1) #CD3E
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000167286'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CD3E

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000149294"],ncols=1) #CD56
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000149294'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #CD56


sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000105374"],ncols=1) #NKG7
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000105374'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #NKG7

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000115523"],ncols=1) #GNLY
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000115523'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #GNLY



# %%
##MAIT - SLC4A10,TRAV1-2
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000144290"],ncols=1) #SLC4A10
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000144290'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #SLC4A10

##gd T  (https://academic.oup.com/jleukbio/article/114/6/630/7223408?login=true)
sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000109906"],ncols=1) #ZBTB16
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000109906'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #ZBTB16

sc.settings.set_figure_params(figsize=(4,4))
sc.pl.umap(onek_T, color = ["ENSG00000204475"],ncols=1) #NCR3
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_T, ['ENSG00000204475'], groupby="leiden_new_T_annotated_further", use_raw=False,size=0) #NCR3



# %%
old_to_new = {
'Naive_CD4_T_cells':'Naive_CD4_T_cells',
'CM_CD4_T_cells':'CM_CD4_T_cells',
'EM_CD8_T_cells':'EM_CD8_T_cells',
'Naive_CD8_T_cells':'Naive_CD8_T_cells',
'CM_CD8_T_cells':'CM_CD8_T_cells',
'Cytotoxic_CD4_T_cells':'Cytotoxic_CD4_T_cells',
'DN_T_cells':'DN_T_cells',
'Further,0':'EM_CD4_T_cells', 
'Further,1':'CM_CD4_T_cells',
'Further,2':'TEMRA',
'Further,3':'Naive_CD4_T_cells',
'Further,4':'Regulatory_CD4_T_cells',
'Further,5':'MAIT_and_GammaDelta_T_cells',
'Further,6':'CD56Dim_NK_cells',
}
onek_T.obs["leiden_new_T_annotated_further_annotated"] = (
onek_T.obs["leiden_new_T_annotated_further"]
.map(old_to_new)
.astype('category')
)

# %%
sc.settings.set_figure_params(figsize=(5,5))
sc.pl.umap(onek_T, color = ["leiden_new_T_annotated_further","leiden_new_T_annotated_further_annotated"],ncols=1)

# %%
onek_T.write("/path/Onek1k/Annotation/Onek1k_female_Tcells.h5ad")

# %%

# %% [markdown]
# ## subset to TBOther and recluster - B

# %%
onek = sc.read("/path/Onek1k/Annotation/Onek1k_female_TBother.h5ad")

# %%
onek_B = onek[onek.obs["leiden_new_TBOther"] == "B"]
sc.pp.neighbors(onek_B, use_rep="X_harmony",n_neighbors=15)
sc.tl.umap(onek_B)
sc.tl.leiden(onek_B, key_added="leiden_new_B",resolution=1)

# %%
sc.pl.umap(
    onek_B,
    color=["leiden_new_B","predicted.celltype.l2"],
    legend_loc='on data',
    legend_fontsize='x-small',
    ncols=2,
)

# %%
import seaborn as sns
ct = pd.crosstab(
    onek_B.obs["predicted.celltype.l2"],
    onek_B.obs["leiden_new_B"]
)
ct_frac = ct.div(ct.sum(axis=0), axis=1)
df = ct_frac.reset_index().melt(
    id_vars="predicted.celltype.l2",
    var_name="leiden_new_B",
    value_name="n_cells"
)
plt.figure(figsize=(0.5 * ct_frac.shape[1], 0.5 * ct_frac.shape[0]))

sns.scatterplot(
    data=df,
    x="leiden_new_B",
    y="predicted.celltype.l2",
    size="n_cells",
    hue="n_cells",
    sizes=(10, 500),
    palette="viridis"
)

#plt.xticks(rotation=90)
plt.xlabel("leiden_new_B")
plt.ylabel("predicted.celltype.l2")
plt.title("Cell-type annotation overlap")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()

# %%
## cluster 9 has B cell markers (cd19, cd20) but not clear if memory or naive
## CD27 levels are similar to memory cell levels
sc.settings.set_figure_params(figsize=(15,5))
sc.pl.violin(onek_B, ['ENSG00000139193'], groupby='leiden_new_B', use_raw=False,size=0)


# %%
old_to_new = {
'0':'Naive_B_cells',
'1':'Naive_B_cells',
'2':'Naive_B_cells',
'3':'Memory_B_cells',
'4':'Memory_B_cells',
'5':'Memory_B_cells',
'6':'Memory_B_cells',
'7':'Memory_B_cells',
'8':'Antibody_Secreting_cells',
'9':'Memory_B_cells',
'10':'Naive_B_cells',
}
onek_B.obs["leiden_new_B_annotated"] = (
onek_B.obs['leiden_new_B']
.map(old_to_new)
.astype('category')
)

# %%
sc.settings.set_figure_params(figsize=(6,6))
sc.pl.umap(onek_B, color = ["leiden_new_B_annotated"],ncols=1)

# %%
onek_B.write("/path/Onek1k/Annotation/Onek1k_female_Bcells.h5ad")

# %%

# %% [markdown]
# ## Bring all annotation together

# %%
onek = sc.read("/path/Onek1k/Annotation/Onek1k_female_TBother.h5ad")

# %%
onek.obs["Celltype_level1"] = onek.obs["leiden_new_TBOther"]
onek.obs["Celltype_level1"] = onek.obs["Celltype_level1"].astype(str)

# %%
onek_T = sc.read("/path/Onek1k/Annotation/Onek1k_female_Tcells.h5ad")
onek_T.obs["leiden_new_T_annotated_further_annotated"] = onek_T.obs["leiden_new_T_annotated_further_annotated"].astype(str)
common_indices = onek.obs.index.intersection(onek_T.obs.index)
onek.obs.loc[common_indices, "Celltype_level1"] = onek_T.obs.loc[common_indices, "leiden_new_T_annotated_further_annotated"]
del onek_T

# %%
onek_B = sc.read("/path/Onek1k/Annotation/Onek1k_female_Bcells.h5ad")
onek_B.obs["leiden_new_B_annotated"] = onek_B.obs["leiden_new_B_annotated"].astype(str)
common_indices = onek.obs.index.intersection(onek_B.obs.index)
onek.obs.loc[common_indices, "Celltype_level1"] = onek_B.obs.loc[common_indices, "leiden_new_B_annotated"]
del onek_B

# %%
onek_Other = sc.read("/path/Onek1k/Annotation/Onek1k_female_Othercells.h5ad")
onek_Other.obs["leiden_new_Other_annotated"] = onek_Other.obs["leiden_new_Other_annotated"].astype(str)
common_indices = onek.obs.index.intersection(onek_Other.obs.index)
onek.obs.loc[common_indices, "Celltype_level1"] = onek_Other.obs.loc[common_indices, "leiden_new_Other_annotated"]
del onek_Other

# %%
onek.uns['leiden_new_colors']

# %%
onek.uns['Celltype_level1_colors'] = onek.uns['leiden_new_colors'][1:25]

# %%
sc.settings.set_figure_params(figsize=(6,6))
sc.pl.umap(onek, color = ["Celltype_level1"],ncols=1)

# %%
onek.obs["Celltype_level1"].value_counts()


# %%
def cluster_small_multiples(adata, clust_key, size=0.7, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

cluster_small_multiples(onek, "Celltype_level1")

# %%
onek.write("/path/Onek1k/Onek1k_female_manual_annotation.h5ad")

# %%
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.set_figure_params(figsize=(5,5), dpi=100, dpi_save=300)
sc.settings.vector_friendly = True
sc.settings.figdir = "/path/Onek1k/Annotation"


# %%
onek.obs["Celltype_level1_forplots"] = onek.obs["Celltype_level1"].str.replace('_', ' ')
onek.uns["Celltype_level1_forplots_colors"] = onek.uns["Celltype_level1_colors"]
sc.pl.umap(onek, 
           color=['Celltype_level1_forplots'],
           ncols=1,save="manual_annotation.pdf")

# %%
onek1.var_names = onek1.var["feature_name"]

# %%
onek1.var.index = onek1.var.index.astype(str)
onek1.var_names_make_unique()

# %%
markers = ["CD3G",
"CD3D",
"CD3E",
"CD4",
"CCR7",
"SELL",
"GIMAP1",
"LEF1",
"IL7R",
"ANXA1",
"FOXP3",
"CTLA4",
"TIGIT",
"RTKN2",
"PRF1",
"GZMB",
"ITGB1",
"CCL5",
"CD8A",
"CD8B",
"GNLY",
"NKG7",
"CCL4",
"SLC4A10",
"GZMA",
"TCL1A",
"PAX5",
"BLK",
"MS4A1",
"CD27",
"TBX21",
"JCHAIN",
"CD38",
"HLA-DRA",
"ITGAX",
"CST3",
"ITGAM",
"CD14",
"FCGR3A",
"ZBTB16",
"NCR3",
"NCAM1",
"GZMK",
"CXCR3",
"MKI67",
          ]

# %%
onek1.obs["Celltype_level1"] = onek1.obs["Celltype_level1"].str.replace('_', ' ')

# %%
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.set_figure_params(figsize=(20,10), dpi=200, dpi_save=250)
sc.settings.vector_friendly = True
sc.settings.figdir = "/path/Onek1k/Annotation"


# %%
dp = sc.pl.dotplot(onek1, markers, groupby="Celltype_level1", 
                   standard_scale="var", use_raw=False,
                   categories_order=[
                       "Naive CD4 T cells", "CM CD4 T cells", "EM CD4 T cells", "Regulatory CD4 T cells",'Cytotoxic CD4 T cells',
                       "Naive CD8 T cells", "CM CD8 T cells", "EM CD8 T cells", "TEMRA", "MAIT and GammaDelta T cells",'DN T cells',
                       "Naive B cells", "Memory B cells", "Antibody Secreting cells",
                       "Classical Monocytes", "Nonclassical Monocytes", "CD56Bright NK cells", "CD56Dim NK cells",'pDC','cDC', 'Proliferating','Erythrocytes', 'Platelet','HSPC'],
                   swap_axes=True,
                   show=False,  # don't render yet
                   return_fig=True)  # return the object


# %%
dp.savefig("/path/Onek1k/Annotation/RNAmarkers.pdf")

# %%

# %%

# %%
### extract genes expressed in each cell type
adata = sc.read("/path/Onek1k/onek1K_raw_eQTL_input.h5ad")

# %%
output_dir = "/path/Onek1k/outputs/1_csvfiles/genessequenced/"

cell_types = adata.obs["Celltype_level1"].unique()

for ct in cell_types:
    print(f"Processing {ct}...")
    
    # Subset cells
    adata_sub = adata[adata.obs["Celltype_level1"] == ct]
    
    X = adata_sub.X
    
    # Handle sparse vs dense matrix
    if hasattr(X, "toarray"):  # sparse
        gene_mask = np.array((X > 1).sum(axis=0)).ravel() > 0
    else:  # dense
        gene_mask = (X > 1).sum(axis=0) > 0
    
    genes = adata_sub.var_names[gene_mask]
    
    # Save to CSV
    df = pd.DataFrame({"gene": genes})
    df.to_csv(f"{output_dir}/{ct}_genes.csv", index=False)

# %%
cell_types = adata.obs["Celltype_level0"].unique()

for ct in cell_types:
    print(f"Processing {ct}...")
    
    # Subset cells
    adata_sub = adata[adata.obs["Celltype_level0"] == ct]
    
    X = adata_sub.X
    
    # Handle sparse vs dense matrix
    if hasattr(X, "toarray"):  # sparse
        gene_mask = np.array((X > 1).sum(axis=0)).ravel() > 0
    else:  # dense
        gene_mask = (X > 1).sum(axis=0) > 0
    
    genes = adata_sub.var_names[gene_mask]
    
    # Save to CSV
    df = pd.DataFrame({"gene": genes})
    df.to_csv(f"{output_dir}/{ct}_genes.csv", index=False)
