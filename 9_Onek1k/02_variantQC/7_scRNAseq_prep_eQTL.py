# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: users/hj10/SLEmap_HJ_py/3
#     language: python
#     name: chl0ag9uaehhss9zb2z0cgfjay91c2vycy9oajewl1nmrw1hcf9isl9wes8zcg..
# ---

# %% [markdown]
# # Filtering cell types and genes
# - oneK1K data has been annotated by Haerin

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# ## 1. Set up

# %%
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import re
import ast
from scipy import sparse
from pathlib import Path

# %%
import IPython
IPython.__version__

# %%
#set up directory
os.chdir('/path')

# %%
os.getcwd()

# %% [markdown]
# ## 2. Load and check data

# ### 2.1. annotated oneK1K data

# %%
#load anontated oneK1K data
onek1K = sc.read("/path/Onek1k/Onek1k_female_manual_annotation.h5ad")
onek1K

# %%
#check whether the X is raw/noramlised value
print(onek1K.X.sum(axis = 1))
print(onek1K.X.max())
print(onek1K.X.min())

# %%
plt.figure(figsize=(6, 6))
sc.pl.umap(onek1K, color = ["Celltype_level1"],ncols=1)

# %% [markdown]
# ### 2.2. raw file od oneK1K

# %%
#load raw file od oneK1K
onek1K_raw = sc.read("/path/OneK1k/original_h5ad/078a26dc-0585-4b94-9252-ee2f2dda5742.h5ad")
onek1K_raw

# %%
#check whether the X is raw/noramlised value
onek1K_raw.X.sum(axis = 1)

# %%
#subset only female donors
onek1K_raw = onek1K_raw[onek1K_raw.obs["sex"] == "female"]
onek1K_raw.obs.donor_id.value_counts()

# %%
# remove immune receptor genes
v1 = 'IG[HKL][VDJ]|AC233755.*|IGH[GMDEA]|IGKC|IGLC|IGLL|TR[ABGD][CVDJ]'
gene_exclusions = set(onek1K_raw.var[onek1K_raw.var.feature_name.str.contains(v1)].index)
print(len(set(onek1K_raw.var.index) - gene_exclusions))
onek1K_raw = onek1K_raw[:,list(set(onek1K_raw.var.index) - gene_exclusions)]
onek1K_raw


# %%
if onek1K.obs_names.equals(onek1K_raw.obs_names):
    onek1K_raw.obsm["X_umap"] = onek1K.obsm["X_umap"].copy()
    onek1K_raw.obs['Celltype_level1'] = onek1K.obs['Celltype_level1'].copy()

# %%
plt.figure(figsize=(6, 6))
sc.pl.umap(onek1K_raw, color = ["Celltype_level1"],ncols=1)

# %%
#add annotation for all cells (bulk like) eQTL
onek1K_raw.obs['Celltype_level0'] = 'All'

# %%
#save the object
onek1K_raw.write("onek1K_raw_eQTL_input.h5ad")