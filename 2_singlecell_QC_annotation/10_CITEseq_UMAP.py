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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import scvi
import re

pd.set_option('display.max_columns', None)
sc.set_figure_params(figsize=(5, 5),dpi=200)

# %%
protein_SLEmap = sc.read("/path/For_annotation_protein.h5ad")

# %%
sc.pp.log1p(protein_SLEmap) #https://discourse.scverse.org/t/totalvi-log-normalization-and-non-negativity/421/3

# %%
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.set_figure_params(figsize=(2, 2),dpi=120, dpi_save=120)
sc.settings.vector_friendly = True
sc.settings.figdir = "/path/"

sc.pl.umap(
    protein_SLEmap,
        color=['anti-human_CD3', 'anti-human_CD4','anti-human_CD8','anti-human_CD19','anti-human_CD14','anti-human_CD16','anti-human_CD11c','anti-human_CD27'],
    # Setting a smaller point size to get prevent overlap
    size=2,ncols=2,use_raw=False,color_map="magma",save="CITE_marker_forpaper.pdf")

