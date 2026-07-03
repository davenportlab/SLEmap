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

# %% [markdown]
# # Split anndata object by ancestry - AFR (102), EUR (81), EAS (22), SAS (62)

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import scvi
import re
#import plotnine as p9

# %config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
# %config InlineBackend.figure_format='retina'
# %matplotlib inline  

pd.set_option('display.max_columns', None)
sc.set_figure_params(figsize=(5, 5),dpi=200)

# %%
SLEmap=sc.read("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/2_YASCP_final/5_QC_afterHLA_afterClusterremoval/QTLite/testrun_input_n281_noOBSM.h5ad")
SLEmap.obs

# %%
ancestry = pd.read_csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/slemap/phenotype/combined_data_tidied_sc_used_w_KING_06_11_25.csv")
ancestry

# %%
obs_index = SLEmap.obs.index
SLEmap.obs = SLEmap.obs.merge(ancestry[["WGS_ID", "Ancestry"]],on="WGS_ID",how="left")
SLEmap.obs.index = obs_index
SLEmap.obs

# %%
SLEmap.obs["Ancestry"].value_counts()

# %%
SLEmap

# %%
SLEmap_EUR = SLEmap[SLEmap.obs["Ancestry"] == "EUR"].copy()
SLEmap_EUR
SLEmap_EUR.write("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/inputs/singlecell_input_n81_EUR.h5ad")

# %%
SLEmap_AFR = SLEmap[SLEmap.obs["Ancestry"] == "AFR"].copy()
SLEmap_AFR.write("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/inputs/singlecell_input_n102_AFR.h5ad")
SLEmap_AFR

# %%
SLEmap_SAS = SLEmap[SLEmap.obs["Ancestry"] == "SAS"].copy()
SLEmap_SAS.write("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/inputs/singlecell_input_n62_SAS.h5ad")
SLEmap_SAS

# %% [markdown]
# ## make files for vcf filtering

# %%
eur_ids = ancestry.loc[ancestry["Ancestry"] == "EUR", "WGS_ID"].dropna().astype(str).drop_duplicates()
n_eur = len(eur_ids)
print(f"Number of EUR samples: {n_eur}")
eur_ids.to_csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/EUR/inputs/EUR_samples.txt", index=False, header=False)


# %%
afr_ids = ancestry.loc[ancestry["Ancestry"] == "AFR", "WGS_ID"].dropna().astype(str).drop_duplicates()
n_afr = len(afr_ids)
print(f"Number of AFR samples: {n_afr}")
afr_ids.to_csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/AFR/inputs/AFR_samples.txt", index=False, header=False)


# %%
sas_ids = ancestry.loc[ancestry["Ancestry"] == "SAS", "WGS_ID"].dropna().astype(str).drop_duplicates()
n_sas = len(sas_ids)
print(f"Number of SAS samples: {n_sas}")
sas_ids.to_csv("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/20_ancestry_coloc/SAS/inputs/SAS_samples.txt", index=False, header=False)


# %%
