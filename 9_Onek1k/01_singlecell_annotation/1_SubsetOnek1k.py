# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
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

# %%
onek = sc.read("/path/OneK1k/original_h5ad/078a26dc-0585-4b94-9252-ee2f2dda5742.h5ad")
onek # 1248980 × 36469

# %%
onek = onek[onek.obs["sex"] == "female"]
onek.obs.donor_id.value_counts()

# %%
## remove immune receptor genes
v1 = 'IG[HKL][VDJ]|AC233755.*|IGH[GMDEA]|IGKC|IGLC|IGLL|TR[ABGD][CVDJ]'
regex_exclusions = set(onek.var[onek.var.feature_name.str.contains(v1)].index)
exclude_gene_list = pd.DataFrame(regex_exclusions, columns=['ensembl_gene_id'])
onek = onek[:,list(set(onek.var.index) - set(exclude_gene_list['ensembl_gene_id']))]
onek

# %%
sc.pp.normalize_total(onek, target_sum = 1e4)
sc.pp.log1p(onek)
#onek.X.expm1().sum(axis = 1)

# %%
## making new umap and neighborhood graph
sc.pp.neighbors(onek, use_rep="X_harmony",n_neighbors=15)
sc.tl.umap(onek)

# %%
sc.tl.leiden(onek, key_added="leiden_new",resolution=2)

# %% [markdown]
# ## Explore Onek1k dataset

# %%
statistics.median(onek.obs["donor_id"].value_counts())

# %%
min(onek.obs["donor_id"].value_counts())

# %%
max(onek.obs["donor_id"].value_counts())

# %%
sc.pl.umap(onek, color = ['cell_type','predicted.celltype.l2'],ncols=1)

# %%
import matplotlib.pyplot as plt
plt.figure(figsize=(15, 6), dpi=80)
# Plotting a basic histogram
plt.hist(onek.obs["nCount_RNA"], bins=1000, color='skyblue', edgecolor='black')

# Adding labels and title
plt.xlabel('nCount_RNA')
plt.ylabel('Frequency')
plt.title('Onek1k: Read count per cell')
ax = plt.gca()
ax.set_xlim([0, 20000])
# Display the plot
plt.show()

# %%
import statistics
statistics.median(onek.obs["nCount_RNA"])


# %%
import matplotlib.pyplot as plt
plt.figure(figsize=(15, 6), dpi=80)
# Plotting a basic histogram
plt.hist(onek.obs["nFeature_RNA"], bins=200, color='skyblue', edgecolor='black')

# Adding labels and title
plt.xlabel('nFeature_RNA')
plt.ylabel('Frequency')
plt.title('Onek1k: Gene count per cell')
ax = plt.gca()
ax.set_xlim([0, 5000])
# Display the plot
plt.show()

# %%
statistics.median(onek.obs["nFeature_RNA"])

# %%
statistics.median(onek.obs.age)

# %%
min(onek.obs.age)

# %%
max(onek.obs.age)

# %%
import matplotlib.pyplot as plt
plt.figure(figsize=(15, 6), dpi=80)
# Plotting a basic histogram
plt.hist(onek.obs["age"], bins=50, color='skyblue', edgecolor='black')

# Adding labels and title
plt.xlabel('age')
plt.ylabel('Frequency')
plt.title('Onek1k: age of females')
ax = plt.gca()
ax.set_xlim([0, 100])
# Display the plot
plt.show()

# %%
onek.obs.self_reported_ethnicity.value_counts()

# %%
onek.obs.disease.value_counts()

# %%
onek.write("/path/Onek1k/Onek1k_female.h5ad")

