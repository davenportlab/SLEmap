# # making input file for HGI eQTL pipeline
# - get all genes back into the anndata object
# - make sure X is raw count
# - make sure there's a column for sampleID and celltype
# - remove healthy controls and individuals without WGS (281 donors)
# - make donor_id match vcf sample names (remove pool from convoluted_samplename)

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

# Get manual annotation
SLEmap=sc.read("/path/Annotation_done.h5ad")

# get raw anndata object
SLEmap_freeze = sc.read("/path/aftersecondtotalVI.h5ad")
SLEmap_freeze = SLEmap_freeze.raw.to_adata()
SLEmap.raw = SLEmap_freeze.copy()
SLEmap.write("/path/Annotation_done.h5ad")

SLEmap = SLEmap.raw.to_adata()
SLEmap.obs["experiment_id"].value_counts()

# ## remove donors without WGS
SLEmap = SLEmap[~SLEmap.obs['experiment_id'].isin(['SLE_donor4_LDP08',
                                                  'SLE_donor3_LDP10',
                                                  'SLE_donor4_LDP12',
                                                  'SLE_donor4_LDP15',
                                                  'SLE_donor3_LDP17',
                                                  'SLE_donor3_LDP30',
                                                  'SLE_donor3_LDP31',
                                                  'SLE_donor4_LDP43',
                                                  'SLE_donor4_LDP53',
                                                  'SLE_donor3_LDP65',
                                                  'SLE_donor3_LDP67',
                                                  'SLE_donor4_LDP69'])]


## remove healthy
SLEmap = SLEmap[~SLEmap.obs['experiment_id'].isin(['SLE_map13436235_LDP29',
                                                           'SLE_WGS13446469_LDP04',
                                                           'SLE_WGS13446470_LDP06',
                                                           'SLE_WGS13446471_LDP09',
                                                           'SLE_WGS13446474_LDP21'])]

# get donor_id match vcf file
SLEmap.obs["donor_id"] =SLEmap.obs["experiment_id"].str.split('_').str[:-1].str.join('_')

## add column for all cells - all cell pseudobulk
SLEmap.obs["Celltype_level0"] = "All"
del SLEmap.obsm
SLEmap.write("/path/testrun_input_n281_noOBSM.h5ad")