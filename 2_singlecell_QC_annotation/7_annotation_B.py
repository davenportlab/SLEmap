import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scvi
import copy
import scanpy as sc
import re
import seaborn as sns

# Read in data with B/T/Other annotation and check annotations

SLEmap = sc.read("For_annotation.h5ad")
protein_SLEmap = sc.read("For_annotation_protein.h5ad")

print(SLEmap.obs.externalID.nunique())
#298

sc.pl.umap(
    SLEmap,
    color=["Azimuth:predicted.celltype.l1","Azimuth:predicted.celltype.l2","VDJ","annotation_BTother"],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=2,
)

# Subset to only the 'B' cells

SLEmap_B = SLEmap[SLEmap.obs["annotation_BTother"] == "B"]
protein_SLEmap_B = protein_SLEmap[protein_SLEmap.obs["annotation_BTother"] == "B"]

print(SLEmap.obs.annotation_BTother.value_counts())
# annotation_BTother
# T        507805
# B         79137
# Other     76491

# Add BCR information (minimally QCed)

annot_info = pd.read_csv("bcr_df_heavy_single_bcr.csv")

obs_tmp = SLEmap_B.obs.copy()

obs_tmp["cell_barcode"] = obs_tmp.index.str.split("-SLE").str[0]

obs_tmp = obs_tmp.reset_index().merge(how = 'left', right = annot_info, left_on = ["cell_barcode", "poolID"], right_on = ["barcode", "pool"]).set_index('index')

obs_tmp = obs_tmp.drop(["cell_barcode", "barcode", "pool"], axis = 1)

SLEmap_B.obs = obs_tmp

# Recluster

sc.pp.neighbors(SLEmap_B, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap_B)

sc.tl.leiden(SLEmap_B, key_added="leiden_totalVI",resolution=1)

c_call_props = SLEmap_B.obs.groupby("leiden_totalVI").c_call.value_counts(normalize = True).reset_index()

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.settings.vector_friendly = True
sc.settings.figdir = "Figures/GEX/b_annotation"



sc.pl.umap(
    SLEmap_B,
    color=["leiden_totalVI", "Celltypist:Immune_All_Low:majority_voting", "recruitment_centre",
           "Azimuth:predicted.celltype.l2"],
    legend_fontsize='xx-small',
    ncols=2,
    save = "_b_overall.pdf",
    wspace = 0.2
)

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(
    SLEmap_B,
    color=["VDJ","c_call", "mut_freq_v"],
    legend_fontsize='xx-small',
    ncols=3,
    save = "_bcr_info.pdf",
)

# Annotation using RNA-seq gene markers

# to use gene symbolsl instead of ensembl ids 
SLEmap_B.var_names = SLEmap_B.var["gene_symbols"]

# Key B cell markers

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
sc.pl.dotplot(SLEmap_B, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

marker_genes_dict = {
    'Naive': ['TCL1A', 'MS4A1'],
    'Memory': ["CD27", "CD86" ],
    'ABC': ["ITGAX", "TBX21", "FCRL3", "FCRL5"],
    'ASC': ["CD38"]
}

sc.pl.dotplot(SLEmap_B, marker_genes_dict, 'leiden_totalVI', dendrogram=False,use_raw=False)

B_markers = ["CD27", "CD38", "CD24", "CR2", "FAS", "CD86", "ITGAX", "TBX21", 
               "SLAMF7", "CXCR5",  "SDC1", "IL10", "IL12A", "EBI3", "FCRL3", "FCRL5",
               "IL4R", "MKI67", "TCL1A", "CLEC2B", "CD19", "MS4A1", "ZEB2", "CD83", "CD69", "CCR7"]

sc.pl.dotplot(SLEmap_B, B_markers, 'leiden_totalVI', dendrogram=False,use_raw=False)

sc.set_figure_params(figsize=(5, 5),dpi=150)

sc.pl.umap(
    SLEmap_B,
    color=["c_call", "leiden_totalVI"],
    legend_fontsize='xx-small',
    ncols=2,
)

sc.pl.stacked_violin(SLEmap_B, B_markers, 'leiden_totalVI', dendrogram=False,use_raw=False)

# Annotation using CITE-seq surface protein markers

protein_SLEmap_B.obs["leiden_totalVI"] = SLEmap_B.obs["leiden_totalVI"] 
protein_SLEmap_B.obsm["X_umap"] = SLEmap_B.obsm["X_umap"] 
protein_SLEmap_B.uns["leiden_totalVI_colors"] = SLEmap_B.uns["leiden_totalVI_colors"] 

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(
    protein_SLEmap_B,
    color=['total_ADTlevels', 'log1p_total_ADTlevels','leiden_totalVI'],
    legend_loc='on data',
    legend_fontsize='xx-small',
    ncols=3,
)

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_B, use_raw=False,color=['anti-human_IgM','anti-human_IgD', 'leiden_totalVI'])

sc.set_figure_params(figsize=(5, 5),dpi=150)
sc.pl.umap(protein_SLEmap_B, use_raw=False,color=["anti-human_CD185_(CXCR5)", "anti-human_CD124_(IL-4R_)", 
                                                  "anti-human_Ig_light_chain_K", "anti-human_Ig_light_chain_G", 
                                                  "anti-human_CD38", "anti-human_CD24", "anti-human_CD11c", "anti-human_CD27"])

# Add new annotation to adata

old_to_new = {
'0':'Naive',
'1':'Naive',
'2':'Naive',
'3':'Naive',
'4':'Memory',
'5':'Memory',
'6':'Naive',
'7':'Naive',
'8':'Memory',
'9':'Naive',
'10':'Memory',
'11':'Naive',
'12':'Naive',
'13':'Antibody secreting cell',
'14':'Memory',
'15':'Memory',
'16':'Naive'
}
SLEmap_B.obs['manual_annotation_l1'] = (
SLEmap_B.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

old_to_new = {
'0':'Naive',
'1':'Naive',
'2':'Naive',
'3':'Naive',
'4':'Unswitched memory',
'5':'Switched memory',
'6':'Naive',
'7':'Naive',
'8':'Switched memory',
'9':'Naive',
'10':'Atypical memory cell',
'11':'Naive',
'12':'Naive',
'13':'Antibody secreting cell',
'14':'Switched memory',
'15':'Unswitched memory',
'16':'Naive'

}
SLEmap_B.obs['manual_annotation_l2'] = (
SLEmap_B.obs['leiden_totalVI']
.map(old_to_new)
.astype('category')
)

marker_genes_dict = {
    'Naive': ['TCL1A', 'MS4A1'],
    'Memory': ["CD27", "CD86" ],
    'ABC': ["ITGAX", "TBX21", "FCRL3", "FCRL5"],
    'ASC': ["CD38"]
}

sc.pl.dotplot(SLEmap_B, marker_genes_dict, 'manual_annotation_l2', dendrogram=False,use_raw=False,
             categories_order = ['Naive', 'Unswitched memory', 'Switched memory',
                                 'Atypical memory cell', 'Antibody secreting cell'])

sc.pl.stacked_violin(SLEmap_B, B_markers, 'manual_annotation_l1', dendrogram=False,use_raw=False)

sc.pl.stacked_violin(SLEmap_B, B_markers, 'manual_annotation_l2', dendrogram=False,use_raw=False)

SLEmap_B.write("manualannotation_B.h5ad")




