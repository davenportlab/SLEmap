import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scvi
import copy
import scanpy as sc
import re
import matplotlib


SLEmap_other = sc.read("/path/manualannotation_Other.h5ad")
SLEmap_T = sc.read("/path/T_annotation.h5ad")
SLEmap_B = sc.read("/path/manualannotation_B.h5ad")

SLEmap = sc.read("/path/For_annotation.h5ad")

SLEmap.obs["Celltype_level1"] = ""
SLEmap.obs["Celltype_level1"] = SLEmap.obs["Celltype_level1"].astype(str)
SLEmap_other.obs['manual_annotation_l1'] = SLEmap_other.obs['manual_annotation_l1'].astype(str)
common_indices = SLEmap.obs.index.intersection(SLEmap_other.obs.index)
SLEmap.obs.loc[common_indices, "Celltype_level1"] = SLEmap_other.obs.loc[common_indices, 'manual_annotation_l1']
SLEmap_T.obs['manual_annotation_l1'] = SLEmap_T.obs['manual_annotation_l1'].astype(str)
common_indices = SLEmap.obs.index.intersection(SLEmap_T.obs.index)
SLEmap.obs.loc[common_indices, "Celltype_level1"] = SLEmap_T.obs.loc[common_indices, 'manual_annotation_l1']
SLEmap_B.obs['manual_annotation_l1'] = SLEmap_B.obs['manual_annotation_l1'].astype(str)
common_indices = SLEmap.obs.index.intersection(SLEmap_B.obs.index)
SLEmap.obs.loc[common_indices, "Celltype_level1"] = SLEmap_B.obs.loc[common_indices, 'manual_annotation_l1']

SLEmap.obs["Celltype_level2"] = ""
SLEmap.obs["Celltype_level2"] = SLEmap.obs["Celltype_level2"].astype(str)
SLEmap_other.obs['manual_annotation_l2'] = SLEmap_other.obs['manual_annotation_l2'].astype(str)
common_indices = SLEmap.obs.index.intersection(SLEmap_other.obs.index)
SLEmap.obs.loc[common_indices, "Celltype_level2"] = SLEmap_other.obs.loc[common_indices, 'manual_annotation_l2']
SLEmap_T.obs['manual_annotation_l2'] = SLEmap_T.obs['manual_annotation_l2'].astype(str)
common_indices = SLEmap.obs.index.intersection(SLEmap_T.obs.index)
SLEmap.obs.loc[common_indices, "Celltype_level2"] = SLEmap_T.obs.loc[common_indices, 'manual_annotation_l2']
SLEmap_B.obs['manual_annotation_l2'] = SLEmap_B.obs['manual_annotation_l2'].astype(str)
common_indices = SLEmap.obs.index.intersection(SLEmap_B.obs.index)
SLEmap.obs.loc[common_indices, "Celltype_level2"] = SLEmap_B.obs.loc[common_indices, 'manual_annotation_l2']

SLEmap.obs["Celltype_level1"].value_counts()
SLEmap.obs["Celltype_level2"].value_counts()

old_to_new = {
    'CD4_naive':'Naive_CD4_T_cells',
    'CD4_TCM':'CM_CD4_T_cells',
    'CD8_naive':'Naive_CD8_T_cells',
    'B_naive':'Naive_B_cells',
    'TEMRA':'TEMRA',
    'CD8_TEM':'EM_CD8_T_cells',
    'CD56Dim_NK':'CD56Dim_NK_cells',
    'CD4_Treg':'Regulatory_CD4_T_cells',
    'Classical_Mono':'Classical_Monocytes',
    'CD4_TEM':'EM_CD4_T_cells',
    'CD8_TCM':'CM_CD8_T_cells',
    'CD4_CTL':'Cytotoxic_CD4_T_cells',
    'B_switched_memory':'Memory_B_cells',
    'MAIT':'MAIT_cells',
    'GD_T':'GammaDelta_T_cells',
    'CD56Bright_NK':'CD56Bright_NK_cells',
    'DN_T':'DN_T_cells',
    'B_unswitched_memory':'Memory_B_cells',
    'B_abc':'Memory_B_cells',
    'Non_classical_Mono':'Nonclassical_Monocytes',
    'cDC':'cDC',
    'HSPC':'HSPC',
    'B_asc':'Antibody_Secreting_cells',
    'pDC':'pDC',
    'ILC':'ILC',
    'Proliferating_NK':'Proliferating_NK_cells',
}
SLEmap.obs['Celltype_level1'] = (
SLEmap.obs['Celltype_level2']
.map(old_to_new)
.astype('category')
)

#unifying names
old_to_new = {
    'CD4_naive':'Naive_CD4_T_cells',
    'CD4_TCM':'CM_CD4_T_cells',
    'CD8_naive':'Naive_CD8_T_cells',
    'B_naive':'Naive_B_cells',
    'TEMRA':'TEMRA',
    'CD8_TEM':'EM_CD8_T_cells',
    'CD56Dim_NK':'CD56Dim_NK_cells',
    'CD4_Treg':'Regulatory_CD4_T_cells',
    'Classical_Mono':'Classical_Monocytes',
    'CD4_TEM':'EM_CD4_T_cells',
    'CD8_TCM':'CM_CD8_T_cells',
    'CD4_CTL':'Cytotoxic_CD4_T_cells',
    'B_switched_memory':'Switched_Memory_B_cells',
    'MAIT':'MAIT_cells',
    'GD_T':'GammaDelta_T_cells',
    'CD56Bright_NK':'CD56Bright_NK_cells',
    'DN_T':'DN_T_cells',
    'B_unswitched_memory':'Unswitched_Memory_B_cells',
    'B_abc':'Atypical_Memory_B_cells',
    'Non_classical_Mono':'Nonclassical_Monocytes',
    'cDC':'cDC',
    'HSPC':'HSPC',
    'B_asc':'Antibody_Secreting_cells',
    'pDC':'pDC',
    'ILC':'ILC',
    'Proliferating_NK':'Proliferating_NK_cells',
}
SLEmap.obs['Celltype_level2'] = (
SLEmap.obs['Celltype_level2']
.map(old_to_new)
.astype('category')
)

SLEmap.obs["Celltype_level1"].value_counts()
SLEmap.obs["Celltype_level2"].value_counts()

SLEmap = sc.read("/path/Annotation_done.h5ad")

SLEmap.obs["Celltype_level1_forplots"] = SLEmap.obs["Celltype_level1"].str.replace('_', ' ')
SLEmap.obs["Celltype_level2_forplots"] = SLEmap.obs["Celltype_level2"].str.replace('_', ' ')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.set_figure_params(figsize=(4,4.5), dpi=100, dpi_save=200)
sc.settings.vector_friendly = True
sc.settings.figdir = "/path"

sc.pl.umap(SLEmap, 
           color=['Celltype_level1_forplots', 'Celltype_level2_forplots'],
           ncols=1,save="manual_annotation.pdf")

# to use gene symbolsl instead of ensembl ids 
SLEmap.var_names = SLEmap.var["gene_symbols"]

markers = ['LEF1','CCR7', 'SELL', 'IL7R','ANXA1','FOXP3', 'CTLA4', 'TIGIT', 'RTKN2', 'PRF1', 'GZMB', 
            'GNLY', 'NKG7','CCL4','SLC4A10','TRAV1-2', #Tcells
           'TCL1A','PAX5','BLK','MS4A1', 'CD27','TBX21','JCHAIN','CD38', #B cells
           'HLA-DRA','ITGAX','CST3','ITGAM', 'CD14','FCGR3A', #Monocytes
           'CD74','HLA-DPA1','GPR183','FLT3','IL3RA', # DC
           'NCAM1','GZMK','CXCR3','CCL5','MKI67', #NK
           'CD34','KIT', #HSPC,ILC
          ]

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.set_figure_params(figsize=(20,10), dpi=200, dpi_save=250)
sc.settings.vector_friendly = True
sc.settings.figdir = "/path"
sc.pl.dotplot(SLEmap, markers, groupby="Celltype_level2_forplots", standard_scale="var",use_raw=False,categories_order=[
    "Naive CD4 T cells", "CM CD4 T cells", "EM CD4 T cells","Cytotoxic CD4 T cells",  "Regulatory CD4 T cells", 
    "Naive CD8 T cells", "CM CD8 T cells","EM CD8 T cells","TEMRA", "DN T cells","MAIT cells", "GammaDelta T cells", 
    "Naive B cells","Unswitched Memory B cells", "Switched Memory B cells","Atypical Memory B cells","Antibody Secreting cells",
    "Classical Monocytes","Nonclassical Monocytes","CD56Bright NK cells", "CD56Dim NK cells","Proliferating NK cells", 
    "cDC", "pDC", "HSPC", "ILC"],swap_axes=True,save="RNAmarkers.pdf")


protein_SLEmap = sc.AnnData(SLEmap.obsm["protein_expression"].copy(), obs=SLEmap.obs)
protein_SLEmap.raw = protein_SLEmap
protein_SLEmap.X = SLEmap.obsm["denoised_protein"]
protein_SLEmap.obsm["X_umap"] = SLEmap.obsm["X_umap"]
sc.pp.log1p(protein_SLEmap)

markers = ['anti-human_CD3', 'anti-human_CD45', 'anti-human_TCR_','anti-human_CD4', 'anti-human_CD45RA', 'anti-human_CD62L', 'anti-human_CD45RO',
            'anti-human_CD95_(Fas)', 'anti-human_CD25', 'anti-human_CD27', 'anti-human_CD127_(IL-7R_)', 'anti-human_CD57_Recombinant',
            'anti-human_CD8', 'anti-human_KLRG1_(MAFA)', 'anti-human_TCR_V_2', 'anti-human_TCR_V_7.2', #T cells
           'anti-human_CD40', 'anti-human_IgD', 'anti-human_IgM', 'anti-human_CD38', # B cells
           'anti-human_HLA-DR', 'anti-human_CD11b',  'anti-human_CD14', 'anti-human_CD16', #mono
           'anti-human_CD56', #NK
           'anti-mouse_human_CD44','anti-human_CD11c','anti-human_CD141_(Thrombomodulin)','anti-human_CD123', 'anti-human_CD1c', #DC
          ]

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.set_figure_params(figsize=(15,10), dpi=200, dpi_save=250)
sc.settings.vector_friendly = True
sc.settings.figdir = "/path"
sc.pl.dotplot(protein_SLEmap, markers, groupby="Celltype_level2_forplots", standard_scale="var",use_raw=False,cmap="Blues",categories_order=[
    "Naive CD4 T cells", "CM CD4 T cells", "EM CD4 T cells","Cytotoxic CD4 T cells",  "Regulatory CD4 T cells", 
    "Naive CD8 T cells", "CM CD8 T cells","EM CD8 T cells","TEMRA", "DN T cells","MAIT cells", "GammaDelta T cells", 
    "Naive B cells","Unswitched Memory B cells", "Switched Memory B cells","Atypical Memory B cells","Antibody Secreting cells",
    "Classical Monocytes","Nonclassical Monocytes","CD56Bright NK cells", "CD56Dim NK cells","Proliferating NK cells", 
    "cDC", "pDC", "HSPC", "ILC"],save="CITEmarkers.pdf")

##CITEseq UMAP
protein_SLEmap = sc.read("/lustre/scratch127/open-targets/Projects/OTAR2064/working/users/hj10/Results/2_YASCP_final/6_annotation/For_annotation_protein.h5ad")
sc.pp.log1p(protein_SLEmap) 

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sc.set_figure_params(figsize=(2, 2),dpi=120, dpi_save=120)
sc.settings.vector_friendly = True
sc.settings.figdir = "/path"

sc.pl.umap(
    protein_SLEmap,
        color=['anti-human_CD3', 'anti-human_CD4','anti-human_CD8','anti-human_CD19','anti-human_CD14','anti-human_CD16','anti-human_CD11c','anti-human_CD27'],
    size=2,ncols=2,use_raw=False,color_map="magma",save="CITE_marker_forpaper.pdf")
