# Import packages and setup

%matplotlib inline
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 400)

# Functions

def get_mut_freq_v(sequence, germline, v_germline_end):
    mutations = 0
    length = 0
    for i, base in enumerate(sequence[:v_germline_end]):
        if base != "." and base != "-" and base != "N" and germline[i] != "N" and germline[i] != "-" and germline[i] !=  "." :
            length += 1
            if base != germline[i]:
                mutations += 1
    if length == 0:
        return np.nan
    else:
        return mutations/length

# Read in IgBLAST BCR data

source_folder = "IgBLAST_output_BCR"

file_list = [x for x in os.listdir(f'{source_folder}') if x.endswith('pass.tsv')]

file_dfs = []

for file in file_list:    
    file_df = pd.read_csv(f'{source_folder}/{file}', sep='\t')
    file_df['pool'] =  file.split("_")[0]
    file_dfs.append(file_df)
    
bcr_df = pd.concat(file_dfs)

print(bcr_df.pool.nunique())
#68

# Read in sample, pool and batch info

batch_info = pd.read_csv("slemap_cellranger_to_pool.txt", sep = "\t", header = None)
batch_info.columns = ["sample", "pool", "irods_loc"]

bcr_df = pd.merge(bcr_df, batch_info)

bcr_df["contig"] = bcr_df["sequence_id"].str.split("_", n = 1).str[1]
bcr_df["chain"] = np.where(bcr_df['locus'].str.contains('IGH'), "Heavy" , "Light")
bcr_df["sample_cell_id"] = bcr_df["sample"] + "_" + bcr_df["cell_id"]

chain_counts_bcr = bcr_df.groupby(['pool','chain'])['sample_cell_id'].agg(number_chains = "value_counts").reset_index()
chain_counts_heavy = chain_counts_bcr.query("chain == 'Heavy'").copy()

chain_counts_heavy["BCR_doublet"] = np.where(chain_counts_heavy['number_chains'] == 1, "No" , "Yes")

chain_counts_heavy[["sample_cell_id", "pool", "BCR_doublet"]]

chain_counts_heavy["barcode"] = chain_counts_heavy["sample_cell_id"].str.rsplit("_", n=1).str[1]

chain_counts_heavy.to_csv("bcr_doublet_barcodes.csv", index = False)

has_single_bcr = chain_counts_heavy.query("number_chains == 1").copy()
has_single_bcr["VDJ"] = "BCR"

# Read in IgBLAST TCR data

source_folder = "IgBLAST_output_TCR"

file_list = [x for x in os.listdir(f'{source_folder}') if x.endswith('pass.tsv')]

file_dfs = []

for file in file_list:    
    file_df = pd.read_csv(f'{source_folder}/{file}', sep='\t')
    file_df['pool'] =  file.split("_")[0]
    file_dfs.append(file_df)
    
tcr_df = pd.concat(file_dfs)

print(tcr_df.pool.nunique())
#69

# Read in sample, pool and batch info

tcr_df = pd.merge(tcr_df, batch_info)

tcr_df["contig"] = tcr_df["sequence_id"].str.split("_", n = 1).str[1]
tcr_df["chain"] = np.where(tcr_df['locus'].str.contains('TRB'), "Beta" , "Alpha")
tcr_df["sample_cell_id"] = tcr_df["sample"] + "_" + tcr_df["cell_id"]

chain_counts_tcr = tcr_df.groupby(['pool', 'chain'])['sample_cell_id'].agg(number_chains = "value_counts").reset_index()
chain_counts_beta = chain_counts_tcr.query("chain == 'Beta'").copy()

chain_counts_beta["tcr_doublet"] = np.where(chain_counts_beta['number_chains'] == 1, "No" , "Yes")

chain_counts_beta[["sample_cell_id", "pool", "tcr_doublet"]]

chain_counts_beta["barcode"] = chain_counts_beta["sample_cell_id"].str.rsplit("_", n=1).str[1]

chain_counts_beta.to_csv("tcr_doublet_barcodes.csv", index = False)

has_single_tcr = chain_counts_beta.query("number_chains == 1").copy()
has_single_tcr["VDJ"] = "TCR"

# Get TCR/BCR doublets

has_trb = list(chain_counts_beta.sample_cell_id)

has_igh = list(chain_counts_heavy.sample_cell_id)

tcr_bcr_doublet = pd.DataFrame(list(set(has_trb) & set(has_igh)))

tcr_bcr_doublet = tcr_bcr_doublet.rename({0:"sample_cell_id"}, axis = 1)

tcr_bcr_doublet["sample"] = tcr_bcr_doublet["sample_cell_id"].str.rsplit("_", n=1).str[0]
tcr_bcr_doublet["barcode"] = tcr_bcr_doublet["sample_cell_id"].str.rsplit("_", n=1).str[1]

tcr_bcr_doublet = pd.merge(tcr_bcr_doublet, batch_info)

tcr_bcr_doublet.to_csv("tcr_bcr_doublet_barcodes.csv", index = False)

# Get cells with a single heavy or beta chain

has_single_tcr = has_single_tcr.drop("tcr_doublet",  axis =1)
has_single_tcr

has_single_bcr = has_single_bcr.drop("BCR_doublet",  axis =1)
has_single_bcr

has_single_bcr_tcr = pd.concat([has_single_bcr, has_single_tcr])

has_single_bcr_or_tcr = has_single_bcr_tcr.copy()

has_single_bcr_or_tcr = has_single_bcr_or_tcr[~has_single_bcr_or_tcr.sample_cell_id.isin(list(set(has_trb) & set(has_igh)))]

has_single_bcr_or_tcr.to_csv("has_tcr_bcr_barcodes.csv", index = False)

check_a = list(has_single_bcr_or_tcr.sample_cell_id)

check_b = list(tcr_bcr_doublet.sample_cell_id)

len(list(set(check_a) & set(check_b)))
#0

# Get single heavy BCR chain info for simple annotation

has_single_bcr_barcode = has_single_bcr_or_tcr.query("VDJ == 'BCR'")['sample_cell_id'].copy()

bcr_df_heavy_single_bcr  = bcr_df.query("(sample_cell_id in @has_single_bcr_barcode) & (locus == 'IGH')").copy()
bcr_df_heavy_single_bcr["sample"] = bcr_df_heavy_single_bcr["sample_cell_id"].str.rsplit("_", n=1).str[0]
bcr_df_heavy_single_bcr["barcode"] = bcr_df_heavy_single_bcr["sample_cell_id"].str.rsplit("_", n=1).str[1]

bcr_df_heavy_single_bcr = bcr_df_heavy_single_bcr.dropna(subset = ["sequence_alignment", "germline_alignment"]).copy()


bcr_df_heavy_single_bcr['mut_freq_v'] = bcr_df_heavy_single_bcr.apply(lambda row : get_mut_freq_v(row['sequence_alignment'],
                 row['germline_alignment'], row["v_germline_end"]), axis = 1)

bcr_df_heavy_single_bcr.to_csv("bcr_df_heavy_single_bcr.csv", index = False)