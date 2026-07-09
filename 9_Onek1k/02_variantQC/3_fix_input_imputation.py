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
# # Code for Preaparing the imputation input
# - listing the snps (allele swap and strand flip) that should be removed based on the qc report from Imputation server

# %%
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import re

# %%
os.getcwd()

# %%
os.chdir('/path/genotype_oneK1K_plink2_v6_correctdata/TOPMed_rerun_2')

# %%
report = pd.read_csv("snps-excluded.txt", sep="\t")
report

# %%
#no NA
report[report['INFO'].isna()]

# %%
report.INFO.value_counts()

# %%
report['INFO'].str.contains('Strand flip', na=False).sum()

# %%
report['INFO'].str.contains('Allele switch', na=False).sum()


# %%
#checking for palindromes - to be removed

# %%
def is_palindromic(ref, alt):
    palindromes = [('A','T'), ('T','A'), ('G','C'), ('C','G')]
    return (ref, alt) in palindromes
report.apply(lambda row: is_palindromic(row['REF'], row['ALT']), axis=1).value_counts()

# %%
report = report[~report.apply(lambda row: is_palindromic(row['REF'], row['ALT']), axis=1)]

# %%
report.shape

# %%
print(report['INFO'].str.contains('Strand flip|Allele switch', na=False).sum())
print(report['INFO'].str.contains('Strand flip', na=False).sum())
print(report['INFO'].str.contains('Allele switch', na=False).sum())

# %%
report.to_csv("snps-excluded_nopalindromic.txt", sep="\t", index = False)

# %%
#take the rsID for strand flipping

# %%
strand_flip_rsids = report.loc[report['INFO'].str.contains('Strand flip', na=False),'ID']
strand_flip_rsids

# %%
strand_flip_rsids.to_csv(
    "Strand-Flip.txt",
    index=False,
    header=False
)

# %%
#take the rsID + allele for allele swap/switch

# %%
input_file = "snps-excluded_nopalindromic.txt"
output_file = "Force-Allele.txt"
with open(input_file) as fin, open(output_file, "w") as fout:
    next(fin)  # skip header

    for line in fin:
        parts = line.strip().split()

        rsid = parts[0]              # rsID
        info = " ".join(parts[5:])   # rebuild INFO field

        # Only allele-switch cases
        if "Allele switch" not in info:
            continue

        # Extract first allele after "Reference Panel:"
        match = re.search(r"Reference Panel:\s*([ACGT])/[ACGT]", info)
        if match:
            ref_allele = match.group(1)
            fout.write(f"{rsid}\t{ref_allele}\n")

# %%
allele_swap_rsids = pd.read_csv("Force-Allele.txt",
                              header = None, sep="\t")[0].tolist()
len(allele_swap_rsids)

