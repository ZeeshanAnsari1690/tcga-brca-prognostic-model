# performing gene sets enrichments analysis


# load import library

import os
import pandas as pd
import numpy as np
import gseapy as gp

# Output directory
output_dir = r"C:\Cancer\TCGA-BRCA\results\GSEA"
os.makedirs(output_dir, exist_ok=True)

# Load annotated DEG results
deg_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\DEG_results_annotated.csv"
)

# keep only proteins coding gene
deg_df = deg_df[
    deg_df["gene_type"] == "protein_coding"
].copy()

print(f"Protein-coding genes: {deg_df.shape[0]}")

# Creating  ranked gene list

rank_df = deg_df[
    ["gene_name", "log2FC"]
].copy()

rank_df = rank_df.dropna()

rank_df = rank_df.drop_duplicates(
    subset="gene_name"
)

rank_df = rank_df.sort_values(
    by="log2FC",
    ascending=False
)

#check data frame
print(rank_df.head())

print(rank_df.tail())

print(rank_df.shape)

# save

rank_file = os.path.join(
    output_dir,
    "Ranked_Genes.rnk"
)

rank_df.to_csv(
    rank_file,
    sep="\t",
    header=False,
    index=False
)

print("Ranked gene list saved.")

# Run GSEA (Preranked)

pre_res = gp.prerank(
    rnk=rank_file,
    gene_sets="MSigDB_Hallmark_2020",
    threads=4,
    permutation_num=1000,
    min_size=15,
    max_size=500,
    seed=42,
    outdir=output_dir,
    verbose=True
)

# check df
print(pre_res.res2d.head())
print(pre_res.res2d.head())

# filters significant pathaway
# Use the standard GSEA criteria from the Broad Institute:
# FDR q-value < 0.25
# Nominal P-value < 0.05

# Copy results
gsea_results = pre_res.res2d.copy()

# Significant pathways
sig_gsea = gsea_results[
    (gsea_results["FDR q-val"] < 0.25) &
    (gsea_results["NOM p-val"] < 0.05)
].copy()

# Sort by normalized enrichment score
sig_gsea = sig_gsea.sort_values(
    "NES",
    ascending=False
)

print(sig_gsea.shape)
print(sig_gsea[["Term", "NES", "NOM p-val", "FDR q-val"]])

# save
sig_gsea.to_csv(
    os.path.join(output_dir, "Significant_Hallmark_GSEA.csv"),
    index=False
)

# =============================================================================
# GSEApy Plot Debugging Script


# continue to complete this plot 