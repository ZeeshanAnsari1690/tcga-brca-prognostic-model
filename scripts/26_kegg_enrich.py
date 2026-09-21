# KEGG analysis---------------------

# important library

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp

# output directory
output_dir = r"C:\Cancer\TCGA-BRCA\results\KEGG"

print("Output directory:", output_dir)
print("Directory exists:", os.path.exists(output_dir))

os.makedirs(output_dir, exist_ok=True)

# loading of DEG file  
deg_df = pd.read_csv(
    "results\Significant_Protein_Coding_DEGs.csv"
)

# creating gene list
gene_list = deg_df["gene_name"].dropna().unique().tolist()

# check gene list
print(f"Number of genes: {len(gene_list)}")
print(gene_list[:10])


# creating KEGG enr 
kegg_results = gp.enrichr(
    gene_list=gene_list,
    gene_sets="KEGG_2021_Human",
    organism="human",
    outdir=None
)

# extract df from data
kegg_df = kegg_results.results.copy()

# verify enr work
print(kegg_df.head())
print(kegg_df.columns.tolist())
print(kegg_df.shape)


# save result
kegg_df.to_csv(
    os.path.join(output_dir, "KEGG.csv"),
    index=False
)

# Check if enrichment found significant pathways
sig_kegg = kegg_df[kegg_df["Adjusted P-value"] < 0.05]

print(f"Significant pathways: {len(sig_kegg)}")


# Significant KEGG pathways
kegg_sig = kegg_df[kegg_df["Adjusted P-value"] < 0.05].copy()

print(f"Significant pathways: {len(kegg_sig)}")

# ranked top 20 gene

top20_kegg = (
    kegg_sig
    .sort_values("Combined Score", ascending=False)
    .head(20)
    .copy()
)

# Create Gene Count
top20_kegg["Gene Count"] = (
    top20_kegg["Overlap"]
    .str.split("/")
    .str[0]
    .astype(int)
)

# Create Gene Ratio
top20_kegg["Gene Ratio"] = (
    top20_kegg["Gene Count"] / len(gene_list)
)

# bubble size
top20_kegg["Bubble Size"] = (
    top20_kegg["Gene Count"] /
    top20_kegg["Gene Count"].max()
) * 1800

# cal of FDR
import numpy as np

top20_kegg["-log10(FDR)"] = -np.log10(
    top20_kegg["Adjusted P-value"]
)

# sorting for plor 
top20_kegg = (
    top20_kegg
    .sort_values("Combined Score", ascending=True)
)


import textwrap

top20_kegg["Term"] = top20_kegg["Term"].apply(
    lambda x: "\n".join(textwrap.wrap(x, width=35))
)


# Verify the processed table
print(
    top20_kegg[
        [
            "Term",
            "Gene Count",
            "Gene Ratio",
            "-log10(FDR)"
        ]
    ]
)

# plot
plt.figure(figsize=(10,8))

scatter = plt.scatter(
    x=top20_kegg["Gene Ratio"],
    y=top20_kegg["Term"],
    s=top20_kegg["Bubble Size"],
    c=top20_kegg["-log10(FDR)"],
    cmap="viridis",
    edgecolors="black",
    linewidth=0.5
)

plt.xlabel("Gene Ratio", fontsize=12)
plt.ylabel("KEGG Pathway", fontsize=12)
plt.title("KEGG Pathway Enrichment Analysis", fontsize=15, pad=15)

cbar = plt.colorbar(scatter)
cbar.set_label("-log10(FDR)", fontsize=11)

plt.xlim(0, top20_kegg["Gene Ratio"].max()*1.10)

plt.tight_layout()

# save figure
plt.savefig(
    os.path.join(output_dir,"KEGG_dotplot.png"),
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(output_dir,"KEGG_dotplot.pdf"),
    bbox_inches="tight"
)

plt.show()