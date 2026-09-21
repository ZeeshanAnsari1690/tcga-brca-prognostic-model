# performing go enrichment analysis
import os
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp 

# creating an output folder
output_dir = "results/GO"
os.makedirs(output_dir, exist_ok=True)

# reading DEG file
deg_file = "results/Significant_Protein_Coding_DEGs.csv"

deg = pd.read_csv(deg_file)

print(deg.head())

# extract gene symbols from the data
gene_list = (
    deg["gene_name"]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
    .tolist()
)

print(f"Total genes for GO enrichment: {len(gene_list)}")

# go enrichment 
go_bp = gp.enrichr(
    gene_list=gene_list,
    gene_sets="GO_Biological_Process_2023",
    organism="human",
    outdir=None
)

# save the results
bp_results = go_bp.results

bp_results.to_csv(
    os.path.join(output_dir, "GO_BP.csv"),
    index=False
)

print(bp_results.head())

# keep only signi variable from data
bp_sig = bp_results[
    bp_results["Adjusted P-value"] < 0.05
]

print(bp_sig.head())

print("=" * 70)
print("GO BIOLOGICAL PROCESS ENRICHMENT")
print("=" * 70)

print("Total genes submitted:", len(gene_list))
print("Significant BP terms:", len(bp_sig))

print("\nTop 20 significant BP terms:")
print(
    bp_sig[
        ["Term", "Overlap", "P-value", "Adjusted P-value", "Combined Score", "Genes"]
    ]
    .sort_values("Adjusted P-value")
    .head(20)
    .to_string(index=False)
)

# top 20 selecting go terms
top20_bp = bp_sig.sort_values(
    "Adjusted P-value"
).head(20)

# dot plot 
# ---------- Publication-quality GO Dot Plot ----------

# Create plotting dataframe
plot_df = top20_bp.copy()

# Count genes contributing to each GO term
plot_df["Gene Count"] = plot_df["Genes"].apply(
    lambda x: len(str(x).split(";"))
)

# Calculate Gene Ratio
TOTAL_GENES = len(gene_list)
plot_df["Gene Ratio"] = plot_df["Gene Count"] / TOTAL_GENES

# Remove GO IDs from labels
plot_df["Term"] = plot_df["Term"].str.replace(
    r"\s*\(GO:\d+\)",
    "",
    regex=True
)

plt.figure(figsize=(10,8))

scatter = plt.scatter(
    x=plot_df["Gene Ratio"],
    y=plot_df["Term"],
    s=plot_df["Gene Count"] * 20,
    c=plot_df["Adjusted P-value"],
    cmap="RdYlBu_r",
    edgecolors="black",
    linewidth=0.5,
    alpha=0.9
)

plt.xlabel("Gene Ratio", fontsize=13)
plt.ylabel("GO Biological Process", fontsize=13)

plt.title(
    "Top 20 Enriched GO Biological Processes",
    fontsize=15,
    fontweight="bold"
)

cbar = plt.colorbar(scatter)
cbar.set_label("Adjusted P-value", fontsize=12)

plt.xticks(fontsize=11)
plt.yticks(fontsize=10)

plt.tight_layout()

plt.savefig(
    os.path.join(output_dir, "GO_BP_dotplot_publication.png"),
    dpi=600,
    bbox_inches="tight"
)

plt.close()


# ---------- Publication-quality GO Bar Plot ----------

plot_df = top20_bp.copy()

plot_df["Term"] = plot_df["Term"].str.replace(
    r"\s*\(GO:\d+\)",
    "",
    regex=True
)

plot_df = plot_df.sort_values("Combined Score")

plt.figure(figsize=(10,8))

bars = plt.barh(
    plot_df["Term"],
    plot_df["Combined Score"],
    color="steelblue",
    edgecolor="black"
)

plt.xlabel("Combined Score", fontsize=13)
plt.ylabel("GO Biological Process", fontsize=13)

plt.title(
    "Top 20 Enriched GO Biological Processes",
    fontsize=15,
    fontweight="bold"
)

plt.xticks(fontsize=11)
plt.yticks(fontsize=10)

plt.grid(
    axis="x",
    linestyle="--",
    alpha=0.4
)

plt.tight_layout()

plt.savefig(
    os.path.join(output_dir, "GO_BP_barplot_publication.png"),
    dpi=600,
    bbox_inches="tight"
)

plt.close()

# GO Cellular Component (CC)
# --------------------------

# run enrich
cc_results = gp.enrichr(
    gene_list=gene_list,
    gene_sets="GO_Cellular_Component_2023",
    organism="human",
    outdir=None
)


# Extract results
cc_df = cc_results.results.copy()

cc_sig = cc_df[
    cc_df["Adjusted P-value"] < 0.05
].copy()

cc_sig = cc_sig.sort_values("Adjusted P-value")

print(f"Significant GO-CC terms: {len(cc_sig)}")


print(len(cc_df))

print(len(cc_sig))

print(cc_df[["Term", "Adjusted P-value"]].head(10))

print(cc_df.shape)


# save
cc_sig.to_csv(
    os.path.join(
        output_dir,
        "GO_CC.csv"
    ),
    index=False
)

# creating plot 
plot_df = cc_sig.copy()

plot_df["Gene Count"] = plot_df["Genes"].apply(
    lambda x: len(str(x).split(";"))
)

plot_df["Gene Ratio"] = (
    plot_df["Gene Count"] / len(gene_list)
)

plot_df["Term"] = plot_df["Term"].str.replace(
    r"\s*\(GO:\d+\)",
    "",
    regex=True
)

top20_cc = plot_df.head(20)

# dot plot
# -----------------------------
# Publication-quality GO-CC Dot Plot
# -----------------------------

import matplotlib.pyplot as plt

# Prepare dataframe
plot_df = cc_sig.copy()

# Gene count
plot_df["Gene Count"] = plot_df["Genes"].apply(
    lambda x: len(str(x).split(";"))
)

# Gene ratio
plot_df["Gene Ratio"] = (
    plot_df["Gene Count"] / len(gene_list)
)

# Clean GO names
plot_df["Term"] = plot_df["Term"].str.replace(
    r"\s*\(GO:\d+\)",
    "",
    regex=True
)

# Top 20 most significant terms
top20_cc = plot_df.sort_values(
    "Adjusted P-value"
).head(20)

# Create figure
plt.figure(figsize=(10,8))

scatter = plt.scatter(
    x=top20_cc["Gene Ratio"],
    y=top20_cc["Term"],
    s=top20_cc["Gene Count"] * 20,
    c=top20_cc["Adjusted P-value"],
    cmap="viridis",          # publication color map
    edgecolors="black",
    linewidth=0.8,
    alpha=0.9
)

plt.xlabel(
    "Gene Ratio",
    fontsize=12
)

plt.ylabel(
    "GO Cellular Component",
    fontsize=12
)

plt.title(
    "Top 20 Enriched GO Cellular Components",
    fontsize=14,
    fontweight="bold"
)

# Colorbar
cbar = plt.colorbar(scatter)
cbar.set_label(
    "Adjusted P-value",
    fontsize=11
)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.tight_layout()

plt.savefig(
    os.path.join(
        output_dir,
        "GO_CC_dotplot_publication.png"
    ),
    dpi=600,
    bbox_inches="tight"
)

plt.close()



# GO_______________________________________MF

# -------------------------------------------------------
# GO Molecular Function (GO-MF) Enrichment Analysis
# -------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp
import os

# =======================================================
# Step 1: Run GO Molecular Function enrichment
# =======================================================

mf_results = gp.enrichr(
    gene_list=gene_list,
    gene_sets="GO_Molecular_Function_2023",
    organism="human",
    outdir=None
)

# =======================================================
# Step 2: Extract enrichment results
# =======================================================

mf_df = mf_results.results.copy()

# Save complete enrichment table
mf_df.to_csv(
    os.path.join(output_dir, "GO_MF.csv"),
    index=False
)

# =======================================================
# Step 3: Keep significant GO terms
# =======================================================

mf_sig = mf_df[
    mf_df["Adjusted P-value"] < 0.05
].copy()

print(f"Significant GO-MF terms: {len(mf_sig)}")

# =======================================================
# Step 4: Select Top 20 GO terms
# =======================================================

top20_mf = (
    mf_sig
    .sort_values("Combined Score", ascending=False)
    .head(20)
    .copy()
)



# =======================================================
# Step 5: Clean GO term names
# =======================================================

top20_mf["Term"] = top20_mf["Term"].str.replace(
    r"\s*\(GO:\d+\)",
    "",
    regex=True
)

# =======================================================
# Step 6: Calculate Gene Count
# =======================================================

top20_mf["Gene Count"] = (
    top20_mf["Genes"]
    .str.split(";")
    .str.len()
)

# =======================================================
# Step 7: Calculate Gene Ratio
# =======================================================

top20_mf["Gene Ratio"] = (
    top20_mf["Gene Count"] / len(gene_list)
)

top20_mf = top20_mf.sort_values("Gene Ratio")

# =======================================================
# Step 8: Bubble Size
# =======================================================

top20_mf["Bubble Size"] = (
    top20_mf["Gene Count"]
    / top20_mf["Gene Count"].max()
) * 1800

# =======================================================
# Step 9: Calculate -log10(FDR)
# =======================================================

top20_mf["-log10(FDR)"] = (
    -np.log10(top20_mf["Adjusted P-value"])
)

# =======================================================
# Step 10: Sort for better visualization
# =======================================================

top20_mf = top20_mf.sort_values("Gene Ratio")

# =======================================================
# Step 11: Create publication-quality dot plot
# =======================================================

plt.figure(figsize=(10,8))

scatter = plt.scatter(
    x=top20_mf["Gene Ratio"],
    y=top20_mf["Term"],
    s=top20_mf["Bubble Size"],
    c=top20_mf["-log10(FDR)"],
    cmap="viridis",
    edgecolors="black",
    linewidth=0.6,
    alpha=0.9
)

# =======================================================
# Step 12: Labels
# =======================================================

plt.xlabel(
    "Gene Ratio",
    fontsize=12
)

plt.ylabel(
    "GO Molecular Function",
    fontsize=12
)

plt.title(
    "Top 20 Enriched GO Molecular Functions",
    fontsize=14,
    fontweight="bold",
    pad=15
    
)

# =======================================================
# Step 13: Colorbar
# =======================================================

cbar = plt.colorbar(scatter)

cbar.set_label(
    "-log10(FDR)",
    fontsize=11
)

# =======================================================
# Step 14: Formatting
# =======================================================

plt.xlim(
    0,
    top20_mf["Gene Ratio"].max()*1.10
)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.tight_layout()

# =======================================================
# Step 15: Save figures
# =======================================================

plt.savefig(
    os.path.join(
        output_dir,
        "GO_MF_dotplot_publication.png"
    ),
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(
        output_dir,
        "GO_MF_dotplot_publication.pdf"
    ),
    bbox_inches="tight"
)

plt.close()

print("GO Molecular Function analysis completed successfully.")