# creating risks group 
import pandas as pd

# Load data
df = pd.read_csv("results\\expression_with_risk_score.csv")

# Create High- and Low-risk groups
cutoff = df["Risk_Score"].median()
df["Risk_Group"] = ["High" if x > cutoff else "Low"
                    for x in df["Risk_Score"]]

print(df["Risk_Group"].value_counts())

# Identify expression columns
gene_cols = [c for c in df.columns if c.startswith("ENSG")]
print(len(gene_cols))

# seprates high risk and low risks
high =df[df["Risk_Group"] == "High"]
low =df[df["Risk_Group"] == "Low"]

# check the samples size
print(high.shape)
print(low.shape)

# # perform DEA
from scipy.stats import ttest_ind
import pandas as pd
import numpy as np

results = []

for gene in gene_cols:

    high_exp = high[gene]
    low_exp = low[gene]

    # Mean expression
    mean_high = high_exp.mean()
    mean_low = low_exp.mean()

    # log2 Fold Change
    log2fc = np.log2((mean_high + 1) / (mean_low + 1))

    # Welch's t-test
    stat, pvalue = ttest_ind(
        high_exp,
        low_exp,
        equal_var=False,
        nan_policy="omit"
    )

    results.append([
        gene,
        mean_high,
        mean_low,
        log2fc,
        pvalue
    ])

deg = pd.DataFrame(
    results,
    columns=[
        "gene_id",
        "Mean_High",
        "Mean_Low",
        "log2FC",
        "P_value"
    ]
)

print(deg.head())


# Adjust p-values (FDR)
from statsmodels.stats.multitest import multipletests

deg["FDR"] = multipletests(
    deg["P_value"],
    method="fdr_bh"
)[1]

# Classify significant genes
deg["Regulation"] = "Not Significant"

deg.loc[
    (deg["FDR"] < 0.05) &
    (deg["log2FC"] > 1),
    "Regulation"
] = "Up"

deg.loc[
    (deg["FDR"] < 0.05) &
    (deg["log2FC"] < -1),
    "Regulation"
] = "Down"

# check dataframe
print(deg.head())

# save file

print(deg["Regulation"].value_counts())
deg.to_csv(
    "results/DEG_results.csv",
    index=False
)

print("DEG analysis completed.")

# load gene type col data in system

# ============================================================
# Step 16: Annotate Differentially Expressed Genes (DEGs)
# ============================================================

import pandas as pd

# -------------------------------
# Load annotation file
# -------------------------------
annotation = pd.read_csv("results/expression_matrix.csv")

# Keep required columns only
annotation = annotation[
    ["gene_id", "gene_name", "gene_type"]
]

# Remove Ensembl version numbers
annotation["gene_id"] = (
    annotation["gene_id"]
    .astype(str)
    .str.split(".")
    .str[0]
)

# Remove duplicate gene IDs
annotation = annotation.drop_duplicates(
    subset="gene_id",
    keep="first"
)

# -------------------------------
# Load DEG results
# -------------------------------
deg = pd.read_csv("results/DEG_results.csv")

# -------------------------------
# Merge annotation with DEG table
# -------------------------------
deg_annotated = deg.merge(
    annotation,
    on="gene_id",
    how="left"
)

# -------------------------------
# Quality control
# -------------------------------
print("="*60)
print("Annotation Summary")
print("="*60)
print(f"Total genes: {deg_annotated.shape[0]}")
print(f"Missing gene names: {deg_annotated['gene_name'].isna().sum()}")
print(f"Missing gene types: {deg_annotated['gene_type'].isna().sum()}")

# Save annotated table
deg_annotated.to_csv(
    "results/DEG_results_annotated.csv",
    index=False
)

# -------------------------------
# Protein-coding genes
# -------------------------------
protein_deg = deg_annotated[
    deg_annotated["gene_type"] == "protein_coding"
].copy()

print(f"\nProtein-coding genes: {protein_deg.shape[0]}")

protein_deg.to_csv(
    "results/Protein_Coding_DEGs.csv",
    index=False
)

# -------------------------------
# Significant protein-coding DEGs
# -------------------------------
significant_deg = protein_deg[
    (protein_deg["FDR"] < 0.05) &
    (abs(protein_deg["log2FC"]) > 1)
].copy()

print(f"Significant protein-coding DEGs: {significant_deg.shape[0]}")
print("\nRegulation Summary")
print(significant_deg["Regulation"].value_counts())



# Save final DEG list
significant_deg.to_csv(
    "results/Significant_Protein_Coding_DEGs.csv",
    index=False
)

print("\nStep 16 completed successfully.")



