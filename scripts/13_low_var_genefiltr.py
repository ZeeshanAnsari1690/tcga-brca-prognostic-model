import os
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.feature_selection import VarianceThreshold


# Input and Output Paths


input_file = r"C:\Cancer\TCGA-BRCA\results\gene_expression_no_constant.csv"

output_folder = r"C:\Cancer\TCGA-BRCA\results"

#loading data

print("=" * 60)
print("Loading Gene Expression Dataset...")
print("=" * 60)

gene_df = pd.read_csv(input_file)

print("Dataset Loaded Successfully")
print("Dataset Shape :", gene_df.shape)

# Preserve patient IDs
patient_ids = gene_df[["cases.submitter_id"]]

# Gene expression only
gene_data = gene_df.drop(columns=["cases.submitter_id"])


# Dataset Information
samples_before = gene_data.shape[0]
genes_before = gene_data.shape[1]

print("\nSamples :", samples_before)
print("Genes   :", genes_before)

#Apply threshold variance


print("\nApplying Variance Threshold (0.25)...")

selector = VarianceThreshold(threshold=0.25)

filtered_array = selector.fit_transform(gene_data)

selected_genes = gene_data.columns[selector.get_support()]

filtered_df = pd.DataFrame(
    filtered_array,
    columns=selected_genes
)

# Add patient IDs back
filtered_df.insert(
    0,
    "cases.submitter_id",
    patient_ids.values.flatten()
)

print("Filtering Completed Successfully")

# Data statistics

samples_after = filtered_df.shape[0]
genes_after = filtered_df.shape[1] - 1

removed_genes = genes_before - genes_after

percent_removed = (removed_genes / genes_before) * 100

# display results

print("\nVariance Filter Summary")
print("=" * 60)

print("Samples Before :", samples_before)
print("Samples After  :", samples_after)

print("\nGenes Before   :", genes_before)
print("Genes After    :", genes_after)

print("\nRemoved Genes  :", removed_genes)

print(f"Percent Removed : {percent_removed:.2f}%")


# save results


output_file = os.path.join(
    output_folder,
    "gene_expression_variance_filtered.csv"
)

filtered_df.to_csv(
    output_file,
    index=False
)

print("\nFiltered Dataset Saved Successfully")


# creating summary


fig = plt.figure(figsize=(12,8))

plt.axis("off")

title = "TCGA-BRCA Low Variance Filtering Summary"

summary = f"""
Input Information
------------------------------------------------------------
Input File                  : {os.path.basename(input_file)}
Variance Threshold          : 0.25

Dataset Summary
------------------------------------------------------------
Samples Before              : {samples_before}
Samples After               : {samples_after}

Genes Before                : {genes_before}
Genes After                 : {genes_after}

Removed Genes               : {removed_genes}
Percent Removed             : {percent_removed:.2f}%

Filtering Status
------------------------------------------------------------
✓ Dataset Loaded Successfully
✓ Variance Calculated
✓ Low Variance Genes Removed
✓ Filtered Dataset Saved
✓ Dataset Ready for Univariate Cox Regression

Output Files
------------------------------------------------------------
gene_expression_variance_filtered.csv
variance_filter_report.txt
"""

plt.text(
    0.02,
    0.98,
    title,
    fontsize=18,
    fontweight="bold",
    va="top"
)

plt.text(
    0.02,
    0.90,
    summary,
    fontsize=12,
    family="monospace",
    va="top"
)

figure_file = os.path.join(
    output_folder,
    "variance_filter_summary.png"
)

plt.savefig(
    figure_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Summary Figure Saved Successfully")