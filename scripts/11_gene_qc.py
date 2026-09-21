# ==========================================================
# Script Name : 11_gene_qc.py
# Project     : TCGA-BRCA Survival Prediction
# Purpose     : Quality Control of Clean Gene Expression Data

# Import Libraries
# ==========================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Input and Output Paths


input_file = r"C:\Cancer\TCGA-BRCA\results\gene_expression_no_constant.csv"
output_folder = r"C:\Cancer\TCGA-BRCA\results"


# Load Gene Expression Matrix


print("=" * 60)
print("Loading Gene Expression Matrix...")
print("=" * 60)


gene_df = pd.read_csv(input_file)

print("Dataset Loaded Successfully")
print("Dataset Shape :", gene_df.shape)

# Separate patient IDs from numeric gene expression

gene_data = gene_df.drop(columns=["cases.submitter_id"])


# Quality Control Checks


missing = gene_data.isnull().sum().sum()

inf = np.isinf(gene_data.to_numpy()).sum()

duplicates = gene_data.columns.duplicated().sum()

constant = (gene_data.nunique() == 1).sum()

minimum = gene_data.min().min()
maximum = gene_data.max().max()
mean = gene_data.mean().mean()
median = gene_data.median().median()


# Display Results


print("\nQuality Control Summary")
print("=" * 60)

print(f"Missing Values        : {missing}")
print(f"Infinite Values       : {inf}")
print(f"Duplicate Gene Names  : {duplicates}")
print(f"Constant Genes        : {constant}")

print("\nExpression Statistics")
print("-" * 30)
print(f"Minimum Expression : {minimum:.2f}")
print(f"Maximum Expression : {maximum:.2f}")
print(f"Mean Expression    : {mean:.2f}")
print(f"Median Expression  : {median:.2f}")


# Save QC Report


report_file = os.path.join(output_folder, "gene_qc_report.txt")

with open(report_file, "w") as f:

    f.write("TCGA-BRCA Gene Expression Quality Control Report\n")
    f.write("=" * 60 + "\n\n")

    f.write("Input Information\n")
    f.write("-" * 30 + "\n")
    f.write(f"Input File : {os.path.basename(input_file)}\n")
    f.write(f"Dataset Shape : {gene_df.shape}\n\n")

    f.write("Quality Control Checks\n")
    f.write("-" * 30 + "\n")
    f.write(f"Missing Values : {missing}\n")
    f.write(f"Infinite Values : {inf}\n")
    f.write(f"Duplicate Gene Names : {duplicates}\n")
    f.write(f"Constant Genes : {constant}\n\n")

    f.write("Expression Statistics\n")
    f.write("-" * 30 + "\n")
    f.write(f"Minimum Expression : {minimum:.2f}\n")
    f.write(f"Maximum Expression : {maximum:.2f}\n")
    f.write(f"Mean Expression : {mean:.2f}\n")
    f.write(f"Median Expression : {median:.2f}\n")

print("\nText QC Report Saved Successfully.")


# Create QC Summary Figure


fig = plt.figure(figsize=(12, 8))
plt.axis("off")

title = "TCGA-BRCA Gene Expression Quality Control Report"

summary = f"""
Input Information
------------------------------------------------------------
Input File                  : {os.path.basename(input_file)}
Dataset Shape               : {gene_df.shape}

Quality Control Checks
------------------------------------------------------------
Missing Values              : {missing}
Infinite Values             : {inf}
Duplicate Gene Names        : {duplicates}
Constant Genes              : {constant}

Expression Statistics
------------------------------------------------------------
Minimum Expression          : {minimum:.2f}
Maximum Expression          : {maximum:.2f}
Mean Expression             : {mean:.2f}
Median Expression           : {median:.2f}

Quality Control Status
------------------------------------------------------------
✓ Dataset Loaded Successfully
✓ Missing Values Checked
✓ Infinite Values Checked
✓ Duplicate Gene Names Checked
✓ Constant Genes Checked
✓ Dataset Ready for Feature Selection

Output Files
------------------------------------------------------------
gene_expression_no_constant.csv
gene_qc_report.txt
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

figure_file = os.path.join(output_folder, "gene_qc_summary.png")

plt.savefig(
    figure_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("QC Summary Figure Saved Successfully.")

print("\nQuality Control Completed Successfully.")