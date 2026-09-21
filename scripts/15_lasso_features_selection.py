"""
=========================================================
STEP 15
LASSO FEATURE SELECTION FOR TCGA-BRCA SURVIVAL ANALYSIS
=========================================================

Objective:
Reduce FDR-significant genes into a smaller biomarker panel
using LASSO feature selection.

Author:
Mohd Jishan Ansari


"""


# # Import Libraries


import os
import warnings
import numpy as np
import pandas as pd

from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv

warnings.filterwarnings("ignore")

# Create Results Folder


os.makedirs("results", exist_ok=True)


# Load Files

# Load Standardized Dataset

scaled_df = pd.read_csv(
    "results/final_dataset_scaled.csv"
)

significant_df = pd.read_csv(
    "results/significant_genes.csv"
)

print()

print("Standardized Dataset Shape :", scaled_df.shape)
print("Significant Genes          :", significant_df.shape)


# Extract Significant Gene Names


gene_list = significant_df["Gene"].tolist()

print()
print("Genes After FDR :", len(gene_list))


# Keep Only Significant Genes Present in Standardized Dataset

available_genes = [
    gene for gene in gene_list
    if gene in scaled_df.columns
]

missing_genes = [
    gene for gene in gene_list
    if gene not in scaled_df.columns
]

print("Genes Found   :", len(available_genes))
print("Genes Missing :", len(missing_genes))


# Filter Gene Expression Matrix


# Filter Significant Gene Expression Matrix

gene_expression = scaled_df[
    ["cases.submitter_id"] + available_genes
]

print()
print("Filtered Gene Matrix :", gene_expression.shape)


# Select Survival Information

clinical_data = scaled_df[
    [
        "cases.submitter_id",
        "survival_time",
        "event"
    ]
]

print("Clinical Data Shape  :", clinical_data.shape)


# Merge Clinical + Gene Expression


# Create Final Dataset

merged_df = scaled_df[
    ["cases.submitter_id", "survival_time", "event"] + available_genes
].copy()

print()
print("Final Dataset :", merged_df.shape)


# Check Missing Values

missing = merged_df.isnull().sum().sum()

print("Missing Values :", missing)

if missing > 0:
    merged_df = merged_df.dropna()

print("Dataset After NA Removal :", merged_df.shape)


# # Prepare Machine Learning Data


# Prepare Machine Learning Data

X = merged_df.drop(
    columns=[
        "cases.submitter_id",
        "survival_time",
        "event"
    ]
)

y_time = merged_df["survival_time"]

y_event = merged_df["event"]

print()
print("=" * 60)
print("FEATURE MATRIX")
print("=" * 60)

print("Samples :", X.shape[0])
print("Genes   :", X.shape[1])

print()
print("Survival Time Shape :", y_time.shape)
print("Event Shape         :", y_event.shape)


# Data already standardized

X_scaled = X.copy()

print()
print("=" * 60)
print("DATA READY FOR LASSO")
print("=" * 60)

print("Samples :", X_scaled.shape[0])
print("Genes   :", X_scaled.shape[1])

print()
print("=" * 60)
print("PREPROCESSING COMPLETED SUCCESSFULLY")
print("=" * 60)