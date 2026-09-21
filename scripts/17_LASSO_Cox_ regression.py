# import required important libraries
import pandas as pd
import numpy as np
import os

from sksurv.util import Surv
from sksurv.linear_model import CoxnetSurvivalAnalysis

from pathlib import Path

# project path 

project_folder = Path("C:/Cancer/TCGA-BRCA")

results_folder = project_folder / "results"

# load std data into the system
data = pd.read_csv(
    results_folder / "final_dataset_scaled.csv"
)

print(data.shape)
print(data.head())

# load significant gene
significant = pd.read_csv(
    results_folder / "significant_genes.csv"
)

print(significant.head())
print(significant.columns.tolist())

# Extract Gene Names

genes = significant["Gene"].tolist()

print("Significant genes:", len(genes))


# keep only available gene 

available_genes = [
    gene for gene in genes
    if gene in data.columns
]

print("Genes Found :", len(available_genes))
print("Genes Missing:", len(genes) - len(available_genes))

# creating features matrix

X = data[available_genes]

# survival object
y = Surv.from_arrays(
    event=data["event"].astype(bool),
    time=data["survival_time"]
)

# check shape of data
print("Samples :", X.shape[0])
print("Genes   :", X.shape[1])

# Fit the LASSO-Cox Model

print("\nStarting LASSO-Cox fitting...")



lasso = CoxnetSurvivalAnalysis(
    l1_ratio=1.0,
    alpha_min_ratio=0.01,
    n_alphas=100,
    max_iter=100000
)

lasso.fit(X, y)

print("LASSO completed successfully!")

# Get coefficient matrix
coef_df = pd.DataFrame(
    lasso.coef_,
    index=available_genes,
    columns=[f"Alpha_{i+1}" for i in range(lasso.coef_.shape[1])]
)

# Save all coefficients
coef_df.to_csv(results_folder / "lasso_coefficients.csv")

# Select genes with non-zero coefficients
selected = coef_df[(coef_df != 0).any(axis=1)]

# Save selected genes
selected.to_csv(results_folder / "lasso_selected_genes.csv")

print(f"Selected genes: {selected.shape[0]}")
print("Files saved successfully!")