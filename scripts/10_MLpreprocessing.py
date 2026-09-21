# ==========================================================

import os
import pandas as pd
import numpy as np

# Load QC dataset
df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_tcga_brca_dataset_QC.csv"
)

print("Dataset Shape:", df.shape)


# seprates clinical data and gene expression data

clinical_columns = [
    "cases.submitter_id",
    "diagnoses.age_at_diagnosis",
    "demographic.sex_at_birth",
    "demographic.race",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.ajcc_pathologic_stage",
    "survival_time",
    "event"
]

# Create clinical dataframe
clinical_df = df[clinical_columns]


# Keep patient ID in gene expression data
gene_columns = [
    col for col in df.columns
    if col not in clinical_columns[1:]   # Exclude everything except cases.submitter_id
]

gene_df = df[gene_columns]

# check gene expression data 
print("\nGene Expression Summary")
print("-------------------------")

print("Missing values :", gene_df.isnull().sum().sum())

# Check infinite values only in gene columns

numeric_gene_df = gene_df.drop(columns=["cases.submitter_id"])

print(
    "Infinite values :",
    np.isinf(numeric_gene_df.to_numpy()).sum()
)

print("Duplicate columns :",
      gene_df.columns.duplicated().sum())

print("\nData Types")
print(numeric_gene_df.dtypes.value_counts())


# save seprates file for future use

output_folder = r"C:\Cancer\TCGA-BRCA\results"

clinical_df.to_csv(
    os.path.join(
        output_folder,
        "clinical_features.csv"
    ),
    index=False
)

print("\nClinical shape:", clinical_df.shape)
print("Gene matrix shape:", gene_df.shape)

gene_df.to_csv(
    os.path.join(output_folder,
                 "gene_expression_matrix.csv"),
    index=False
)

print("\nDatasets saved successfully!")