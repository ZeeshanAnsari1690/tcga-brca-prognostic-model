import pandas as pd
from pathlib import Path

# Project folder and merging both clinical file as well as raseq
project_folder = Path("C:/Cancer/TCGA-BRCA")

# Results folder
results_folder = project_folder / "results"

# Load cleaned RNA-seq data
expression_df = pd.read_csv(results_folder / "expression_matrix_cleaned.csv")

# Load cleaned clinical data
clinical_df = pd.read_csv(results_folder / "clinical_data_cleaned.csv")

# Check dimensions
print("RNA-seq shape:", expression_df.shape)
print("Clinical shape:", clinical_df.shape)

# Display the first 10 column names
print(expression_df.columns[:10])

# First RNA-seq sample UUID
print(expression_df.columns[3])

# First clinical case UUID
print(clinical_df["cases.submitter_id"].iloc[0])

# Check if the RNA-seq UUID exists in the clinical table
print(expression_df.columns[3] in clinical_df["cases.submitter_id"].values)

# # load the sample sheet
import pandas as pd

sample_sheet = pd.read_csv(
    project_folder / "Data" / "clinical_data" / "gdc_sample_sheet.2026-07-22.tsv",
    sep="\t"
)
print(sample_sheet.columns)
print(sample_sheet.head())

# # identify which identifier used as expression matrix
rna_id = expression_df.columns[3]

print("RNA ID:", rna_id)

print("File ID match:",
      rna_id in sample_sheet["File ID"].values)

print("Sample ID match:",
      rna_id in sample_sheet["Sample ID"].values)

print("File Name match:",
      rna_id in sample_sheet["File Name"].values)


# # creating the mapping of CASE ID to file ID of 
# Create File ID -> Case ID mapping

file_to_case = dict(
    zip(sample_sheet["File ID"], sample_sheet["Case ID"])
)

# # Check the mapping
print(list(file_to_case.items())[:5])

# checking barcode id of case id 
print(sample_sheet.columns)
print(sample_sheet[["File ID", "Case ID"]].head())

print(clinical_df.columns)

# # Rename RNA-seq sample columns
expression_df = expression_df.rename(columns=file_to_case)
print(expression_df.columns[:10])


# # # Check the first few column names
print(expression_df.columns[:10])


# # check the columns name 
print(clinical_df.columns.tolist())

# # Set gene_id as the row index
expression_df = expression_df.set_index("gene_id")

# # Remove annotation columns
expression_df = expression_df.drop(columns=["gene_name", "gene_type"])

# # Transpose the matrix
expression_df = expression_df.T

# # Convert index into a column
expression_df.reset_index(inplace=True)

# # Rename the first column
expression_df.rename(columns={"index": "cases.submitter_id"}, inplace=True)

# # Check the result
print(expression_df.head())
print(expression_df.shape)

# creating one RNA file for each patients and rebuilding RNA expression matrix
selected_samples = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\selected_representative_samples.csv"
)

selected_file_ids = set(selected_samples["File ID"])

metadata_cols = ["gene_id", "gene_name", "gene_type"]

selected_columns = metadata_cols + [
    col for col in expression_df.columns
    if col in selected_file_ids
]

expression_df = expression_df[selected_columns]

print(expression_df.shape)




# # VERIFY PATIENT MATCHING


clinical_ids = set(clinical_df["cases.submitter_id"])
expression_ids = set(expression_df["cases.submitter_id"])

matched_ids = clinical_ids & expression_ids
missing_rnaseq = clinical_ids - expression_ids
missing_clinical = expression_ids - clinical_ids

print("\n" + "="*60)
print("PATIENT MATCHING SUMMARY")
print("="*60)

print(f"Clinical Patients        : {len(clinical_ids)}")
print(f"RNA-seq Patients         : {len(expression_ids)}")
print(f"Matched Patients         : {len(matched_ids)}")
print(f"Clinical Only            : {len(missing_rnaseq)}")
print(f"RNA-seq Only             : {len(missing_clinical)}")

# # Save unmatched patients
pd.DataFrame({
    "cases.submitter_id": sorted(missing_rnaseq)
}).to_csv(
    results_folder / "patients_without_rnaseq.csv",
    index=False
)

pd.DataFrame({
    "cases.submitter_id": sorted(missing_clinical)
}).to_csv(
    results_folder / "patients_without_clinical.csv",
    index=False
)

print("\nLists saved in results folder.")




# # # Merge RNA and clinical data
merged_df = pd.merge(
    clinical_df,
    expression_df,
    on="cases.submitter_id",
    how="inner"
)

print("Merged shape:", merged_df.shape)
print(merged_df.head())

# # save final csv file of data sets
import os

print(os.getcwd())
print(os.path.exists("results"))

# # save final file into the system

merged_df.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_tcga_brca_dataset.csv",
    index=False
)

print("Saved successfully!")