# creating new clean data scripts
# loading the free files
import pandas as pd

# Load RNA expression matrix (before transpose)
expression_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\expression_matrix_cleaned.csv"
)

# Load representative samples
selected_samples = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\selected_representative_samples.csv"
)

# Load cleaned clinical data
clinical_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\clinical_data_cleaned.csv"
)

print("Expression matrix:", expression_df.shape)
print("Representative samples:", selected_samples.shape)
print("Clinical data:", clinical_df.shape)

# keep only selected file RNA
# Representative File IDs
selected_file_ids = set(selected_samples["File ID"])

# Metadata columns
metadata_cols = ["gene_id", "gene_name", "gene_type"]

# Keep only representative RNA-seq files
selected_columns = metadata_cols + [
    col for col in expression_df.columns
    if col in selected_file_ids
]

expression_df = expression_df[selected_columns]

print("Filtered RNA matrix:", expression_df.shape)

# print("representative patients:",selected_samples["case ID"].nunique())

# rename file IDs to case IDs
# File ID → Case ID mapping
file_to_case = dict(
    zip(
        selected_samples["File ID"],
        selected_samples["Case ID"]
    )
)

expression_df.rename(
    columns=file_to_case,
    inplace=True
)

print(expression_df.columns[:10])


# transpos
expression_df = expression_df.set_index("gene_id")

expression_df = expression_df.drop(
    columns=["gene_name", "gene_type"]
)

expression_df = expression_df.T

expression_df.reset_index(inplace=True)

expression_df.rename(
    columns={"index": "cases.submitter_id"},
    inplace=True
)

print(expression_df.shape)
print(expression_df.head())


# Merge RNAseq and clinical data
final_df = pd.merge(
    clinical_df,
    expression_df,
    on="cases.submitter_id",
    how="inner"
)

print("Final shape:", final_df.shape)
print("Unique patients:", final_df["cases.submitter_id"].nunique())
print("Duplicate patients:",
      final_df["cases.submitter_id"].duplicated().sum())

# save it
final_df.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_tcga_brca_dataset_clean.csv",
    index=False
)

print("Dataset saved successfully!")