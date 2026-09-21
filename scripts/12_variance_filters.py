import os
import pandas as pd
from sklearn.feature_selection import VarianceThreshold

# Load gene expression matrix
gene_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\gene_expression_matrix.csv"
)

print("Original Shape:", gene_df.shape)

# Preserve patient IDs
patient_ids = gene_df[["cases.submitter_id"]]

# Gene expression only
gene_data = gene_df.drop(columns=["cases.submitter_id"])

# Remove constant genes (variance = 0)
selector = VarianceThreshold(threshold=0)

filtered_data = selector.fit_transform(gene_data)

filtered_gene_names = gene_data.columns[selector.get_support()]

filtered_df = pd.DataFrame(
    filtered_data,
    columns=filtered_gene_names
)

# Add patient IDs back
filtered_df.insert(
    0,
    "cases.submitter_id",
    patient_ids.values.flatten()
)

print("Filtered Shape:", filtered_df.shape)

# Save filtered dataset
output_path = r"C:\Cancer\TCGA-BRCA\results"

filtered_df.to_csv(
    os.path.join(output_path, "gene_expression_no_constant.csv"),
    index=False
)

print("Constant genes removed successfully!")

print("Genes Removed:",
      gene_data.shape[1] - filtered_data.shape[1])