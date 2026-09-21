import pandas as pd
import numpy as np

# read the raw expression matrix
expression_df =pd.read_csv(r"C:\Cancer\TCGA-BRCA\results\expression_matrix_raw.csv")

# inspects the first few rows of the dataframe
print(expression_df.head())
print(expression_df.shape)
print(expression_df.columns[:10])  # print the first 10 columns

# remove the star summary columns from the dataframe
expression_df = expression_df[expression_df["gene_id"].str.startswith("ENSG", na=False)]

# validate the shape of the cleaned expression dataframe
print(expression_df.shape)
print(expression_df.head())

# check for any missing values in the dataframe
missing_values = expression_df.isnull().sum().sum()
print(missing_values) # print the total number of missing values

# check the gene type 
print(expression_df["gene_type"].value_counts())

# check the duplicate gene names
duplicate_gene_names = expression_df["gene_id"].duplicated().sum()

print(f"duplicate gene IDs: {duplicate_gene_names}")

# remove ensembl gene IDs version numbers from the gene_id column
expression_df["gene_id"] = expression_df["gene_id"].str.split(".").str[0]
print(expression_df.head())

# save the cleaned expression matrix to a new CSV file
output_file = r"C:\Cancer\TCGA-BRCA\results\expression_matrix_cleaned.csv"
expression_df.to_csv(output_file, index=False)
print("Cleaned expression matrix saved successfully!")
print(f"location: {output_file}")
