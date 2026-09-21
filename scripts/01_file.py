import pandas as pd

# Path to one TCGA RNA-seq file
file_path = r"C:\Cancer\TCGA-BRCA\Data\RNA_seq\fc0e74d2-ab94-48b1-8e31-808bde7a3fe6\3f3b6a7b-94f0-40a3-9f16-32bab2854678.rna_seq.augmented_star_gene_counts.tsv"

# Read the TSV file
df = pd.read_csv(file_path, sep="\t", comment="#")

# Show the first 5 rows
print(df.head())

# Show the dimensions (rows, columns)
print("\nShape:")
print(df.shape)

# Show column names
print("\nColumns:")
print(df.columns)

# Show data types
print("\nData Types:")
print(df.dtypes)
# Show first 10 genes
print("\nFirst 10 genes:")
print(df[['gene_id', 'gene_name', 'gene_type', 'unstranded']].head(10))