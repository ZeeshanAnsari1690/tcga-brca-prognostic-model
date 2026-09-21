# import libraries
import pandas as pd
from pathlib import Path

# define the rna-seq data directory
rna_folder= Path(r"C:\Cancer\TCGA-BRCA\Data\RNA_seq")

# fine the all the rna-seq files in the directory
tsv_files = list(rna_folder.glob("*/*.tsv"))

print(f"number of rna-seq files found: {len(tsv_files)}")

# read the each rna-seq file and concatenate them into a single dataframe

expression_list = []
genes = None

for i, file in enumerate(tsv_files):

    print(f"Reading file {i+1}/{len(tsv_files)}")

    df = pd.read_csv(file, sep="\t", comment="#")

    df = df[
        ["gene_id", "gene_name", "gene_type", "unstranded"]
    ]

    sample_name = file.parent.name

    df.rename(
        columns={"unstranded": sample_name},
        inplace=True
    )

    if i == 0:
     genes = df[["gene_id", "gene_name", "gene_type"]]

    expression_list.append(df[sample_name])

# concatenate the expression dataframes into a single dataframe
expression_df = pd.concat(expression_list, axis=1)

# add the gene information to the expression dataframe
expression_df = pd.concat([genes, expression_df], axis=1)

# validate the shape of the expression dataframe
print(expression_df.shape)

# Save the complete raw expression matrix
output_file = r"C:\Cancer\TCGA-BRCA\results\expression_matrix_raw.csv"

expression_df.to_csv(output_file, index=False)

print(f"Expression matrix saved successfully!")
print(f"Location: {output_file}")

# Load the saved expression matrix
expression_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\expression_matrix_raw.csv")

