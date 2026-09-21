import pandas as pd
from pathlib import Path
# check path of folder , which contain desire file in it
rna_folder = Path(r"C:\Cancer\TCGA-BRCA\Data\RNA_seq")

# glob()glob() Think of it as:, "Search for files."
tsv_files = list(rna_folder.glob("*/*.tsv"))

# find out number of file present in the all folder
print("Number of RNA-seq files found:", len(tsv_files))

len(tsv_files)
# Show the first two RNA-seq files
for file in tsv_files[:2]:
    print(file)

#  create an empty variable 
merged_df=None

# start  the loop for all file 
for i, file in enumerate(tsv_files):

    # read one file 
    df = pd.read_csv(file, sep="\t", comment="#")

    # keep the required columns 
    df= df[
        ["gene_id", "gene_name", "gene_type", "unstranded"]
    ]
    # rename the sample column
    sample_name = file.parent.name
    print(f"Processing file {i+1}/{len(tsv_files)}: {sample_name}")
    df.rename(columns={"unstranded": sample_name}, inplace=True)
    # merge the dataframes
if merged_df is None:
    merged_df = df
    print(f"first file:{merged_df.shape}")
else:
    merged_df = pd.merge(
        merged_df,
        df,
        on=["gene_id", "gene_name", "gene_type"]
    )

    print(f"after merging file{i+1}/{len(tsv_files)}: {merged_df.shape}")
    # now end thwe loop and save the final merged dataframe to a CSV file

print(f"Final merged dataframe shape: {merged_df.shape}")    

merged_df.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\expression_matrix.csv",
    index=False
)
# here we can print a message to indicate that the expression matrix has been saved successfully
print("Expression matrix saved successfully!")    
