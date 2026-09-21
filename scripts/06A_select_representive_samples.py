import pandas as pd

# Load sample sheet
sample_sheet = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\Data\clinical_data\gdc_sample_sheet.2026-07-22.tsv",
    sep="\t"
)

# Select one sample per patient
selected_samples = sample_sheet.drop_duplicates(
    subset="Case ID",
    keep="first"
)

print("Original rows:", len(sample_sheet))
print("Selected rows:", len(selected_samples))
print("Unique patients:", selected_samples["Case ID"].nunique())

# Save representative samples
selected_samples.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\selected_representative_samples.csv",
    index=False
)

print("Representative sample file saved successfully!")


