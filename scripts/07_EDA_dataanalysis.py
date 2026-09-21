import pandas as pd

# read the csv file 

import pandas as pd

df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_tcga_brca_dataset_clean.csv"
)

# overview of data structure 

print("\nDataset Shape")
print(df.shape)

print("\nFirst 10 columns")
print(df.columns[:10])

print("\nData Types")
print(df.dtypes.head(15))

# # check missing values
missing = (
    df.isnull()
      .sum()
      .sort_values(ascending=False)
)

print(missing.head(20))

# # Clinical variable summary

print(df["event"].value_counts())

print(df["survival_time"].describe())

# Categorical variables
print(df["demographic.sex_at_birth"].value_counts())

print(df["demographic.race"].value_counts())

print(df["diagnoses.ajcc_pathologic_stage"].value_counts())

#  check duplicate
print(
    "Duplicate patients:",
    df["cases.submitter_id"].duplicated().sum()
)

# # Check for negative survival time

negative_survival = df[df["survival_time"] < 0]

print("Patients with negative survival time:")
print(negative_survival[[
    "cases.submitter_id",
    "survival_time",
    "event"
]])

print("Number of negative survival patients:", len(negative_survival))

# Find patient with negative survival time

negative_survival = df[df["survival_time"] < 0]

print(negative_survival[
    [
        "cases.submitter_id",
        "survival_time",
        "event",
        "demographic.days_to_death",
        "diagnoses.days_to_last_follow_up"
    ]
])

# Remove patients with negative survival time

df = df[df["survival_time"] >= 0]

print("Dataset shape after removing negative survival:")
print(df.shape)
# lets verify
print("Dataset shape:", df.shape)
print("Missing survival_time:", df["survival_time"].isna().sum())
print("Valid survival_time:", df["survival_time"].notna().sum())

# # Remove patient with negative survival time
df = df[df["survival_time"] >= 0]

# Save QC-passed dataset
# df.to_csv(
#     r"C:\Cancer\TCGA-BRCA\results\final_tcga_brca_dataset_QC.csv",
#     index=False
# )

# print("QC dataset saved successfully!")
# print("Final dataset shape:", df.shape)

print("Shape before filter:", df.shape)

print("Missing survival_time:",
      df["survival_time"].isna().sum())

print("Negative survival_time:",
      (df["survival_time"] < 0).sum())


