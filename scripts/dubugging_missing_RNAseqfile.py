import pandas as pd

# Load files
clinical_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\clinical_data_cleaned.csv"
)

rna_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\gene_expression_variance_filtered.csv"
)

# Keep only clinical patients that have RNA-seq
clinical_matched = clinical_df[
    clinical_df["cases.submitter_id"].isin(rna_df["cases.submitter_id"])
]

print("Clinical matched:", clinical_matched.shape)
print("RNA patients:", rna_df.shape)

#merging of clinical matches of each data sets

final_df = pd.merge(
    clinical_matched,
    rna_df,
    on="cases.submitter_id",
    how="inner"
)

print(final_df.shape)

print(final_df["cases.submitter_id"].nunique())

# Save
final_df.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_dataset_for_model.csv",
    index=False
)

print(final_df.isnull().sum().sum())

# checking final missing values
missing = final_df.isnull().sum()

print(missing[missing > 0].sort_values(ascending=False))

import pandas as pd

# Load your merged dataset
df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_dataset_for_model.csv"
)

# Check missing age
print("Missing age before:",
      df["diagnoses.age_at_diagnosis"].isna().sum())

# Calculate median age
median_age = df["diagnoses.age_at_diagnosis"].median()

print("Median age:", median_age)



# Fill missing age values
df["diagnoses.age_at_diagnosis"] = (
    df["diagnoses.age_at_diagnosis"]
    .fillna(median_age)
)

# Convert age from days to years
df["diagnoses.age_at_diagnosis"] = (
    df["diagnoses.age_at_diagnosis"] / 365.25
)

df["diagnoses.age_at_diagnosis"] = (
    df["diagnoses.age_at_diagnosis"].round(2)
)

# Verify
print("Missing age after:",
      df["diagnoses.age_at_diagnosis"].isna().sum())



print("Median age:", median_age)

# Fill missing age values
df["diagnoses.age_at_diagnosis"] = (
    df["diagnoses.age_at_diagnosis"]
    .fillna(median_age)
)

# Verify
print("Missing age after:",
      df["diagnoses.age_at_diagnosis"].isna().sum())


# Fill missing age values
df["diagnoses.age_at_diagnosis"] = (
    df["diagnoses.age_at_diagnosis"]
    .fillna(median_age)
)




# Remove original survival columns


df.drop(
    columns=[
        "demographic.days_to_death",
        "diagnoses.days_to_last_follow_up"
    ],
    inplace=True
)

print("Columns removed successfully.")

# Save updated dataset
df.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_datasetA_model.csv",
    index=False
)
print(df.head())
print("completed successfully.")


# # First, check the unique values (for categorical variable)

print(df["demographic.sex_at_birth"].unique())

print(df["demographic.race"].unique())

print(df["diagnoses.ajcc_pathologic_stage"].unique())

# missing satges of patients
print(df["diagnoses.ajcc_pathologic_stage"].value_counts(dropna=False))


# # Convert categorical variables into dummy variables

# # Create a copy
model_df = df.copy()

# # Encode categorical variables
model_df = pd.get_dummies(
    model_df,
    columns=[
        "demographic.sex_at_birth",
        "demographic.race"
    ],
    drop_first=True
)

print(model_df[:])


# Clinical columns first
clinical_cols = [
    "cases.submitter_id",
    "diagnoses.age_at_diagnosis",
    "demographic.sex_at_birth_male",
    "demographic.race_asian",
    "demographic.race_black or african american",
    "demographic.race_not reported",
    "demographic.race_white",
    "survival_time",
    "event"
]

# Remaining columns (mostly genes)
gene_cols = [col for col in model_df.columns if col not in clinical_cols]

# Reorder
model_df = model_df[clinical_cols + gene_cols]

# Save

model_df.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_datasetB_model.csv",
    index=False
)

print("Complete encoded dataset saved successfully.")

