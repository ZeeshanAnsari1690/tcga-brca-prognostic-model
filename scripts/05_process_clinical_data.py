# This script processes clinical data for analysis. It reads in raw clinical data files, performs necessary cleaning and transformations, and outputs a processed dataset ready for downstream analysis.

import pandas as pd
from pathlib import Path

# Define the clinical data directory
clinical_folder = Path(r"C:\Cancer\TCGA-BRCA\Data\clinical_Data")

# Find all the clinical data files in the directory
clinical_files = list(clinical_folder.glob("*.tsv"))

# Print the number of clinical data files found
print(f"Number of clinical data files found: {len(clinical_files)}")

# print the names of the clinical data files found

for file in clinical_files:
    print(f"Found clinical data file: {file.name}")


    # read each clinical data file and perform necessary cleaning and transformations
clinical_file = Path (r"C:\Cancer\TCGA-BRCA\Data\clinical_data\clinical.tsv")
clinical_df = pd.read_csv(clinical_file, sep="\t", comment="#",low_memory=False)

print(clinical_df.shape)
print(clinical_df.head())
print(clinical_df.columns.tolist())

# lets check the data types of the columns
keywords = ["age_at_diagnosis", "gender", "race", "death","vital","follow","stage","grade","age","gender","race","case_id"]
for col in keywords:
    if col in clinical_df.columns:
        print(f"{col}: {clinical_df[col].dtype}")
for col in keywords:
    if any(keyword in col.lower() for keyword in keywords):
        print(col)

        # find need columns for analysis from large clinical data file

for col in clinical_df.columns:

    if "case_id" in col.lower():
        print(col)

for col in clinical_df.columns:
    if "vital" in col.lower():
        print(col)

for col in clinical_df.columns:
    if "death" in col.lower():
        print(col)

for col in clinical_df.columns:
    if "follow" in col.lower():
        print(col)

for col in clinical_df.columns:
    if "stage" in col.lower():
        print(col)

for col in clinical_df.columns:
    if "age_at_diagnosis" in col.lower():
        print(col)

        # find out two columns for gender and race from the clinical data file
        for col in clinical_df.columns:
            if "sex" in col.lower():
                print(col)

        for col in clinical_df.columns:
            if "race" in col.lower():
                print(col)

                # creating a new dataframe with selected columns for analysis
clinical_df_selected = clinical_df[["cases.submitter_id", "diagnoses.age_at_diagnosis", "demographic.sex_at_birth", 
"demographic.race", "demographic.vital_status", "demographic.days_to_death", 
"diagnoses.days_to_last_follow_up", "diagnoses.ajcc_pathologic_stage"
    ]
]

print(clinical_df_selected.head())
print(clinical_df_selected.shape)

# check duplicate case_id in the selected dataframe and unique case_id
print("total case_ids:", clinical_df_selected["cases.submitter_id"].shape[0])
print("duplicate case_ids:", clinical_df_selected.duplicated(subset="cases.submitter_id").sum())
print("unique case_ids:", clinical_df_selected["cases.submitter_id"].nunique())

# creating one row per case_id by keeping the first occurrence of each case_id
clinical_df_selected = clinical_df_selected.drop_duplicates(subset="cases.submitter_id", keep="first")
print(clinical_df_selected.shape)
print(clinical_df_selected["cases.submitter_id"].nunique())

# checking the data have any missing values in the selected dataframe
print(clinical_df_selected.isnull().sum())

# checking unusual missing values in the selected dataframe
print(clinical_df_selected["demographic.vital_status"].value_counts())
print(clinical_df_selected["demographic.days_to_death"].head(10))
            
print(clinical_df_selected["diagnoses.days_to_last_follow_up"].head(10))



# lets verify the data types of the columns in the selected dataframe
import numpy as np

clinical_df_selected["demographic.days_to_death"] = (
    clinical_df_selected["demographic.days_to_death"]
    .astype(str)
    .str.strip()
    .replace("--", np.nan)
)

clinical_df_selected["diagnoses.days_to_last_follow_up"] = (
    clinical_df_selected["diagnoses.days_to_last_follow_up"]
    .astype(str)
    .str.strip()
    .replace("--", np.nan)
)

# converts the columns to numeric types
clinical_df_selected["demographic.days_to_death"] = pd.to_numeric(
    clinical_df_selected["demographic.days_to_death"],
    errors="coerce"
)

clinical_df_selected["diagnoses.days_to_last_follow_up"] = pd.to_numeric(
    clinical_df_selected["diagnoses.days_to_last_follow_up"],
    errors="coerce"
)

# print the data types of the columns in the selected dataframe after conversion
print(clinical_df_selected.isnull().sum())
print(clinical_df_selected.dtypes)

# converts the "diagnoses.age_at_diagnosis" column to numeric type, coercing errors to NaN
clinical_df_selected["diagnoses.age_at_diagnosis"] = pd.to_numeric(
    clinical_df_selected["diagnoses.age_at_diagnosis"],
    errors="coerce"
)

print(clinical_df_selected.dtypes)

# create the survival time ,every columns

import numpy as np

clinical_df_selected["survival_time"] = np.where(
    clinical_df_selected["demographic.vital_status"] == "Dead",
    clinical_df_selected["demographic.days_to_death"],
    clinical_df_selected["diagnoses.days_to_last_follow_up"]
)

clinical_df_selected["event"] = np.where(
    clinical_df_selected["demographic.vital_status"] == "Dead",
    1,
    0
)

# check the survival time and event columns
print(clinical_df_selected[["cases.submitter_id", "demographic.vital_status", "demographic.days_to_death",
     "diagnoses.days_to_last_follow_up", "survival_time", "event"]].head(10))

# quality check of dataframe
dead_patients = clinical_df_selected[
    clinical_df_selected["demographic.vital_status"] == "Dead"
]

print(dead_patients[[
    "cases.submitter_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.days_to_last_follow_up",
    "survival_time",
    "event"
]].head(10))

# save clean clinical data as csv file 
clinical_df_selected.to_csv(
    "results/clinical_data_cleaned.csv",
    index=False
)
# Let's check whether your clinical.tsv already contains the TCGA barcode
for col in clinical_df.columns:
    if "submitter" in col.lower():
        print(col)
    # 
for col in clinical_df.columns:
    if "case" in col.lower():
        print(col)