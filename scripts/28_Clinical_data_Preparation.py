

import pandas as pd
import numpy as np

# ==========================================================
# Load Dataset
# ==========================================================

df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\clinical_characteristics_dataset.csv"
)

print("="*60)
print("Original Dataset")
print("="*60)
print(df.shape)

# ==========================================================
# Replace Missing Values
# ==========================================================

missing_values = [
    "'--",
    "--",
    "not reported",
    "Not Reported",
    "unknown",
    "Unknown",
    ""
]

df.replace(missing_values, np.nan, inplace=True)

# ==========================================================
# Clean T Stage
# ==========================================================

df["diagnoses.ajcc_pathologic_t"] = (
    df["diagnoses.ajcc_pathologic_t"]
    .replace({
        "T1a":"T1",
        "T1b":"T1",
        "T1c":"T1",

        "T2a":"T2",
        "T2b":"T2",

        "T3a":"T3",

        "T4a":"T4",
        "T4b":"T4",
        "T4c":"T4",
        "T4d":"T4",

        "Tis (DCIS)":"Tis",
        "Tis (LCIS)":"Tis"
    })
)

# ==========================================================
# Clean N Stage
# ==========================================================

df["diagnoses.ajcc_pathologic_n"] = (
    df["diagnoses.ajcc_pathologic_n"]
    .replace({
        "N0 (mol+)":"N0",
        "N0 (i+)":"N0",
        "N0 (i-)":"N0",

        "N1a":"N1",
        "N1b":"N1",
        "N1c":"N1",
        "N1mi":"N1",

        "N2a":"N2",

        "N3a":"N3",
        "N3b":"N3"
    })
)

# ==========================================================
# Clean M Stage
# ==========================================================

df["diagnoses.ajcc_pathologic_m"] = (
    df["diagnoses.ajcc_pathologic_m"]
    .replace({
        "cM0 (i+)":"M0"
    })
)

# ==========================================================
# Clean AJCC Stage
# ==========================================================

df["diagnoses.ajcc_pathologic_stage"] = (
    df["diagnoses.ajcc_pathologic_stage"]
    .str.replace("Stage ", "", regex=False)
)

# ==========================================================
# Convert Age
# ==========================================================

df["diagnoses.age_at_diagnosis"] = pd.to_numeric(
    df["diagnoses.age_at_diagnosis"],
    errors="coerce"
)

# ==========================================================
# Display Summary
# ==========================================================

print("\nUnique AJCC Stage")
print(sorted(df["diagnoses.ajcc_pathologic_stage"].dropna().unique()))

print("\nUnique T Stage")
print(sorted(df["diagnoses.ajcc_pathologic_t"].dropna().unique()))

print("\nUnique N Stage")
print(sorted(df["diagnoses.ajcc_pathologic_n"].dropna().unique()))

print("\nUnique M Stage")
print(sorted(df["diagnoses.ajcc_pathologic_m"].dropna().unique()))

print("\nMissing Values")
print(df[
    [
        "diagnoses.ajcc_pathologic_stage",
        "diagnoses.ajcc_pathologic_t",
        "diagnoses.ajcc_pathologic_n",
        "diagnoses.ajcc_pathologic_m"
    ]
].isna().sum())

# ==========================================================
# Save
# ==========================================================

output = r"C:\Cancer\TCGA-BRCA\results\clinical_characteristics_dataset_clean.csv"

df.to_csv(output,index=False)

print("\n"+"="*60)
print("Clinical variables cleaned successfully.")
print("Saved to:")
print(output)
print("="*60)

