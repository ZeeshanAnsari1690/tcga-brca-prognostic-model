import pandas as pd
from sklearn.preprocessing import StandardScaler

# Load dataset
df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_datasetB_model.csv"
)

print("Original shape:", df.shape)

# Keep non-feature columns separately
patient_id = df["cases.submitter_id"]
survival_time = df["survival_time"]
event = df["event"]

# Features to scale

X = df.drop(columns=[
    "cases.submitter_id",
    "survival_time",
    "event",                                                            # target, not a feature
    "demographic.vital_status",                                                # redundant with event
    "diagnoses.ajcc_pathologic_stage"
])

print("Feature matrix:", X.shape)

# Standardize
scaler = StandardScaler()

X_scaled = pd.DataFrame(
    scaler.fit_transform(X),
    columns=X.columns
)

# Rebuild dataset
scaled_df = pd.concat(
    [
        patient_id,
        survival_time,
        event,
        X_scaled
    ],
    axis=1
)

print("Scaled dataset:", scaled_df.shape)

# Save
scaled_df.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_dataset_scaled.csv",
    index=False
)

print("Scaling completed successfully!")