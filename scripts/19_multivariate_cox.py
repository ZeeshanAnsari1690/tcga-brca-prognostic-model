import pandas as pd
from pathlib import Path

from sksurv.util import Surv
from sksurv.linear_model import CoxPHSurvivalAnalysis

project_folder = Path("C:/Cancer/TCGA-BRCA")
results_folder = project_folder / "results"

# Load patient expression data
data = pd.read_csv(results_folder / "final_dataset_scaled.csv")

# Load LASSO-selected genes
lasso_genes = pd.read_csv(
    results_folder / "lasso_final_selected_genes.csv"
)

print(data.shape)
print(lasso_genes.shape)

# Select genes
genes = lasso_genes["Gene"].tolist()

available_genes = [g for g in genes if g in data.columns]

print("Available genes:", available_genes)

X = data[available_genes]

print("Expression matrix shape:", X.shape)

# Create survival object
y = Surv.from_arrays(
    event=data["event"].astype(bool),
    time=data["survival_time"]
)

# Fit multivariate Cox model
cox = CoxPHSurvivalAnalysis()

cox.fit(X, y)

print("Multivariate Cox model fit successfully")

coef_df = pd.DataFrame({
    "Gene": available_genes,
    "Coefficient": cox.coef_
})

print(coef_df)

# -----------------------------
# Calculate Risk Score
# -----------------------------

risk_scores = cox.predict(X)



from sksurv.metrics import concordance_index_censored

# ---------------------------------
# C-index
# ---------------------------------

c_index_result = concordance_index_censored(
    data["event"].astype(bool).values,
    data["survival_time"].values,
    risk_scores
)

c_index = c_index_result[0]

print("\n" + "="*60)
print("Concordance Index (C-index)")
print("="*60)

print(f"C-index: {c_index:.4f}")


cindex_df = pd.DataFrame({
    "Metric": ["C-index"],
    "Value": [c_index]
})

# save the C-index result to a CSV file


cindex_df.to_csv(
    results_folder / "c_index_result.csv",
    index=False
)

print("\nC-index result saved:")
print(results_folder / "c_index_result.csv")



# Add risk score to dataset
data = data.copy()
data["Risk_Score"] = risk_scores

print("\nFirst 5 Risk Scores:")
print(data[["Risk_Score"]].head())

# Save dataset with risk scores
data.to_csv(
    results_folder / "dataset_with_risk_scores.csv",
    index=False
)

print("\nDataset with risk scores saved successfully!")


# Read the saved file back
# Read the saved file
check = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\dataset_with_risk_scores.csv"
)

print("Saved file shape:", check.shape)




# save  the coefficients

coef_df = pd.DataFrame({
    "Gene": available_genes,
    "Coefficient": cox.coef_
})

coef_df.to_csv(
    results_folder / "multivariate_cox_coefficients.csv",
    index=False
)

print(coef_df)



