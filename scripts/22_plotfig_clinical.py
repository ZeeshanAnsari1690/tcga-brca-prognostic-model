# import important library 
from scipy.stats import zscore
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# read data from csv file 
# Read risk score file  ----------------for risks score plot

df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\dataset_with_risk_scores.csv"
)

# Check columns
print(df.columns)

# Display first rows
print(df.head())

# Sort patients
# Sort by risk score
df = df.sort_values("Risk_Score").reset_index(drop=True)

# Median cutoff
median = df["Risk_Score"].median()

# Assign risk groups
df["Risk_Group"] = np.where(
    df["Risk_Score"] >= median,
    "High",
    "Low"
)

print(df.head())

# count patients
low_count = (df["Risk_Group"] == "Low").sum()
high_count = (df["Risk_Group"] == "High").sum()

print("Low-risk:", low_count)
print("High-risk:", high_count)

# creating plot by clinical
plt.figure(figsize=(11,4.5))

# Plot low-risk patients
plt.scatter(
    df.index[df["Risk_Group"]=="Low"],
    df.loc[df["Risk_Group"]=="Low","Risk_Score"],
    color="#2C7BB6",
    s=20,
    label="Low Risk"
)

# Plot high-risk patients
plt.scatter(
    df.index[df["Risk_Group"]=="High"],
    df.loc[df["Risk_Group"]=="High","Risk_Score"],
    color="#D7191C",
    s=20,
    label="High Risk"
)

# Connect points
plt.plot(
    df.index,
    df["Risk_Score"],
    color="gray",
    linewidth=0.8,
    alpha=0.6
)

# Median cutoff
plt.axvline(
    x=low_count,
    color="black",
    linestyle="--",
    linewidth=1.5
)

plt.xlabel("Patients (ordered by Risk Score)", fontsize=12)
plt.ylabel("Risk Score", fontsize=12)
plt.title(
    "Risk Score Distribution",
    fontsize=14,
    fontweight="bold"
)

plt.legend(frameon=False)

plt.grid(False)

plt.tight_layout()

# saving figures
plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Risk_Score_Distribution.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Risk_Score_Distribution.tiff",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

# survival status of patients

# Sort patients by Risk Score


df = df.sort_values("Risk_Score").reset_index(drop=True)

# Median cutoff
median = df["Risk_Score"].median()

df["Risk_Group"] = np.where(
    df["Risk_Score"] >= median,
    "High",
    "Low"
)

# Number of low-risk patients
low_count = (df["Risk_Group"] == "Low").sum()

# ==========================================================
# Separate alive and dead patients
# event:
# 0 = Alive (Censored)
# 1 = Dead
# ==========================================================

alive = df[df["event"] == 0]
dead = df[df["event"] == 1]

# ==========================================================
# Create Figure
# ==========================================================

fig, ax = plt.subplots(figsize=(12, 5))

# Alive patients
ax.scatter(
    alive.index,
    alive["survival_time"],
    color="#4C72B0",
    s=18,
    label="Alive"
)

# Dead patients
ax.scatter(
    dead.index,
    dead["survival_time"],
    color="#C44E52",
    s=18,
    label="Dead"
)

# Median cutoff
ax.axvline(
    x=low_count,
    color="black",
    linestyle="--",
    linewidth=1.5
)

# ==========================================================
# Labels
# ==========================================================

ax.set_xlabel(
    "Patients (ordered by Risk Score)",
    fontsize=13
)

ax.set_ylabel(
    "Overall Survival (Days)",
    fontsize=13
)

ax.set_title(
    "Survival Status",
    fontsize=16,
    fontweight="bold"
)

# ==========================================================
# Beautify
# ==========================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(labelsize=11)

ax.legend(
    frameon=False,
    fontsize=11,
    loc="upper right"
)

# Low / High labels
ymax = df["survival_time"].max()

ax.text(
    low_count/2,
    ymax*1.02,
    "Low Risk",
    ha="center",
    fontsize=12,
    color="#4C72B0",
    fontweight="bold"
)

ax.text(
    low_count+(len(df)-low_count)/2,
    ymax*1.02,
    "High Risk",
    ha="center",
    fontsize=12,
    color="#C44E52",
    fontweight="bold"
)

plt.tight_layout()

# ==========================================================
# Save
# ==========================================================

plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Survival_Status.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Survival_Status.tiff",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

# heatmap for gene by risk score

risk = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\dataset_with_risk_scores.csv"
)

expr = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\gene_expression_variance_filtered.csv"
)


print(expr.columns[:5])



# Remove duplicate patient IDs
expr = expr.drop_duplicates(subset="cases.submitter_id")
risk = risk.drop_duplicates(subset="cases.submitter_id")

risk =risk[["cases.submitter_id","Risk_Score"]]

# Merge
merged = expr.merge(
    risk,
    on="cases.submitter_id",
    how="inner",
    validate="one_to_one"
)

print("Expression:", expr.shape)
print("Risk:", risk.shape)
print("Merged:", merged.shape)

# Save
merged.to_csv(
    r"C:\Cancer\TCGA-BRCA\results\expression_with_risk_score.csv",
    index=False
)

# Sort patients by Risk Score
# -----------------------------------------
merged = merged.sort_values("Risk_Score").reset_index(drop=True)

# -----------------------------------------
# Signature genes
# -----------------------------------------
genes = [
    "ENSG00000004846",
    "ENSG00000124343",
    "ENSG00000269495"
]

# -----------------------------------------
# Extract expression matrix
# -----------------------------------------
heatmap = merged[genes].copy()

# -----------------------------------------
# Gene-wise Z-score normalization
# (each gene normalized across all patients)
# -----------------------------------------
heatmap = heatmap.apply(
    lambda x: (x - x.mean()) / x.std(),
    axis=0
)

# Heatmap requires genes as rows
heatmap = heatmap.T

# -----------------------------------------
# Plot
# -----------------------------------------
plt.figure(figsize=(12, 2.8))

ax = sns.heatmap(
    heatmap,
    cmap="RdBu_r",
    center=0,
    vmin=-2,
    vmax=2,
    xticklabels=False,
    yticklabels=True,
    linewidths=0,
    cbar_kws={
        "label": "Z-score",
        "shrink": 0.8
    }
)

# Median risk separator
plt.axvline(
    x=len(merged)//2,
    color="black",
    linestyle="--",
    linewidth=1.5
)

plt.xlabel("Patients (Low Risk → High Risk)", fontsize=12)
plt.ylabel("Genes", fontsize=12)
plt.title("Expression of Prognostic Signature Genes", fontsize=14)

plt.yticks(rotation=0)

plt.tight_layout()

plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Heatmap.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Heatmap.tiff",
    dpi=600,
    bbox_inches="tight"
)

plt.show()