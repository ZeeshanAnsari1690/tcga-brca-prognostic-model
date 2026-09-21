"""
==========================================================
Project : TCGA-BRCA Survival Prediction
Script  : 14_univariate_cox.py

Objective:
Perform Univariate Cox Proportional Hazards Regression
for each gene to identify genes significantly associated
with patient overall survival.

Input:
    results/gene_expression_variance_filtered.csv
    results/clinical_features.csv

Output:
    results/univariate_cox_results.csv
    results/significant_genes.csv
    results/univariate_cox_report.txt
    results/volcano_plot.png

Author : Your Name 
 Mohd Jishan
==========================================================
"""
# Import Libraries



import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lifelines

from lifelines import CoxPHFitter

from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


# Create Output Directory


os.makedirs("results", exist_ok=True)

# Load Datasets

import pandas as pd

merged_df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_datasetB_model.csv"
)

print(merged_df.columns[:15])




# # Identify Gene Columns


exclude_columns = [
    "cases.submitter_id",
    "diagnoses.age_at_diagnosis",
    "demographic.sex_at_birth_male",
    "demographic.race_asian",
    "demographic.race_black or african american",
    "demographic.race_not reported",
    "demographic.race_white",
    "demographic.vital_status",
    "diagnoses.ajcc_pathologic_stage",
    "survival_time",
    "event"
]

gene_columns = [
    col
    for col in merged_df.columns
    if col not in exclude_columns
]

print("Genes to Analyze :", len(gene_columns))


non_gene = [
    col for col in merged_df.columns
    if not col.startswith("ENSG")
]

print(non_gene)
print(len(non_gene))


# # Run Univariate Cox Regression

print("\nRunning Univariate Cox Regression...")

results = []

cox = CoxPHFitter()

for i, gene in enumerate(gene_columns):

    temp = merged_df[
        ["survival_time", "event", gene]
    ].copy()

    temp = temp.dropna()

    try:

        cox.fit(
            temp,
            duration_col="survival_time",
            event_col="event"
        )

        summary = cox.summary.loc[gene]

        results.append({

            "Gene": gene,

            "Hazard_Ratio":
            summary["exp(coef)"],

            "Lower_95CI":
            summary["exp(coef) lower 95%"],

            "Upper_95CI":
            summary["exp(coef) upper 95%"],

            "Coefficient":
            summary["coef"],

            "P_value":
            summary["p"]

        })

    except Exception:
        continue

    if (i + 1) % 1000 == 0:
        print(f"{i+1} genes completed...")


# # Create Results DataFrame


results_df = pd.DataFrame(results)

results_df = results_df.sort_values(
    by="P_value"
)


# Multiple testing correction (Benjamini-Hochberg FDR)

results_df["FDR"] = multipletests(
    results_df["P_value"],
    method="fdr_bh"
)[1]


# # Save Complete Results


results_df.to_csv(
    "results/univariate_cox_results.csv",
    index=False
)

print("\nComplete Results Saved.")


# # Significant Genes

significant = results_df[
    results_df["FDR"] < 0.05
]

significant.to_csv(
    "results/significant_genes.csv",
    index=False
)

print("Significant Genes :", len(significant))


# # Volcano Plot


results_df["logHR"] = np.log2(
    results_df["Hazard_Ratio"]
)

results_df["minus_log10_p"] = -np.log10(
    results_df["P_value"]
)

plt.figure(figsize=(10,7))

plt.scatter(
    results_df["logHR"],
    results_df["minus_log10_p"],
    s=8,
    alpha=0.6
)

plt.axhline(
    -np.log10(0.05),
    color="red",
    linestyle="--"
)

plt.xlabel("log2(Hazard Ratio)")
plt.ylabel("-log10(P-value)")
plt.title("Univariate Cox Volcano Plot")

plt.tight_layout()

plt.savefig(
    "results/volcano_plot01.png",
    dpi=300
)

plt.close()

print("Volcano Plot Saved.")


# Report

with open(
    "results/univariate_cox_report01.txt",
    "w"
) as f:

    f.write("="*50 + "\n")
    f.write("UNIVARIATE COX REGRESSION REPORT\n")
    f.write("="*50 + "\n\n")

    f.write(f"Patients : {merged_df.shape[0]}\n")
    f.write(f"Genes Analyzed : {len(gene_columns)}\n")
    f.write(f"Successful Models : {len(results_df)}\n")
    f.write(f"Significant Genes (FDR<0.05): {len(significant)}\n")

print("Report Saved.")

# Finish

print("\n" + "=" * 60)
print("UNIVARIATE COX ANALYSIS COMPLETED SUCCESSFULLY")
print("=" * 60)