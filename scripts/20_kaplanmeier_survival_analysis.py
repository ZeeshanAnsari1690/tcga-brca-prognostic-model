import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


# Define folders


project_folder = Path("C:/Cancer/TCGA-BRCA")
results_folder = project_folder / "results"


# Load dataset with risk scores


data = pd.read_csv(
    results_folder / "dataset_with_risk_scores.csv"
)

print(data.shape)


# Calculate median risk score


median_score = data["Risk_Score"].median()

print("Median Risk Score:", median_score)


# Divide into High- and Low-risk groups


data["Risk_Group"] = data["Risk_Score"].apply(
    lambda x: "High" if x >= median_score else "Low"
)

print(data["Risk_Group"].value_counts())

# Save grouped dataset
data.to_csv(
    results_folder / "dataset_with_risk_groups.csv",
    index=False
)

# --------------------------------------------------
# Split data
# --------------------------------------------------

high = data[data["Risk_Group"] == "High"]
low = data[data["Risk_Group"] == "Low"]

# --------------------------------------------------
# Kaplan–Meier model


from lifelines.plotting import add_at_risk_counts

kmf_high = KaplanMeierFitter()
kmf_low = KaplanMeierFitter()

fig, ax = plt.subplots(figsize=(8,6))

kmf_high.fit(
    durations=high["survival_time"],
    event_observed=high["event"],
    label="High Risk"
)

kmf_low.fit(
    durations=low["survival_time"],
    event_observed=low["event"],
    label="Low Risk"
)

# Plot survival curves with 95% confidence intervals
kmf_high.plot_survival_function(
    ax=ax,
    ci_show=True,
    linewidth=2.5,
    color="#D55E00"   # reddish-orange
)

kmf_low.plot_survival_function(
    ax=ax,
    ci_show=True,
    linewidth=2.5,
    color="#0072B2"   # blue
)


# Log-rank test


results = logrank_test(
    high["survival_time"],
    low["survival_time"],
    event_observed_A=high["event"],
    event_observed_B=low["event"]
)

# Add Number at Risk Table


add_at_risk_counts(
    kmf_high,
    kmf_low,
    ax=ax
)


# Figure formatting


ax.set_title(
    "Kaplan–Meier Overall Survival Curve",
    fontsize=16,
    fontweight="bold"
)

ax.set_xlabel(
    "Survival Time (Days)",
    fontsize=13
)

ax.set_ylabel(
    "Overall Survival Probability",
    fontsize=13
)

ax.tick_params(axis='both', labelsize=11)

ax.grid(True, linestyle="--", alpha=0.4)

ax.text(
    0.05,
    0.10,
    f"Log-rank p = {results.p_value:.3e}",
    transform=ax.transAxes,
    fontsize=12,
    bbox=dict(facecolor="white", edgecolor="black", alpha=0.8)
)

plt.tight_layout()

plt.savefig(
    results_folder / "Kaplan_Meier_Publication.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()

print("\nLog-rank p-value:", results.p_value)