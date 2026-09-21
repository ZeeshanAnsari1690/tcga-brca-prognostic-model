import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

from sksurv.util import Surv
from sksurv.metrics import cumulative_dynamic_auc

# ----------------------------------------
# Define folders
# ----------------------------------------

project_folder = Path("C:/Cancer/TCGA-BRCA")
results_folder = project_folder / "results"

# ----------------------------------------
# Load data
# ----------------------------------------

data = pd.read_csv(
    results_folder / "dataset_with_risk_scores.csv"
)

print(data["survival_time"].min())
print(data["survival_time"].max())

# ----------------------------------------
# Create survival object

y = Surv.from_arrays(
    event=data["event"].astype(bool),
    time=data["survival_time"]
)

# ----------------------------------------
# Risk scores
# ----------------------------------------

risk_scores = data["Risk_Score"].values

# ----------------------------------------
# Evaluation times (days)
# ----------------------------------------

times = [365, 1095, 1825]

# ----------------------------------------
# Calculate time-dependent AUC
# ----------------------------------------

auc, mean_auc = cumulative_dynamic_auc(
    y,
    y,
    risk_scores,
    times
)

# ----------------------------------------
# Print AUC values
# ----------------------------------------

for t, a in zip(times, auc):
    print(f"AUC at {t} days: {a:.3f}")

print(f"\nMean AUC: {mean_auc:.3f}")

# ----------------------------------------
# Plot
# ----------------------------------------

plt.figure(figsize=(7,5))

plt.plot(
    [1, 3, 5],
    auc,
    marker="o",
    linewidth=2.5
)

plt.xticks([1, 3, 5])

plt.ylim(0.5, 1.0)

plt.xlabel("Time (Years)", fontsize=12)
plt.ylabel("AUC", fontsize=12)
plt.title("Time-dependent ROC Analysis", fontsize=14)

plt.grid(True, alpha=0.3)

plt.savefig(
    results_folder / "Time_Dependent_ROC.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()