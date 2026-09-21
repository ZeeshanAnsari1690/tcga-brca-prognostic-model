import pandas as pd
import numpy as np

from pathlib import Path

from sksurv.util import Surv
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

from sklearn.model_selection import KFold


# ============================================================
# PATHS
# ============================================================

project_folder = Path("C:/Cancer/TCGA-BRCA")
results_folder = project_folder / "results"


# ============================================================
# LOAD DATA
# ============================================================

data = pd.read_csv(
    results_folder / "final_dataset_scaled.csv"
)

print("=" * 70)
print("C-INDEX VALIDATION")
print("=" * 70)

print("Dataset shape:", data.shape)


# ============================================================
# FINAL THREE GENES
# ============================================================

genes = [
    "ENSG00000004846",   # ABCB5
    "ENSG00000124343",   # XG
    "ENSG00000269495"    # AC011483.2
]

X = data[genes].copy()

y = Surv.from_arrays(
    event=data["event"].astype(bool).values,
    time=data["survival_time"].values
)

event = data["event"].astype(bool).values
time = data["survival_time"].values


# ============================================================
# ORIGINAL MODEL
# ============================================================

cox = CoxPHSurvivalAnalysis()

cox.fit(X, y)

original_risk = cox.predict(X)

original_cindex = concordance_index_censored(
    event,
    time,
    original_risk
)[0]

print("\nOriginal apparent C-index:")
print(f"{original_cindex:.4f}")


# ============================================================
# BOOTSTRAP OPTIMISM CORRECTION
# ============================================================

B = 1000

rng = np.random.default_rng(42)

optimism_values = []

n = len(data)

print("\nRunning bootstrap...")
print(f"Number of bootstrap samples: {B}")

for b in range(B):

    # Bootstrap sample
    boot_idx = rng.integers(
        0,
        n,
        size=n
    )

    # Bootstrap dataset
    X_boot = X.iloc[boot_idx]
    y_boot = Surv.from_arrays(
        event=event[boot_idx],
        time=time[boot_idx]
    )

    # Fit model in bootstrap sample
    boot_cox = CoxPHSurvivalAnalysis()

    boot_cox.fit(
        X_boot,
        y_boot
    )

    # Prediction in bootstrap sample
    risk_boot = boot_cox.predict(X_boot)

    c_boot = concordance_index_censored(
        event[boot_idx],
        time[boot_idx],
        risk_boot
    )[0]

    # Prediction in ORIGINAL dataset
    risk_original = boot_cox.predict(X)

    c_test = concordance_index_censored(
        event,
        time,
        risk_original
    )[0]

    # Optimism
    optimism = c_boot - c_test

    optimism_values.append(
        optimism
    )

    if (b + 1) % 100 == 0:
        print(
            f"Bootstrap {b + 1}/{B}"
        )


# ============================================================
# OPTIMISM-CORRECTED C-INDEX
# ============================================================

mean_optimism = np.mean(
    optimism_values
)

corrected_cindex = (
    original_cindex
    - mean_optimism
)


print("\n" + "=" * 70)
print("BOOTSTRAP RESULTS")
print("=" * 70)

print(
    f"Apparent C-index       : {original_cindex:.4f}"
)

print(
    f"Mean optimism          : {mean_optimism:.4f}"
)

print(
    f"Optimism-corrected C-index : {corrected_cindex:.4f}"
)


# ============================================================
# BOOTSTRAP DISTRIBUTION
# ============================================================

bootstrap_cindex = []

for b in range(B):

    boot_idx = rng.integers(
        0,
        n,
        size=n
    )

    X_boot = X.iloc[boot_idx]

    y_boot = Surv.from_arrays(
        event=event[boot_idx],
        time=time[boot_idx]
    )

    boot_cox = CoxPHSurvivalAnalysis()

    boot_cox.fit(
        X_boot,
        y_boot
    )

    risk_boot = boot_cox.predict(
        X_boot
    )

    c_boot = concordance_index_censored(
        event[boot_idx],
        time[boot_idx],
        risk_boot
    )[0]

    bootstrap_cindex.append(
        c_boot
    )


# ============================================================
# 95% BOOTSTRAP CI
# ============================================================

lower = np.percentile(
    bootstrap_cindex,
    2.5
)

upper = np.percentile(
    bootstrap_cindex,
    97.5
)


print(
    f"Bootstrap 95% CI       : "
    f"{lower:.4f} - {upper:.4f}"
)


# ============================================================
# SAVE RESULTS
# ============================================================

results = pd.DataFrame({
    "Metric": [
        "Apparent C-index",
        "Mean optimism",
        "Optimism-corrected C-index",
        "Bootstrap 95% CI lower",
        "Bootstrap 95% CI upper"
    ],

    "Value": [
        original_cindex,
        mean_optimism,
        corrected_cindex,
        lower,
        upper
    ]
})


output_file = (
    results_folder /
    "Bootstrap_C_index_results.csv"
)

results.to_csv(
    output_file,
    index=False
)


print("\nResults saved to:")
print(output_file)

print("=" * 70)



# ============================================================
# 5-FOLD CROSS-VALIDATED C-INDEX
# ============================================================

print("\n" + "=" * 70)
print("5-FOLD CROSS-VALIDATED C-INDEX")
print("=" * 70)

kf = KFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

cv_risk = np.zeros(
    len(data)
)

for fold, (train_idx, test_idx) in enumerate(
    kf.split(X),
    start=1
):

    print(
        f"Fold {fold}/5"
    )

    X_train = X.iloc[train_idx]
    X_test = X.iloc[test_idx]

    y_train = Surv.from_arrays(
        event=event[train_idx],
        time=time[train_idx]
    )

    # Fit only on training data
    fold_cox = CoxPHSurvivalAnalysis()

    fold_cox.fit(
        X_train,
        y_train
    )

    # Predict ONLY held-out patients
    cv_risk[test_idx] = fold_cox.predict(
        X_test
    )


# ============================================================
# CROSS-VALIDATED C-INDEX
# ============================================================

cv_cindex = concordance_index_censored(
    event,
    time,
    cv_risk
)[0]


print("\n" + "=" * 70)
print("CROSS-VALIDATION RESULT")
print("=" * 70)

print(
    f"5-fold cross-validated C-index: "
    f"{cv_cindex:.4f}"
)