"""
==========================================================
Project : TCGA-BRCA Survival Prediction
Script  : 18_lasso_feature_selection.py

Objective:
Perform LASSO-Cox feature selection using
5-fold Cross Validation.

Input:
    results/final_dataset_scaled.csv
    results/significant_genes.csv

Output:
    results/lasso_cv_results.csv
    results/lasso_final_selected_genes.csv

Author : Mohd Jishan
==========================================================
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.model_selection import KFold

from sksurv.util import Surv
from sksurv.linear_model import CoxnetSurvivalAnalysis

from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings(
    "ignore",
    category=ConvergenceWarning
)

print("=" * 60)
print("LASSO FEATURE SELECTION")
print("=" * 60)

# ---------------------------------------------------------
# Project folders
# ---------------------------------------------------------

project_folder = Path("C:/Cancer/TCGA-BRCA")
results_folder = project_folder / "results"

# ---------------------------------------------------------
# Load Data
# ---------------------------------------------------------

data = pd.read_csv(
    results_folder / "final_dataset_scaled.csv"
)

significant = pd.read_csv(
    results_folder / "significant_genes.csv"
)

genes = significant["Gene"].tolist()

available_genes = [
    g for g in genes
    if g in data.columns
]

print(f"Available genes : {len(available_genes)}")

X = data[available_genes]

y = Surv.from_arrays(
    event=data["event"].astype(bool),
    time=data["survival_time"]
)

print(f"Samples : {X.shape[0]}")
print(f"Genes   : {X.shape[1]}")
print(X.shape)
print(data["event"].value_counts())

# ---------------------------------------------------------
# # Fit Initial LASSO Path
# ---------------------------------------------------------

print("\nFitting initial LASSO path...")

lasso = CoxnetSurvivalAnalysis(
    l1_ratio=1.0,
    n_alphas=30,
    alpha_min_ratio=0.3,
    max_iter=100000
)

lasso.fit(X, y)

alphas = lasso.alphas_

print("LASSO path completed.")
print("Number of alpha values:", len(alphas))

# ---------------------------------------------------------
# Cross Validation
# ---------------------------------------------------------

cv = KFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

cv_scores = []

print("\nStarting 5-Fold Cross Validation...\n")

for alpha in alphas:

    fold_scores = []

    for train_idx, test_idx in cv.split(X):

        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]

        y_train = y[train_idx]
        y_test = y[test_idx]

        model = CoxnetSurvivalAnalysis(
            l1_ratio=1.0,
            alphas=[alpha],
            max_iter=100000
        )

        try:

            model.fit(X_train, y_train)

            score = model.score(X_test, y_test)

            fold_scores.append(score)

        except Exception as e:

            print(f"Alpha {alpha:.6f} failed : {e}")
            continue

    if len(fold_scores) == cv.get_n_splits():

        mean_score = np.mean(fold_scores)

        cv_scores.append(mean_score)

        print(f"Alpha={alpha:.6f}  CV={mean_score:.4f}")

    else:

        cv_scores.append(np.nan)

print("\nCross Validation Completed!")

# ---------------------------------------------------------
# Best Alpha
# ---------------------------------------------------------

cv_scores = np.array(cv_scores)

valid = ~np.isnan(cv_scores)

if np.sum(valid) == 0:
    raise RuntimeError(
        "No valid alpha values found during cross-validation."
    )

best_index = np.argmax(cv_scores[valid])

best_alpha = alphas[valid][best_index]

best_score = cv_scores[valid][best_index]

print(f"Best Alpha : {best_alpha:.6f}")
print(f"Best Score : {best_score:.4f}")

# ---------------------------------------------------------
# Save CV Results
# ---------------------------------------------------------

cv_results = pd.DataFrame({

    "Alpha": alphas,
    "CV_Score": cv_scores

})

cv_results.to_csv(
    results_folder / "lasso_cv_results.csv",
    index=False
)

print("CV results saved.")

# ---------------------------------------------------------
# Final LASSO Model
# ---------------------------------------------------------

print("\nFitting final LASSO model...")

final_model = CoxnetSurvivalAnalysis(
    l1_ratio=1.0,
    alphas=[best_alpha],
    max_iter=100000
)

final_model.fit(X, y)

final_model.fit(X, y)

coef = final_model.coef_.ravel()

selected_mask = coef != 0

selected_genes = pd.DataFrame({

    "Gene": np.array(available_genes)[selected_mask],
    "Coefficient": coef[selected_mask]

})

selected_genes.to_csv(
    results_folder / "lasso_final_selected_genes.csv",
    index=False
)

print("\nSelected genes:", len(selected_genes))

print(selected_genes.head())

print("\nSelected genes saved.")

print("\n" + "=" * 60)
print("LASSO FEATURE SELECTION COMPLETED SUCCESSFULLY")
print("=" * 60)