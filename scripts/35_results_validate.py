import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os

# ============================================================
# 1. FILE
# ============================================================

file = "results/Clinical_Characteristics_Table.csv"

df = pd.read_csv(file)

print("=" * 70)
print("CLINICAL CHARACTERISTICS BY RISK GROUP")
print("=" * 70)

print("Columns:")
print(df.columns.tolist())


# ============================================================
# 2. CLEAN COLUMN NAMES
# ============================================================

df.columns = df.columns.str.strip()

# Keep original rows
data = df.copy()


# ============================================================
# 3. EXTRACT PERCENTAGE
# ============================================================

def extract_percent(value):

    if pd.isna(value):
        return np.nan

    value = str(value)

    match = re.search(r"\(([-+]?\d*\.?\d+)%\)", value)

    if match:
        return float(match.group(1))

    return np.nan


# ============================================================
# 4. EXTRACT COUNT
# ============================================================

def extract_count(value):

    if pd.isna(value):
        return np.nan

    value = str(value)

    match = re.search(r"^(\d+)", value)

    if match:
        return int(match.group(1))

    return np.nan


# ============================================================
# 5. EXTRACT MEDIAN AGE
# ============================================================

def extract_age(value):

    if pd.isna(value):
        return np.nan

    value = str(value)

    match = re.match(r"([\d.]+)", value)

    if match:
        return float(match.group(1))

    return np.nan


# ============================================================
# 6. BUILD CLEAN CLINICAL DATA
# ============================================================

records = []

current_group = None

for _, row in data.iterrows():

    characteristic = row["Characteristic"]

    # --------------------------------------------------------
    # If characteristic is present, identify the group
    # --------------------------------------------------------

    if pd.notna(characteristic):

        characteristic = str(characteristic).strip()

        if characteristic.startswith("Age"):
            current_group = "Age"

        elif characteristic.startswith("Male"):
            current_group = "Sex"

        elif characteristic.startswith("White"):
            current_group = "Race"

        elif characteristic.startswith("Black"):
            current_group = "Race"

        elif characteristic.startswith("Asian"):
            current_group = "Race"

        elif characteristic.startswith("Race Not"):
            current_group = "Race"

        elif characteristic == "AJCC Stage":
            current_group = "AJCC Stage"

        elif characteristic == "T Stage":
            current_group = "T Stage"

        elif characteristic == "N Stage":
            current_group = "N Stage"

        elif characteristic == "M Stage":
            current_group = "M Stage"

        elif characteristic == "Vital Status":
            current_group = "Vital Status"


    # --------------------------------------------------------
    # AGE
    # --------------------------------------------------------

    if current_group == "Age":

        records.append({
            "Group": "Age",
            "Category": "Age",
            "Low": extract_age(row["Low Risk"]),
            "High": extract_age(row["High Risk"]),
            "P": row["P value"]
        })


    # --------------------------------------------------------
    # SEX
    # --------------------------------------------------------

    elif current_group == "Sex":

        records.append({
            "Group": "Sex",
            "Category": "Male",
            "Low": extract_percent(row["Low Risk"]),
            "High": extract_percent(row["High Risk"]),
            "P": row["P value"]
        })


    # --------------------------------------------------------
    # RACE
    # --------------------------------------------------------

    elif current_group == "Race":

        category = str(row["Characteristic"]).strip()

        if category.startswith("White"):
            category = "White"

        elif category.startswith("Black"):
            category = "Black"

        elif category.startswith("Asian"):
            category = "Asian"

        elif category.startswith("Race Not"):
            category = "Not reported"

        records.append({
            "Group": "Race",
            "Category": category,
            "Low": extract_percent(row["Low Risk"]),
            "High": extract_percent(row["High Risk"]),
            "P": row["P value"]
        })


    # --------------------------------------------------------
    # STAGE / VITAL STATUS
    # --------------------------------------------------------

    elif current_group in [
        "AJCC Stage",
        "T Stage",
        "N Stage",
        "M Stage",
        "Vital Status"
    ]:

        # Category is stored in Low Risk / High Risk
        low_text = str(row["Low Risk"])

        # Extract text before colon
        category_match = re.match(
            r"([^:]+):",
            low_text
        )

        if category_match:

            category = category_match.group(1).strip()

            records.append({
                "Group": current_group,
                "Category": category,
                "Low": extract_percent(row["Low Risk"]),
                "High": extract_percent(row["High Risk"]),
                "P": row["P value"]
            })


# ============================================================
# 7. CLEAN DATAFRAME
# ============================================================

clean = pd.DataFrame(records)

clean["P"] = pd.to_numeric(
    clean["P"],
    errors="coerce"
)

print("\n" + "=" * 70)
print("CLEANED DATA")
print("=" * 70)

print(clean.to_string(index=False))


# ============================================================
# 8. CATEGORY ORDER
# ============================================================

category_order = {

    "Sex": [
        "Male"
    ],

    "Race": [
        "White",
        "Black or African American",
        "Asian",
        "Race Not Reported"
    ],

    "AJCC Stage": [
        "Stage I",
        "Stage II",
        "Stage III",
        "Unknown"
    ],

    "T Stage": [
        "T1",
        "T2",
        "T3",
        "T4",
        "TX"
    ],

    "N Stage": [
        "N0",
        "N1",
        "N2",
        "N3",
        "NX"
    ],

    "M Stage": [
        "M0",
        "M1",
        "MX"
    ],

    "Vital Status": [
        "Alive",
        "Dead"
    ]
}


# ============================================================
# 9. CREATE FIGURE
# ============================================================

fig, axes = plt.subplots(
    4,
    2,
    figsize=(16, 20)
)

axes = axes.flatten()


# ============================================================
# 10. AGE PANEL
# ============================================================

ax = axes[0]

age_data = clean[
    clean["Group"] == "Age"
]

if not age_data.empty:

    low_age = age_data["Low"].iloc[0]
    high_age = age_data["High"].iloc[0]
    p_age = age_data["P"].iloc[0]

    ax.scatter(
        [0],
        [low_age],
        s=180,
        label="Low Risk"
    )

    ax.scatter(
        [1],
        [high_age],
        s=180,
        label="High Risk"
    )

    ax.plot(
        [0, 1],
        [low_age, high_age],
        linewidth=2
    )

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        ["Low Risk", "High Risk"]
    )

    ax.set_ylabel(
        "Median age at diagnosis (years)"
    )

    ax.set_title(
        "Age at Diagnosis",
        fontweight="bold"
    )

    ax.text(
        0.98,
        1.03,
        f"P = {p_age:.4f}",
        transform=ax.transAxes,
        ha="right",
        fontweight="bold"
    )

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3
    )


# ============================================================
# 11. CATEGORICAL PANELS
# ============================================================

panel_groups = [
    "Sex",
    "Race",
    "AJCC Stage",
    "T Stage",
    "N Stage",
    "M Stage",
    "Vital Status"
]


for ax, group in zip(
    axes[1:],
    panel_groups
):

    plot = clean[
        clean["Group"] == group
    ].copy()

    # --------------------------------------------------------
    # Apply category order
    # --------------------------------------------------------

    if group in category_order:

        order = category_order[group]

        plot["Category"] = pd.Categorical(
            plot["Category"],
            categories=order,
            ordered=True
        )

        plot = plot.sort_values(
            "Category"
        )

    # --------------------------------------------------------
    # Remove missing percentage rows
    # --------------------------------------------------------

    plot = plot.dropna(
        subset=["Low", "High"],
        how="all"
    )

    if plot.empty:

        ax.text(
            0.5,
            0.5,
            "No data available",
            ha="center",
            va="center"
        )

        ax.axis("off")

        continue


    y = np.arange(
        len(plot)
    )

    height = 0.34


    # --------------------------------------------------------
    # LOW RISK
    # --------------------------------------------------------

    ax.barh(
        y - height / 2,
        plot["Low"].fillna(0),
        height=height,
        label="Low Risk"
    )


    # --------------------------------------------------------
    # HIGH RISK
    # --------------------------------------------------------

    ax.barh(
        y + height / 2,
        plot["High"].fillna(0),
        height=height,
        label="High Risk"
    )


    # --------------------------------------------------------
    # DATA LABELS
    # --------------------------------------------------------

    for i, (_, r) in enumerate(
        plot.iterrows()
    ):

        if pd.notna(r["Low"]):

            ax.text(
                r["Low"] + 0.8,
                i - height / 2,
                f"{r['Low']:.1f}%",
                va="center",
                fontsize=8
            )

        if pd.notna(r["High"]):

            ax.text(
                r["High"] + 0.8,
                i + height / 2,
                f"{r['High']:.1f}%",
                va="center",
                fontsize=8
            )


    # --------------------------------------------------------
    # P-VALUE
    # --------------------------------------------------------

    p_values = plot["P"].dropna()

    if len(p_values) > 0:

        p = p_values.iloc[0]

        if p < 0.001:
            p_text = "P < 0.001"
        else:
            p_text = f"P = {p:.4f}"

        if p < 0.05:
            p_text += "  *"

        ax.text(
            0.98,
            1.03,
            p_text,
            transform=ax.transAxes,
            ha="right",
            fontweight="bold"
        )


    # --------------------------------------------------------
    # FORMAT
    # --------------------------------------------------------

    ax.set_yticks(y)

    ax.set_yticklabels(
        plot["Category"]
    )

    ax.set_xlabel(
        "Patients (%)"
    )

    ax.set_title(
        group,
        fontweight="bold"
    )

    ax.grid(
        axis="x",
        linestyle="--",
        alpha=0.3
    )

    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_xlim(
        0,
        100
    )


# ============================================================
# 12. LEGEND
# ============================================================

handles, labels = axes[1].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc="upper center",
    ncol=2,
    frameon=False,
    fontsize=11,
    bbox_to_anchor=(0.5, 0.985)
)


# ============================================================
# 13. TITLE
# ============================================================

fig.suptitle(
    "Clinicopathological Characteristics by Prognostic Risk Group",
    fontsize=20,
    fontweight="bold",
    y=0.998
)


# ============================================================
# 14. FOOTNOTE
# ============================================================

fig.text(
    0.5,
    0.005,
    "* P < 0.05. Percentages represent the proportion of patients "
    "within each risk group (n = 515 per group).",
    ha="center",
    fontsize=10,
    style="italic"
)


# ============================================================
# 15. LAYOUT
# ============================================================

plt.tight_layout(
    rect=[
        0,
        0.025,
        1,
        0.965
    ]
)


# ============================================================
# 16. SAVE
# ============================================================

output_dir = "results/clinical"

os.makedirs(
    output_dir,
    exist_ok=True
)

output_file = os.path.join(
    output_dir,
    "Clinicopathological_Characteristics_by_Risk_Group.png"
)

plt.savefig(
    output_file,
    dpi=600,
    bbox_inches="tight"
)

plt.show()

print("\n" + "=" * 70)
print("FIGURE SAVED")
print("=" * 70)

print(output_file)