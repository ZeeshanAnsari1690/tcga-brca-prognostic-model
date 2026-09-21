# ==========================================================
# TCGA-BRCA Data Visualization
# ==========================================================
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Improve figure quality
plt.style.use("ggplot")
sns.set(font_scale=1.1)

# Load QC dataset
df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\final_tcga_brca_dataset_QC.csv"
)

# print("Dataset Shape:", df.shape)

# # # Event distribution

# plt.figure(figsize=(6,5))

# sns.countplot(
#     data=df,
#     x="event"
# )

# plt.title("Distribution of Survival Events")
# plt.xlabel("Event")
# plt.ylabel("Number of Patients")

# plt.xticks([0,1], ["Alive (0)", "Dead (1)"])

# plt.tight_layout()

# plt.savefig(
#     r"C:\Cancer\TCGA-BRCA\figures\event_distribution.png",
#     dpi=300
# )

# plt.show()

# # # Survival Time Distribution
# # # ==========================================================

# plt.figure(figsize=(8,5))

# sns.histplot(
#     df["survival_time"],
#     bins=30,
#     kde=True
# )

# plt.title("Survival Time Distribution")
# plt.xlabel("Survival Time (Days)")
# plt.ylabel("Frequency")

# figure_path = r"C:\Cancer\TCGA-BRCA\figures"

# os.makedirs(figure_path, exist_ok=True)

# plt.tight_layout()
# plt.savefig(os.path.join(figure_path,"02_survival_time_distribution.png"))
# plt.close()

# Sex Distribution

plt.figure(figsize=(6,5))

sns.countplot(data=df, x="demographic.sex_at_birth")

plt.title("Sex Distribution")

figure_path = r"C:\Cancer\TCGA-BRCA\figures"

os.makedirs(figure_path, exist_ok=True)

plt.tight_layout()
plt.savefig(os.path.join(figure_path,"03_sex_distribution.png"))
plt.close()

# 
plt.figure(figsize=(8,5))

sns.countplot(
    data=df,
    x="demographic.race",
    order=df["demographic.race"].value_counts().index
)

plt.xticks(rotation=45)

plt.title("Race Distribution")

figure_path = r"C:\Cancer\TCGA-BRCA\figures"

os.makedirs(figure_path, exist_ok=True)


plt.tight_layout()
plt.savefig(os.path.join(figure_path,"03_race_distribution.png"))
plt.close()

# Age distribution


fig, ax = plt.subplots(1,2,figsize=(12,5))

sns.histplot(df["diagnoses.age_at_diagnosis"], bins=25, kde=True, ax=ax[0])
ax[0].set_title("Age Distribution")

sns.boxplot(y=df["diagnoses.age_at_diagnosis"], ax=ax[1])
ax[1].set_title("Age Boxplot")

figure_path = r"C:\Cancer\TCGA-BRCA\figures"

os.makedirs(figure_path, exist_ok=True)

plt.tight_layout()
plt.savefig(os.path.join(figure_path,"03_age_distribution.png"))
plt.close()

# Age in years
# ==========================================================
# Task 3: Age Distribution (Years)
# ==========================================================

# Convert age from days to years
df["age_years"] = df["diagnoses.age_at_diagnosis"] / 365.25

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Histogram
sns.histplot(
    df["age_years"].dropna(),
    bins=25,
    kde=True,
    ax=ax[0]
)

ax[0].set_title("Age Distribution")
ax[0].set_xlabel("Age at Diagnosis (Years)")
ax[0].set_ylabel("Number of Patients")

# Boxplot
sns.boxplot(
    y=df["age_years"].dropna(),
    ax=ax[1]
)

ax[1].set_title("Age Boxplot")
ax[1].set_ylabel("Age at Diagnosis (Years)")

figure_path = r"C:\Cancer\TCGA-BRCA\figures"

os.makedirs(figure_path, exist_ok=True)



plt.tight_layout()
plt.savefig(
    os.path.join(figure_path, "03_age_distribution.png"),
    dpi=300
)
plt.show()
plt.close()


# RNA-seq preprocessing for machine learning
