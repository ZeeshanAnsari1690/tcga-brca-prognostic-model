# import libraries

import pandas as pd
import numpy as np

from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu

# load clean datasets

df = pd.read_csv(
    r"C:\Cancer\TCGA-BRCA\results\clinical_characteristics_dataset_clean.csv"
)

print("="*60)
print("Clinical Characteristics Dataset")
print("="*60)
print(df.shape)

for col in df.columns:
    print(col)

# creating low/high risk group


high = df[df["Risk_Group"]=="High"].copy()
low  = df[df["Risk_Group"]=="Low"].copy()

print("\nHigh Risk :",len(high))
print("Low Risk  :",len(low))

# Helper Functions  for analysis


# 1.percentage
def percentage(n,total):
    return f"{n} ({100*n/total:.1f}%)"

# 2.chia/square test
def categorical_pvalue(data,column):

    temp=data[[column,"Risk_Group"]].dropna()

    table=pd.crosstab(
        temp[column],
        temp["Risk_Group"]
    )

    if table.shape==(2,2):

        if (table.values<5).any():

            _,p=fisher_exact(table)

        else:

            _,p,_,_=chi2_contingency(table)

    else:

        _,p,_,_=chi2_contingency(table)

    return p

# 3.Age stats
def age_statistics():

    age=df["diagnoses.age_at_diagnosis"]

    overall=f"{age.median():.1f} ({age.quantile(0.25):.1f}-{age.quantile(0.75):.1f})"

    age_low=low["diagnoses.age_at_diagnosis"]

    low_text=f"{age_low.median():.1f} ({age_low.quantile(0.25):.1f}-{age_low.quantile(0.75):.1f})"

    age_high=high["diagnoses.age_at_diagnosis"]

    high_text=f"{age_high.median():.1f} ({age_high.quantile(0.25):.1f}-{age_high.quantile(0.75):.1f})"

    _,p=mannwhitneyu(
        age_low,
        age_high,
        alternative="two-sided"
    )

    return overall,low_text,high_text,p

# creating tabel
table=[]


# Age Row
overall,low_age,high_age,p=age_statistics()

table.append({

"Characteristic":"Age (years), median (IQR)",

"Overall":overall,

"Low Risk":low_age,

"High Risk":high_age,

"P value":round(p,4)

})

# Function for Categorical Variables

def add_variable(column, title):

    p = categorical_pvalue(df, column)

    categories = sorted(df[column].dropna().unique())

    # Detect binary variables (0/1 or False/True)
    if set(categories).issubset({0, 1, False, True}):

        positive = 1 if 1 in categories else True

        overall = df[column].eq(positive).sum()
        low_n = low[column].eq(positive).sum()
        high_n = high[column].eq(positive).sum()

        table.append({
            "Characteristic": title,
            "Overall": percentage(overall, len(df)),
            "Low Risk": percentage(low_n, len(low)),
            "High Risk": percentage(high_n, len(high)),
            "P value": round(p, 4)
        })

    else:
        for i, cat in enumerate(categories):

            overall = df[column].eq(cat).sum()
            low_n = low[column].eq(cat).sum()
            high_n = high[column].eq(cat).sum()

            table.append({

                "Characteristic": title if i == 0 else "",
                "Overall": f"{cat}: {percentage(overall, len(df))}",
                "Low Risk": percentage(low_n, len(low)),
                "High Risk": percentage(high_n, len(high)),
                "P value": round(p, 4) if i == 0 else ""

            })


        # run analysis--------------------------

def simplify_ajcc_stage(stage):

    if pd.isna(stage):
        return np.nan

    stage = str(stage).strip().upper()

    if stage.startswith("I") and not stage.startswith("II"):
        return "Stage I"

    elif stage.startswith("II") and not stage.startswith("III"):
        return "Stage II"

    elif stage.startswith("III"):
        return "Stage III"

    elif stage.startswith("IV"):
        return "Stage IV"

    elif stage in ["X", "STAGE X"]:
        return "Unknown"

    else:
        return "Unknown"

                        # creating new col for better stages representation

df["AJCC_Stage_Group"] = df["diagnoses.ajcc_pathologic_stage"].apply(simplify_ajcc_stage)

low["AJCC_Stage_Group"] = low["diagnoses.ajcc_pathologic_stage"].apply(simplify_ajcc_stage)

high["AJCC_Stage_Group"] = high["diagnoses.ajcc_pathologic_stage"].apply(simplify_ajcc_stage)






        # Sex
# Sex
add_variable(
    "demographic.sex_at_birth_male",
    "Male, n (%)"
)

# Race
add_variable(
    "demographic.race_white",
    "White, n (%)"
)

add_variable(
    "demographic.race_black or african american",
    "Black or African American, n (%)"
)

add_variable(
    "demographic.race_asian",
    "Asian, n (%)"
)

add_variable(
    "demographic.race_not reported",
    "Race Not Reported, n (%)"
)

# Clinical variables
add_variable(
    "AJCC_Stage_Group",
    "AJCC Stage"
)

add_variable(
    "diagnoses.ajcc_pathologic_t",
    "T Stage"
)

add_variable(
    "diagnoses.ajcc_pathologic_n",
    "N Stage"
)

add_variable(
    "diagnoses.ajcc_pathologic_m",
    "M Stage"
)

add_variable(
    "demographic.vital_status",
    "Vital Status"
)
# Export Results

results=pd.DataFrame(table)

csv_file=r"C:\Cancer\TCGA-BRCA\results\Clinical_Characteristics_Table.csv"

xlsx_file=r"C:\Cancer\TCGA-BRCA\results\Clinical_Characteristics_Table.xlsx"

results.to_csv(csv_file,index=False)

results.to_excel(xlsx_file,index=False)

print("="*60)
print("Clinical Characteristics Table Generated")
print("="*60)
print(results.head(20))

print("\nSaved:")
print(csv_file)
print(xlsx_file)
