# loading signi lib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

deg = pd.read_csv("results/DEG_results_annotated.csv")

# creating ploting catogory
deg["Category"] = "Not Significant"

deg.loc[
    (deg["FDR"] < 0.05) &
    (deg["log2FC"] > 1),
    "Category"
] = "Up"

deg.loc[
    (deg["FDR"] < 0.05) &
    (deg["log2FC"] < -1),
    "Category"
] = "Down"

# Calculate −log10(FDR)
deg["minusLog10FDR"] = -np.log10(deg["FDR"] + 1e-300)

# Creating  the volcano plot
plt.figure(figsize=(10,8))

colors = {
    "Not Significant":"lightgray",
    "Up":"red",
    "Down":"blue"
}

for category, color in colors.items():

    subset = deg[deg["Category"] == category]

    plt.scatter(
        subset["log2FC"],
        subset["minusLog10FDR"],
        c=color,
        s=12,
        alpha=0.7,
        label=category
    )

# Threshold lines
plt.axvline(1, color="black", linestyle="--", linewidth=1)
plt.axvline(-1, color="black", linestyle="--", linewidth=1)
plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)

plt.xlabel("log2 Fold Change", fontsize=12)
plt.ylabel("-log10(FDR)", fontsize=12)
plt.title("Volcano Plot of Differentially Expressed Genes", fontsize=14)

plt.legend()

plt.tight_layout()

plt.savefig(
    "results/Volcano_plot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# label top gene
top_genes = deg.nsmallest(10, "FDR")

for _, row in top_genes.iterrows():
    plt.text(
        row["log2FC"],
        row["minusLog10FDR"],
        row["gene_name"],
        fontsize=8
    )
# Heatmap of Top 50 DEGs

#  imp lib
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# load data

expr = pd.read_csv("results/expression_with_risk_score.csv")
deg = pd.read_csv("results/Significant_Protein_Coding_DEGs.csv")

# selecting top 50 DEGs
top50 = (
    deg.assign(absLog2FC=deg["log2FC"].abs())
       .sort_values(["FDR", "absLog2FC"],
                    ascending=[True, False])
       .head(50)
)

# ext expression matrix
top50_genes = top50["gene_id"].tolist()

gene_matrix = expr[top50_genes].copy()

# check avail gene

available_genes = [g for g in top50_genes if g in expr.columns]

print("Available genes:", len(available_genes))

gene_matrix = expr[available_genes].copy()



# rename gene by using symbols
mapping = dict(zip(top50["gene_id"], top50["gene_name"]))

gene_matrix.columns = [
    mapping.get(gene, gene)
    for gene in gene_matrix.columns
]

# add risks groups

median_score = expr["Risk_Score"].median()

expr["Risk_Group"] = np.where(
    expr["Risk_Score"] > median_score,
    "High",
    "Low"
)

# Add Risk_Group to gene matrix
gene_matrix["Risk_Group"] = expr["Risk_Group"]

# patients sorting 


# Add Risk Score
gene_matrix["Risk_Score"] = expr["Risk_Score"]

# Sort by increasing Risk Score
gene_matrix = gene_matrix.sort_values("Risk_Score")

# Remove Risk Score before plotting
gene_matrix = gene_matrix.drop(columns="Risk_Score")

# prep heat matrix
heatmap_data = gene_matrix.drop(columns="Risk_Group")

heatmap_data = heatmap_data.T

# noramalization
heatmap_scaled = heatmap_data.apply(
    lambda x: (x - x.mean()) / (x.std() + 1e-8),
    axis=1
)

print("Heatmap shape:", heatmap_scaled.shape)

print("\nTotal NaN values:")
print(heatmap_scaled.isna().sum().sum())

print("\nNaN values per gene:")
print(heatmap_scaled.isna().sum(axis=1))

print("\nAny row completely NaN?")
print(heatmap_scaled.isna().all(axis=1).sum())

print(heatmap_scaled.dtypes.unique())



# ---------------------------------------
#  Heatmap


# Add a High-/Low-risk annotation bar
risk_colors = gene_matrix["Risk_Group"].map({
    "Low": "#3B82F6",
    "High": "#DC2626"
})

g = sns.clustermap(
    heatmap_scaled,
    row_cluster=False,
    col_cluster=False,
    cmap="RdBu_r",
    center=0,
    vmin=-2,
    vmax=2,
    figsize=(14, 9),
    dendrogram_ratio=(0.02, 0.02),
    colors_ratio=0.02,
    xticklabels=False,
    yticklabels=True,
    cbar_pos=(0.98, 0.25, 0.02, 0.40) ,   # Right side
    cbar_kws={"label": "Scaled Expression (Z-score)"}

)

# add legend

from matplotlib.patches import Patch

handles = [
    Patch(facecolor="#3B82F6", label="Low Risk"),
    Patch(facecolor="#DC2626", label="High Risk")
]

g.ax_heatmap.legend(
    handles=handles,
    title="Risk Group",
    bbox_to_anchor=(1.18,1.02),
    loc="upper left",
    frameon=False
)

g.ax_heatmap.tick_params(
    axis="y",
    labelsize=8
)

for label in g.ax_heatmap.get_yticklabels():
    label.set_fontstyle("italic")
    label.set_fontsize(9)


g.fig.suptitle(
    "Top 50 Differentially Expressed Genes",
    fontsize=18,
    fontweight="bold",
    y=1.02
)

g.savefig(
    "results/Heatmap_Top50_DEGs.png",
    dpi=600,
    bbox_inches="tight"
)

