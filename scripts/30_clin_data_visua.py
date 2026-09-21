import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------------------------------
# Clinical Heatmap
# -------------------------------------------------------

# Create heatmap by using data of clinical characteristics
import pandas as pd

heatmap = pd.DataFrame(
    {
        "Low Risk":[59.2,1.7,64.1,17.7,7.6,17.5,10.7,7.8,11.7],
        "High Risk":[58.8,0.4,74.8,15.1,4.1,24.5,12.2,13.6,17.5],
        "P value":[0.0864,0.0640,0.0003,0.3127,0.0237,0.0301,0.8774,0.0091,0.0104]
    },
    index=[
        "Age (years)",
        "Male (%)",
        "White (%)",
        "Black (%)",
        "Asian (%)",
        "Stage III (%)",
        "T3 (%)",
        "N2 (%)",
        "Dead (%)"
    ]
)

# -------------------------------------------------------
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Figure
# -----------------------------
fig = plt.figure(figsize=(9,8))

# Heatmap area
ax = plt.axes([0.22,0.12,0.48,0.75])

# P-value area
ax2 = plt.axes([0.73,0.12,0.18,0.75])

# -----------------------------
# Row-wise normalization
# -----------------------------
# -----------------------------
# Heatmap values
# -----------------------------
values = heatmap[["Low Risk","High Risk"]].values

# -----------------------------
# Heatmap
# -----------------------------
im = ax.imshow(
    values,
    cmap="coolwarm",
    aspect="auto"
)

# -----------------------------
# Cell values
# -----------------------------
for i in range(values.shape[0]):
    for j in range(2):

        ax.text(
            j,
            i,
            f"{values[i,j]:.1f}",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color="black"
        )

# -----------------------------
# Axis labels
# -----------------------------
ax.set_xticks([0,1])
ax.set_xticklabels(
    ["Low Risk","High Risk"],
    fontsize=11,
    fontweight="bold"
)

ax.set_yticks(range(len(heatmap.index)))
ax.set_yticklabels(
    heatmap.index,
    fontsize=12,
    fontweight="bold"
)

# -----------------------------
# White borders
# -----------------------------
ax.set_xticks([-0.5,0.5,1.5],minor=True)
ax.set_yticks(np.arange(-0.5,len(heatmap.index),1),minor=True)

ax.grid(
    which="minor",
    color="white",
    linewidth=2
)

ax.tick_params(which="minor",bottom=False,left=False)

# -----------------------------
# Remove spines
# -----------------------------
for s in ax.spines.values():
    s.set_visible(False)

# -----------------------------
# P-value column
# -----------------------------
# -----------------------------
# P-value column
# -----------------------------
ax2.set_xlim(0,1)
ax2.set_ylim(-0.5, len(heatmap)-0.5)
ax2.invert_yaxis()
ax2.axis("off")

# Header
ax2.text(
    0.5,
    -0.7,
    "P-value",
    fontsize=14,
    fontweight="bold",
    ha="center"
)

# Values
for i, p in enumerate(heatmap["P value"]):

    ax2.text(
        0.5,
        i,
        f"{p:.4f}",
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold" if p < 0.05 else "normal",
        color="#B22222" if p < 0.05 else "black"
    )
# -----------------------------
# Vertical separator
# -----------------------------
fig.lines.append(
    plt.Line2D(
        [0.71,0.71],
        [0.12,0.87],
        transform=fig.transFigure,
        color="#777777",
        linewidth=1.5
    )
)

# -----------------------------
# Title
# -----------------------------
plt.suptitle(
    " Clinicopathological Characteristics by Risk Group",
    fontsize=17,
    fontweight="bold",
    y=0.96
)

# -----------------------------
# Footnote
# -----------------------------
plt.figtext(
    0.22,
    0.04,
    "Bold red P-values indicate statistically significant differences between Low-risk and High-risk groups (P < 0.05).",
    fontsize=10,
    style="italic"
)

# -----------------------------
# Save
# -----------------------------
plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Clinical_Heatmap_Publication.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    r"C:\Cancer\TCGA-BRCA\results\Clinical_Heatmap_Publication.pdf",
    bbox_inches="tight"
)

plt.show()