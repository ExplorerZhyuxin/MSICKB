# -*- coding: utf-8 -*-
"""
09_M36_plot_pooled_mutation_forest.py

Plot pooled TCGA mutation support results for revised hub genes
as a forest-style odds ratio figure.

Input:
- pan_tcga_revised_hub_mutation_validation.csv

Outputs:
- m36_pooled_mutation_forest.png
- m36_pooled_mutation_forest.pdf
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =========================
# 0. Paths
# =========================
BASE_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M36_TCGA_validation"
input_file = os.path.join(BASE_DIR, "pan_tcga_revised_hub_mutation_validation.csv")

out_png = os.path.join(BASE_DIR, "m36_pooled_mutation_forest.png")
out_pdf = os.path.join(BASE_DIR, "m36_pooled_mutation_forest.pdf")

# =========================
# 1. Load data
# =========================
df = pd.read_csv(input_file)

required_core = ["Gene", "OR_haldane", "CI_lower", "CI_upper", "FDR", "MSIH_rate"]
missing_core = [c for c in required_core if c not in df.columns]
if missing_core:
    raise ValueError(f"Missing required columns: {missing_core}")

# Backward/forward compatibility for comparison-group rate column
if "Non_MSIH_rate" in df.columns:
    comp_rate_col = "Non_MSIH_rate"
    comp_label = "non-MSI-H"
elif "MSS_rate" in df.columns:
    comp_rate_col = "MSS_rate"
    comp_label = "MSS/MSI-L"
else:
    raise ValueError("Missing comparison-group rate column: expected 'Non_MSIH_rate' or 'MSS_rate'.")

# sort by OR for display
df = df.sort_values("OR_haldane", ascending=True).reset_index(drop=True)

# log scale safety
df = df[(df["OR_haldane"] > 0) & (df["CI_lower"] > 0) & (df["CI_upper"] > 0)].copy()

# direction for coloring
df["Direction"] = np.where(df["OR_haldane"] >= 1, "Higher in MSI-H", f"Lower in MSI-H")
color_map = {"Higher in MSI-H": "#d7301f", f"Lower in MSI-H": "#3182bd"}
colors = df["Direction"].map(color_map)

# y positions
y = np.arange(len(df))

# =========================
# 2. Plot
# =========================
fig, ax = plt.subplots(figsize=(10.2, 5.8), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# CI lines
for i, row in df.iterrows():
    ax.plot([row["CI_lower"], row["CI_upper"]], [i, i], color="gray", lw=1.2, zorder=1)

# points
ax.scatter(df["OR_haldane"], y, s=65, c=colors, edgecolor="black", linewidth=0.5, zorder=2)

# reference line
ax.axvline(1, color="black", linestyle="--", linewidth=1)

# y-axis
ax.set_yticks(y)
ax.set_yticklabels(df["Gene"], fontsize=11)

# x-axis
ax.set_xscale("log")
ax.set_xlabel(f"Odds ratio for mutation prevalence (MSI-H vs {comp_label})", fontsize=12)
ax.set_title("TCGA molecular support: mutation enrichment of revised hub genes", fontsize=14, pad=12)

# right-side annotation
x_text = df["CI_upper"].max() * 1.18
for i, row in df.iterrows():
    label = f'{row["MSIH_rate"]:.1f}% vs {row[comp_rate_col]:.1f}%   FDR={row["FDR"]:.3g}'
    ax.text(x_text, i, label, va="center", fontsize=9)

# x limits
xmin = max(df["CI_lower"].min() / 1.6, 0.1)
xmax = df["CI_upper"].max() * 3.0
ax.set_xlim(xmin, xmax)

# legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Higher in MSI-H',
           markerfacecolor=color_map["Higher in MSI-H"], markeredgecolor='black', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Lower in MSI-H',
           markerfacecolor=color_map[f"Lower in MSI-H"], markeredgecolor='black', markersize=8)
]
ax.legend(handles=legend_elements, frameon=False, loc="upper left", fontsize=10)

# style
ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.3)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(out_png, bbox_inches="tight", facecolor="white")
plt.savefig(out_pdf, bbox_inches="tight", facecolor="white")
plt.close()

print("Done.")
print(out_png)
print(out_pdf)
