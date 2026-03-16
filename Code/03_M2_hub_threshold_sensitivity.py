# -*- coding: utf-8 -*-
"""
03_M2_hub_threshold_sensitivity.py

Purpose:
1. Calculate gene degree in the primary simple gene–cancer network
2. Evaluate alternative hub thresholds: k >= 2, 3, 4, 5
3. Summarize hub counts, proportions, and gene lists
4. Test association between hub status and universal-gene status
5. Export concise tables and one sensitivity plot for M2 response

Input:
- primary_simple_gene_cancer_edges.xlsx
  Required columns:
    - Gene
    - Cancer

Output:
- m2_threshold_summary.xlsx
- m2_hub_gene_lists.xlsx
- m2_hub_universal_association.xlsx
- m2_threshold_sensitivity_plot.png
- m2_response_ready_summary.txt
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

# =========================
# 0. Paths
# =========================
input_file = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\0_data\primary_simple_gene_cancer_edges.xlsx"
output_dir = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M2"
os.makedirs(output_dir, exist_ok=True)

# =========================
# 1. Load data
# =========================
df = pd.read_excel(input_file)
df.columns = [str(c).strip() for c in df.columns]

required_cols = {"Gene", "Cancer"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# Keep only needed columns and rename for internal consistency
df = df[["Gene", "Cancer"]].copy()
df.columns = ["Gene", "Cancer_Type"]

# Ensure unique simple edges
df = df.drop_duplicates()

# =========================
# 2. Calculate gene degree
# =========================
degree_df = (
    df.groupby("Gene")["Cancer_Type"]
    .nunique()
    .reset_index(name="Degree")
    .sort_values(["Degree", "Gene"], ascending=[False, True])
    .reset_index(drop=True)
)

total_genes = degree_df["Gene"].nunique()

# Universal gene: connected to >=2 cancer types
degree_df["Universal"] = degree_df["Degree"] >= 2

# Descriptive statistics
degree_median = degree_df["Degree"].median()
degree_mean = degree_df["Degree"].mean()
degree_max = degree_df["Degree"].max()

# Degree frequency table
degree_freq_df = (
    degree_df["Degree"]
    .value_counts()
    .sort_index()
    .reset_index()
)
degree_freq_df.columns = ["Degree", "Gene_count"]

# =========================
# 3. Threshold sensitivity
# =========================
thresholds = [2, 3, 4, 5]

summary_rows = []
hub_list_frames = []
assoc_rows = []

for k in thresholds:
    temp = degree_df.copy()
    temp["Hub"] = temp["Degree"] >= k

    hub_df = temp[temp["Hub"]].copy().sort_values(["Degree", "Gene"], ascending=[False, True])
    nonhub_df = temp[~temp["Hub"]].copy()

    hub_count = len(hub_df)
    hub_prop = hub_count / total_genes if total_genes > 0 else np.nan

    # 2x2 contingency table
    #                Universal   Not Universal
    # Hub
    # Non-hub
    a = int(((temp["Hub"] == True) & (temp["Universal"] == True)).sum())
    b = int(((temp["Hub"] == True) & (temp["Universal"] == False)).sum())
    c = int(((temp["Hub"] == False) & (temp["Universal"] == True)).sum())
    d = int(((temp["Hub"] == False) & (temp["Universal"] == False)).sum())

    contingency = [[a, b], [c, d]]

    try:
        odds_ratio, p_value = fisher_exact(contingency, alternative="two-sided")
    except Exception:
        odds_ratio, p_value = np.nan, np.nan

    if np.isinf(odds_ratio):
        odds_ratio_str = "Inf"
    elif pd.notna(odds_ratio):
        odds_ratio_str = f"{odds_ratio:.4f}"
    else:
        odds_ratio_str = "NA"

    hub_universal_count = int(hub_df["Universal"].sum())
    hub_universal_prop = hub_universal_count / hub_count if hub_count > 0 else np.nan
    nonhub_universal_count = int(nonhub_df["Universal"].sum())
    nonhub_count = len(nonhub_df)
    nonhub_universal_prop = nonhub_universal_count / nonhub_count if nonhub_count > 0 else np.nan

    hub_gene_list = "; ".join(hub_df["Gene"].tolist()) if hub_count > 0 else ""

    summary_rows.append({
        "Threshold": f"Degree >= {k}",
        "k": k,
        "Hub_count": hub_count,
        "Hub_proportion": round(hub_prop, 4),
        "Hub_percentage": round(hub_prop * 100, 2),
        "Interpretation": f"{hub_count}/{total_genes} genes ({hub_prop*100:.2f}%)",
        "Hub_genes": hub_gene_list
    })

    assoc_rows.append({
        "Threshold": f"Degree >= {k}",
        "k": k,
        "Hub_count": hub_count,
        "Hub_universal_count": hub_universal_count,
        "Hub_universal_proportion": round(hub_universal_prop, 4) if pd.notna(hub_universal_prop) else np.nan,
        "Nonhub_count": nonhub_count,
        "Nonhub_universal_count": nonhub_universal_count,
        "Nonhub_universal_proportion": round(nonhub_universal_prop, 4) if pd.notna(nonhub_universal_prop) else np.nan,
        "Fisher_table_[HubUniversal,HubNonUniversal;NonhubUniversal,NonhubNonUniversal]": str(contingency),
        "Odds_ratio": odds_ratio_str,
        "P_value": p_value
    })

    if hub_count > 0:
        temp_hub_list = hub_df.copy()
        temp_hub_list.insert(0, "Threshold", f"Degree >= {k}")
        hub_list_frames.append(temp_hub_list)

summary_df = pd.DataFrame(summary_rows)
assoc_df = pd.DataFrame(assoc_rows)
hub_lists_df = pd.concat(hub_list_frames, ignore_index=True) if hub_list_frames else pd.DataFrame()

# =========================
# 4. Overlap / retention summary
# =========================
threshold_to_genes = {
    k: set(degree_df.loc[degree_df["Degree"] >= k, "Gene"].tolist())
    for k in thresholds
}

overlap_rows = []
pairs = [(2, 3), (3, 4), (4, 5), (3, 5)]
for k1, k2 in pairs:
    s1 = threshold_to_genes[k1]
    s2 = threshold_to_genes[k2]
    shared = sorted(list(s1 & s2))
    retention = len(shared) / len(s1) if len(s1) > 0 else np.nan
    overlap_rows.append({
        "From_threshold": f"Degree >= {k1}",
        "To_threshold": f"Degree >= {k2}",
        "From_count": len(s1),
        "To_count": len(s2),
        "Shared_count": len(shared),
        "Retention_from_first_threshold": round(retention, 4) if pd.notna(retention) else np.nan,
        "Shared_genes": "; ".join(shared)
    })

overlap_df = pd.DataFrame(overlap_rows)

# Supplement-ready combined table
supp_table_df = summary_df[[
    "Threshold", "k", "Hub_count", "Hub_proportion", "Hub_percentage", "Interpretation", "Hub_genes"
]].sort_values("k")


# =========================
# 5. Save Excel outputs
# =========================
summary_file = os.path.join(output_dir, "m2_threshold_summary.xlsx")
hub_list_file = os.path.join(output_dir, "m2_hub_gene_lists.xlsx")
assoc_file = os.path.join(output_dir, "m2_hub_universal_association.xlsx")
degree_freq_file = os.path.join(output_dir, "m2_degree_frequency.xlsx")

with pd.ExcelWriter(summary_file, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="threshold_summary", index=False)
    overlap_df.to_excel(writer, sheet_name="overlap_summary", index=False)
    supp_table_df.to_excel(writer, sheet_name="supplement_ready_table", index=False)

with pd.ExcelWriter(hub_list_file, engine="openpyxl") as writer:
    hub_lists_df.to_excel(writer, sheet_name="hub_gene_lists", index=False)
    degree_df.to_excel(writer, sheet_name="all_gene_degrees", index=False)

assoc_df.to_excel(assoc_file, index=False)
degree_freq_df.to_excel(degree_freq_file, index=False)

# =========================
# 6. Plot: threshold sensitivity
# =========================
plot_file = os.path.join(output_dir, "m2_threshold_sensitivity_plot.png")
plot_file_pdf = os.path.join(output_dir, "m2_threshold_sensitivity_plot.pdf")

plot_df = summary_df.sort_values("k").copy()

plt.style.use("default")
fig, ax1 = plt.subplots(figsize=(7.2, 4.8), dpi=300)
fig.patch.set_facecolor("white")
ax1.set_facecolor("white")

x = np.arange(len(plot_df))
bar_color = "#4C78A8"
line_color = "#D94F45"

bars = ax1.bar(
    x,
    plot_df["Hub_count"],
    color=bar_color,
    alpha=0.88,
    width=0.62,
    edgecolor="black",
    linewidth=0.6
)

ax1.set_xlabel("Degree threshold", fontsize=12)
ax1.set_ylabel("Number of hub genes", fontsize=12, color=bar_color)
ax1.set_xticks(x)
ax1.set_xticklabels([f"≥{k}" for k in plot_df["k"]], fontsize=11)
ax1.tick_params(axis="y", labelcolor=bar_color, labelsize=10)
ax1.tick_params(axis="x", labelsize=11)

# Clean spines
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Mild grid for readability
ax1.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.3)
ax1.set_axisbelow(True)

# Add count labels
for rect, value in zip(bars, plot_df["Hub_count"]):
    ax1.text(
        rect.get_x() + rect.get_width() / 2,
        rect.get_height() + 0.35,
        str(value),
        ha="center",
        va="bottom",
        fontsize=10,
        color="black"
    )

# Secondary axis for percentages
ax2 = ax1.twinx()
ax2.plot(
    x,
    plot_df["Hub_percentage"],
    color=line_color,
    marker="o",
    markersize=6.5,
    linewidth=2.2
)
ax2.set_ylabel("Percentage of all genes (%)", fontsize=12, color=line_color)
ax2.tick_params(axis="y", labelcolor=line_color, labelsize=10)
ax2.spines["top"].set_visible(False)

# Set upper limit with headroom
ax1.set_ylim(0, max(plot_df["Hub_count"]) + 3)
ax2.set_ylim(0, max(plot_df["Hub_percentage"]) + 3)

# Add percentage labels
for xi, value in zip(x, plot_df["Hub_percentage"]):
    ax2.text(
        xi,
        value + 2.0,
        f"{value:.1f}%",
        color=line_color,
        ha="center",
        va="bottom",
        fontsize=9.5
    )

plt.title("Sensitivity of hub-gene counts to alternative degree thresholds", fontsize=13, pad=12)

fig.tight_layout()
plt.savefig(plot_file, bbox_inches="tight", facecolor="white")
plt.savefig(plot_file_pdf, bbox_inches="tight", facecolor="white")
plt.close()


# =========================
# 7. Response-ready text
# =========================
txt_file = os.path.join(output_dir, "m2_response_ready_summary.txt")

k2_count = int(summary_df.loc[summary_df["k"] == 2, "Hub_count"].values[0])
k3_count = int(summary_df.loc[summary_df["k"] == 3, "Hub_count"].values[0])
k4_count = int(summary_df.loc[summary_df["k"] == 4, "Hub_count"].values[0])
k5_count = int(summary_df.loc[summary_df["k"] == 5, "Hub_count"].values[0])

k3_pct = float(summary_df.loc[summary_df["k"] == 3, "Hub_percentage"].values[0])

k3_assoc = assoc_df.loc[assoc_df["k"] == 3].iloc[0]


k3_genes = sorted(list(threshold_to_genes[3]))
k4_genes = sorted(list(threshold_to_genes[4]))
k5_genes = sorted(list(threshold_to_genes[5]))

response_text = f"""
M2 response-ready summary

In the revised primary simple gene–cancer network, gene degree was highly right-skewed
(median = {degree_median:.1f}, mean = {degree_mean:.2f}, max = {degree_max}).

Sensitivity analysis across alternative degree thresholds showed:
- degree >= 2: {k2_count} genes
- degree >= 3: {k3_count} genes
- degree >= 4: {k4_count} genes
- degree >= 5: {k5_count} genes

Thus, genes with degree >= 3 accounted for {k3_count}/{total_genes} genes ({k3_pct:.2f}%),
placing this cutoff within the sparse upper tail of the empirical degree distribution.

Hub genes at degree >= 3:
{", ".join(k3_genes)}

Hub genes at degree >= 4:
{", ".join(k4_genes)}

Hub genes at degree >= 5:
{", ".join(k5_genes)}

Hub-versus-universal association:
- degree >= 3: OR = {k3_assoc['Odds_ratio']}, p = {k3_assoc['P_value']:.4g}

Suggested interpretation:
We adopted degree >= 3 as an operational hub threshold in the revised manuscript because
it captures a small upper-tail subset of highly connected genes while retaining sufficient
interpretability for downstream functional analysis. Stricter thresholds (degree >= 4 and
degree >= 5) yielded smaller but overlapping core hub sets, indicating that the main
interpretation of high-connectivity genes is not driven by a single arbitrary cutoff.
""".strip()

with open(txt_file, "w", encoding="utf-8") as f:
    f.write(response_text)

# =========================
# 8. Console summary
# =========================
print("Done.")
print(f"Total genes: {total_genes}")
print(f"Degree summary: median={degree_median:.1f}, mean={degree_mean:.2f}, max={degree_max}")
print("Threshold counts:")
print(summary_df[["Threshold", "Hub_count", "Hub_percentage"]].to_string(index=False))
print(f"\nFiles saved to: {output_dir}")
print(f"- {summary_file}")
print(f"- {hub_list_file}")
print(f"- {assoc_file}")
print(f"- {degree_freq_file}")
print(f"- {plot_file}")
print(f"- {txt_file}")
