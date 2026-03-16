# -*- coding: utf-8 -*-
"""
05_M34_publication_bias_sensitivity.py

Purpose:
PMID-based sensitivity analysis for literature-induced bias in the revised primary simple gene–cancer network.

Input:
- study_level_unique_gene_cancer_pmid.xlsx

Outputs:
- m34_gene_publication_bias_summary.xlsx
- m34_hub_vs_weighted_overlap.xlsx
- m34_pubfreq_vs_degree_scatter.png
- m34_pubfreq_vs_degree_scatter.pdf
- m34_publication_bias_summary.txt
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# =========================
# 0. Paths
# =========================
input_file = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\0_data\study_level_unique_gene_cancer_pmid.xlsx"
output_dir = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M34_publication_bias"
os.makedirs(output_dir, exist_ok=True)

# =========================
# 1. Load study-level data
# =========================
df = pd.read_excel(input_file)
df.columns = [str(c).strip() for c in df.columns]

required_cols = {"Source", "harmonized_cancer_name", "PMID", "include_in_primary_analysis"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# keep primary-analysis records only
df = df[df["include_in_primary_analysis"].astype(str).str.strip().str.lower() == "yes"].copy()

# keep essential columns
df = df[["Source", "harmonized_cancer_name", "PMID"]].copy()
df.columns = ["Gene", "Cancer", "PMID"]

# clean
df["Gene"] = df["Gene"].astype(str).str.strip()
df["Cancer"] = df["Cancer"].astype(str).str.strip()
df["PMID"] = df["PMID"].astype(str).str.strip()

# unique study-level records
df = df.drop_duplicates(subset=["Gene", "Cancer", "PMID"]).reset_index(drop=True)

# =========================
# 2. Build edge-level study counts
# =========================
edge_study = (
    df.groupby(["Gene", "Cancer"])["PMID"]
    .nunique()
    .reset_index(name="Edge_study_count")
)

# =========================
# 3. Build gene-level summary
# =========================
gene_degree = (
    edge_study.groupby("Gene")["Cancer"]
    .nunique()
    .reset_index(name="Simple_degree")
)

gene_pub = (
    df.groupby("Gene")["PMID"]
    .nunique()
    .reset_index(name="Unique_PMID_count")
)

gene_weighted = (
    edge_study.groupby("Gene")["Edge_study_count"]
    .sum()
    .reset_index(name="Weighted_degree_sum")
)

gene_edge_mean = (
    edge_study.groupby("Gene")["Edge_study_count"]
    .mean()
    .reset_index(name="Mean_edge_study_count")
)

summary = gene_degree.merge(gene_pub, on="Gene", how="outer") \
                     .merge(gene_weighted, on="Gene", how="outer") \
                     .merge(gene_edge_mean, on="Gene", how="outer")

summary = summary.fillna(0)

summary["Hub_degree_ge_3"] = summary["Simple_degree"] >= 3
summary["Normalized_degree_log"] = summary["Simple_degree"] / np.log2(summary["Unique_PMID_count"] + 1)
summary["Normalized_degree_sqrt"] = summary["Simple_degree"] / np.sqrt(summary["Unique_PMID_count"])

# ranks
summary["Rank_by_simple_degree"] = summary["Simple_degree"].rank(method="min", ascending=False).astype(int)
summary["Rank_by_weighted_degree"] = summary["Weighted_degree_sum"].rank(method="min", ascending=False).astype(int)
summary["Rank_by_normalized_degree_log"] = summary["Normalized_degree_log"].rank(method="min", ascending=False).astype(int)

summary = summary.sort_values(
    by=["Simple_degree", "Weighted_degree_sum", "Unique_PMID_count", "Gene"],
    ascending=[False, False, False, True]
).reset_index(drop=True)

# =========================
# 4. Compare top simple hubs vs weighted/normalized rankings
# =========================
hub_genes = summary.loc[summary["Hub_degree_ge_3"], "Gene"].tolist()
hub_genes_sorted = sorted(hub_genes)

n_hubs = len(hub_genes_sorted)

top_weighted = summary.sort_values(
    by=["Weighted_degree_sum", "Simple_degree", "Unique_PMID_count", "Gene"],
    ascending=[False, False, False, True]
).head(n_hubs)["Gene"].tolist()

top_normalized = summary.sort_values(
    by=["Normalized_degree_log", "Simple_degree", "Weighted_degree_sum", "Gene"],
    ascending=[False, False, False, True]
).head(n_hubs)["Gene"].tolist()

overlap_weighted = sorted(set(hub_genes_sorted).intersection(set(top_weighted)))
overlap_normalized = sorted(set(hub_genes_sorted).intersection(set(top_normalized)))

overlap_df = pd.DataFrame({
    "Comparison": [
        "Simple hubs (degree >= 3)",
        f"Top {n_hubs} by weighted_degree_sum",
        f"Top {n_hubs} by normalized_degree_log"
    ],
    "Genes": [
        ", ".join(hub_genes_sorted),
        ", ".join(top_weighted),
        ", ".join(top_normalized)
    ]
})

overlap_stats_df = pd.DataFrame({
    "Metric": [
        "Number of simple hubs",
        "Overlap: simple hubs vs top weighted",
        "Overlap: simple hubs vs top normalized"
    ],
    "Value": [
        n_hubs,
        len(overlap_weighted),
        len(overlap_normalized)
    ],
    "Genes": [
        ", ".join(hub_genes_sorted),
        ", ".join(overlap_weighted),
        ", ".join(overlap_normalized)
    ]
})

# =========================
# 5. Correlation analysis
# =========================
rho_pub, p_pub = spearmanr(summary["Simple_degree"], summary["Unique_PMID_count"])
rho_weighted, p_weighted = spearmanr(summary["Simple_degree"], summary["Weighted_degree_sum"])

# =========================
# 6. Save Excel outputs
# =========================
excel_file = os.path.join(output_dir, "m34_gene_publication_bias_summary.xlsx")
overlap_file = os.path.join(output_dir, "m34_hub_vs_weighted_overlap.xlsx")

with pd.ExcelWriter(excel_file, engine="openpyxl") as writer:
    summary.to_excel(writer, sheet_name="gene_level_summary", index=False)
    edge_study.to_excel(writer, sheet_name="edge_study_counts", index=False)

with pd.ExcelWriter(overlap_file, engine="openpyxl") as writer:
    overlap_df.to_excel(writer, sheet_name="gene_lists", index=False)
    overlap_stats_df.to_excel(writer, sheet_name="overlap_stats", index=False)

# =========================
# 7. Scatter plot: publication frequency vs degree
# =========================
plot_png = os.path.join(output_dir, "m34_pubfreq_vs_degree_scatter.png")
plot_pdf = os.path.join(output_dir, "m34_pubfreq_vs_degree_scatter.pdf")

plot_df = summary.copy()
hub_df = plot_df[plot_df["Hub_degree_ge_3"]].copy()
nonhub_df = plot_df[~plot_df["Hub_degree_ge_3"]].copy()

fig, ax = plt.subplots(figsize=(8.6, 6.2), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

ax.scatter(
    nonhub_df["Unique_PMID_count"],
    nonhub_df["Simple_degree"],
    s=42,
    color="#9ecae1",
    edgecolor="white",
    linewidth=0.5,
    alpha=0.85,
    label="Non-hub genes"
)

ax.scatter(
    hub_df["Unique_PMID_count"],
    hub_df["Simple_degree"],
    s=75,
    color="#e34a33",
    edgecolor="black",
    linewidth=0.6,
    alpha=0.95,
    label="Hub genes (degree ≥ 3)"
)

# label hub genes
for _, row in hub_df.iterrows():
    ax.text(
        row["Unique_PMID_count"] + 0.15,
        row["Simple_degree"] + 0.03,
        row["Gene"],
        fontsize=8.5
    )

ax.set_xlabel("Unique supporting PMIDs per gene", fontsize=12)
ax.set_ylabel("Simple degree (# cancer types)", fontsize=12)
ax.set_title("Publication frequency versus gene degree", fontsize=13, pad=12)

stats_text = f"Spearman ρ = {rho_pub:.2f}\np = {p_pub:.2e}"
ax.text(
    0.94, 0.20, stats_text,
    transform=ax.transAxes,
    ha="right", va="top",
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.9)
)

ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(frameon=False, fontsize=10, loc="lower right")

plt.tight_layout()
plt.savefig(plot_png, bbox_inches="tight", facecolor="white")
plt.savefig(plot_pdf, bbox_inches="tight", facecolor="white")
plt.close()

# =========================
# 8. Response-ready summary
# =========================
summary_txt = os.path.join(output_dir, "m34_publication_bias_summary.txt")

top_simple = summary.sort_values(
    by=["Simple_degree", "Weighted_degree_sum", "Unique_PMID_count", "Gene"],
    ascending=[False, False, False, True]
).head(15)

top_weighted_df = summary.sort_values(
    by=["Weighted_degree_sum", "Simple_degree", "Unique_PMID_count", "Gene"],
    ascending=[False, False, False, True]
).head(15)

top_norm_df = summary.sort_values(
    by=["Normalized_degree_log", "Simple_degree", "Weighted_degree_sum", "Gene"],
    ascending=[False, False, False, True]
).head(15)

text = f"""
Publication-bias sensitivity summary

Analysis scope:
- Input dataset: study_level_unique_gene_cancer_pmid.xlsx
- Restricted to include_in_primary_analysis = Yes
- Study-level unit: unique (gene, harmonized cancer, PMID)

Primary simple-network hubs (degree >= 3; n={n_hubs}):
- {", ".join(hub_genes_sorted)}

Correlation analyses:
- Spearman correlation between simple degree and unique PMID count:
  rho = {rho_pub:.3f}, p = {p_pub:.3g}
- Spearman correlation between simple degree and weighted degree sum:
  rho = {rho_weighted:.3f}, p = {p_weighted:.3g}

Top {n_hubs} genes by publication-weighted support:
- {", ".join(top_weighted)}

Top {n_hubs} genes by publication-normalized degree:
- {", ".join(top_normalized)}

Overlap with simple-network hubs:
- Simple hubs vs top weighted: {len(overlap_weighted)}/{n_hubs}
  Genes: {", ".join(overlap_weighted)}
- Simple hubs vs top normalized: {len(overlap_normalized)}/{n_hubs}
  Genes: {", ".join(overlap_normalized)}

Interpretation:
- The correlation analyses quantify the extent to which publication frequency is associated with gene degree.
- The weighted and normalized comparisons provide a sensitivity check for literature-induced bias rather than a replacement for the primary simple-network definition.

Top genes by simple degree:
{top_simple[['Gene','Simple_degree','Unique_PMID_count','Weighted_degree_sum','Normalized_degree_log']].to_string(index=False)}

Top genes by weighted degree:
{top_weighted_df[['Gene','Simple_degree','Unique_PMID_count','Weighted_degree_sum','Normalized_degree_log']].to_string(index=False)}

Top genes by normalized degree:
{top_norm_df[['Gene','Simple_degree','Unique_PMID_count','Weighted_degree_sum','Normalized_degree_log']].to_string(index=False)}
""".strip()

with open(summary_txt, "w", encoding="utf-8") as f:
    f.write(text)

print("Done.")
print(f"Files saved to: {output_dir}")
print(f"- {excel_file}")
print(f"- {overlap_file}")
print(f"- {plot_png}")
print(f"- {plot_pdf}")
print(f"- {summary_txt}")
