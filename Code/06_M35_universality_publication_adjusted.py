# -*- coding: utf-8 -*-
"""
06_M35_universality_publication_adjusted.py

Purpose:
Assess whether the observed enrichment of universal genes among hubs
remains after adjusting for gene publication frequency.

Input:
- study_level_unique_gene_cancer_pmid.xlsx

Outputs:
- m35_universality_gene_summary.xlsx
- m35_universality_stats.txt
- m35_pubcount_by_universality.png
- m35_pubcount_by_universality.pdf
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings("ignore")

# =========================
# 0. Paths
# =========================
input_file = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\0_data\study_level_unique_gene_cancer_pmid.xlsx"
output_dir = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M35_universality_bias"
os.makedirs(output_dir, exist_ok=True)

# =========================
# 1. Load and restrict to primary analysis
# =========================
df = pd.read_excel(input_file)
df.columns = [str(c).strip() for c in df.columns]

required_cols = {"Source", "harmonized_cancer_name", "PMID", "include_in_primary_analysis"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df = df[df["include_in_primary_analysis"].astype(str).str.strip().str.lower() == "yes"].copy()
df = df[["Source", "harmonized_cancer_name", "PMID"]].copy()
df.columns = ["Gene", "Cancer", "PMID"]

df["Gene"] = df["Gene"].astype(str).str.strip()
df["Cancer"] = df["Cancer"].astype(str).str.strip()
df["PMID"] = df["PMID"].astype(str).str.strip()

# unique study-level records
df = df.drop_duplicates(subset=["Gene", "Cancer", "PMID"]).reset_index(drop=True)

# =========================
# 2. Gene-level summary
# =========================
edge_df = df[["Gene", "Cancer"]].drop_duplicates().copy()

degree_df = (
    edge_df.groupby("Gene")["Cancer"]
    .nunique()
    .reset_index(name="Simple_degree")
)

pub_df = (
    df.groupby("Gene")["PMID"]
    .nunique()
    .reset_index(name="Unique_PMID_count")
)

gene_df = degree_df.merge(pub_df, on="Gene", how="left")
gene_df["Unique_PMID_count"] = gene_df["Unique_PMID_count"].fillna(0).astype(int)

gene_df["Universal"] = (gene_df["Simple_degree"] >= 2).astype(int)
gene_df["Hub_degree_ge_3"] = (gene_df["Simple_degree"] >= 3).astype(int)
gene_df["log2_pubcount_plus1"] = np.log2(gene_df["Unique_PMID_count"] + 1)

gene_df = gene_df.sort_values(
    by=["Simple_degree", "Unique_PMID_count", "Gene"],
    ascending=[False, False, True]
).reset_index(drop=True)

# =========================
# 3. Fisher exact test
# =========================
a = int(((gene_df["Hub_degree_ge_3"] == 1) & (gene_df["Universal"] == 1)).sum())
b = int(((gene_df["Hub_degree_ge_3"] == 1) & (gene_df["Universal"] == 0)).sum())
c = int(((gene_df["Hub_degree_ge_3"] == 0) & (gene_df["Universal"] == 1)).sum())
d = int(((gene_df["Hub_degree_ge_3"] == 0) & (gene_df["Universal"] == 0)).sum())

contingency = [[a, b], [c, d]]
fisher_or, fisher_p = fisher_exact(contingency)

# =========================
# 4. Publication count comparison
# =========================
pub_univ = gene_df.loc[gene_df["Universal"] == 1, "Unique_PMID_count"]
pub_nonuniv = gene_df.loc[gene_df["Universal"] == 0, "Unique_PMID_count"]

mw_stat, mw_p = mannwhitneyu(pub_univ, pub_nonuniv, alternative="two-sided")

# =========================
# 5. Logistic regression adjusted for publication count
#    Outcome: Universal
#    Predictors: Hub status + log2(pubcount+1)
# =========================
X = gene_df[["Hub_degree_ge_3", "log2_pubcount_plus1"]].copy()
X = sm.add_constant(X)
y = gene_df["Universal"].copy()

logit_method = "statsmodels_logit"
logit_summary_text = ""
logit_result_rows = []

try:
    model = sm.Logit(y, X)
    result = model.fit(disp=False)

    params = result.params
    conf = result.conf_int()
    pvals = result.pvalues

    for var in params.index:
        logit_result_rows.append({
            "Variable": var,
            "Coefficient": params[var],
            "OR": np.exp(params[var]),
            "CI_lower_OR": np.exp(conf.loc[var, 0]),
            "CI_upper_OR": np.exp(conf.loc[var, 1]),
            "P_value": pvals[var]
        })

    logit_summary_text = result.summary2().as_text()

except Exception as e:
    # fallback to penalized logistic regression
    logit_method = f"penalized_logistic_regression_sklearn_fallback ({str(e)})"

    X2 = gene_df[["Hub_degree_ge_3", "log2_pubcount_plus1"]].copy()
    y2 = gene_df["Universal"].copy()

    clf = LogisticRegression(
        penalty="l2",
        C=1.0,
        solver="liblinear",
        max_iter=1000
    )
    clf.fit(X2, y2)

    coef = clf.coef_[0]
    intercept = clf.intercept_[0]

    result_map = {
        "const": intercept,
        "Hub_degree_ge_3": coef[0],
        "log2_pubcount_plus1": coef[1]
    }

    for var, val in result_map.items():
        logit_result_rows.append({
            "Variable": var,
            "Coefficient": val,
            "OR": np.exp(val),
            "CI_lower_OR": np.nan,
            "CI_upper_OR": np.nan,
            "P_value": np.nan
        })

    logit_summary_text = (
        "Standard logistic regression encountered separation or convergence issues.\n"
        "A penalized logistic regression (L2 regularization, sklearn) was therefore used as a sensitivity model.\n"
        f"Intercept: {intercept:.4f}\n"
        f"Hub_degree_ge_3 coef: {coef[0]:.4f}, OR={np.exp(coef[0]):.4f}\n"
        f"log2_pubcount_plus1 coef: {coef[1]:.4f}, OR={np.exp(coef[1]):.4f}\n"
    )

logit_df = pd.DataFrame(logit_result_rows)

# =========================
# 6. Save tables
# =========================
excel_file = os.path.join(output_dir, "m35_universality_gene_summary.xlsx")

with pd.ExcelWriter(excel_file, engine="openpyxl") as writer:
    gene_df.to_excel(writer, sheet_name="gene_level_data", index=False)
    pd.DataFrame(contingency, index=["Hub", "Non_hub"], columns=["Universal", "Non_universal"]).to_excel(
        writer, sheet_name="fisher_contingency"
    )
    logit_df.to_excel(writer, sheet_name="logistic_results", index=False)

# =========================
# 7. Plot publication count by universality
# =========================
plot_png = os.path.join(output_dir, "m35_pubcount_by_universality.png")
plot_pdf = os.path.join(output_dir, "m35_pubcount_by_universality.pdf")

fig, ax = plt.subplots(figsize=(6.8, 5.8), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

data_to_plot = [
    gene_df.loc[gene_df["Universal"] == 0, "Unique_PMID_count"].values,
    gene_df.loc[gene_df["Universal"] == 1, "Unique_PMID_count"].values
]

box = ax.boxplot(
    data_to_plot,
    patch_artist=True,
    widths=0.55,
    labels=["Non-universal\n(degree = 1)", "Universal\n(degree ≥ 2)"]
)

colors = ["#9ecae1", "#fb6a4a"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.85)
    patch.set_edgecolor("black")

for item in ["whiskers", "caps", "medians"]:
    for artist in box[item]:
        artist.set_color("black")
        artist.set_linewidth(1.0)

ax.set_ylabel("Unique supporting PMIDs per gene", fontsize=12)
ax.set_title("Publication frequency by universality status", fontsize=13, pad=12)
ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.3)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

annot = f"Mann–Whitney U p = {mw_p:.2e}"
ax.text(
    0.98, 0.97, annot,
    transform=ax.transAxes,
    ha="right", va="top",
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="gray", alpha=0.9)
)

plt.tight_layout()
plt.savefig(plot_png, bbox_inches="tight", facecolor="white")
plt.savefig(plot_pdf, bbox_inches="tight", facecolor="white")
plt.close()

# =========================
# 8. Response-ready summary
# =========================
txt_file = os.path.join(output_dir, "m35_universality_stats.txt")

hub_total = int((gene_df["Hub_degree_ge_3"] == 1).sum())
hub_univ = int(((gene_df["Hub_degree_ge_3"] == 1) & (gene_df["Universal"] == 1)).sum())
nonhub_total = int((gene_df["Hub_degree_ge_3"] == 0).sum())
nonhub_univ = int(((gene_df["Hub_degree_ge_3"] == 0) & (gene_df["Universal"] == 1)).sum())

hub_pct = 100 * hub_univ / hub_total if hub_total > 0 else np.nan
nonhub_pct = 100 * nonhub_univ / nonhub_total if nonhub_total > 0 else np.nan

text = f"""
Universality and publication-frequency analysis summary

Definitions:
- Universal gene: associated with >= 2 cancer types in the revised primary simple network
- Hub gene: degree >= 3

Counts:
- Hub genes universal: {hub_univ}/{hub_total} ({hub_pct:.1f}%)
- Non-hub genes universal: {nonhub_univ}/{nonhub_total} ({nonhub_pct:.1f}%)

Fisher exact test:
- Contingency table = [[{a}, {b}], [{c}, {d}]]
- Odds ratio = {fisher_or:.4g}
- P value = {fisher_p:.3g}

Publication-frequency comparison:
- Mann–Whitney U p value for unique PMID count between universal and non-universal genes = {mw_p:.3g}

Adjusted model:
- Outcome: Universal gene status
- Predictors: hub status (degree >= 3) + log2(unique PMID count + 1)
- Method: {logit_method}

Model coefficients / odds ratios:
{logit_df.to_string(index=False)}

Interpretation:
- This model evaluates whether the observed enrichment of universal genes among hubs is partly explained by publication frequency.
- Because universal status and hub status are both derived from degree, the regression should be interpreted as a publication-adjusted sensitivity analysis rather than an independent causal model.
""".strip()

with open(txt_file, "w", encoding="utf-8") as f:
    f.write(text)

print("Done.")
print(f"Files saved to: {output_dir}")
print(f"- {excel_file}")
print(f"- {txt_file}")
print(f"- {plot_png}")
print(f"- {plot_pdf}")
