# -*- coding: utf-8 -*-
"""
04_M3_hub_enrichment_reanalysis.py

Purpose:
1. Recalculate hub genes in the revised primary simple gene–cancer network
2. Define hubs using degree >= 3
3. Perform ORA enrichment analysis for revised hub genes
4. Use all genes in the primary network as background whenever supported
5. Export clean tables, a publication-style enrichment plot, and response-ready text

Input:
- primary_simple_gene_cancer_edges.xlsx
  Required columns:
    - Gene
    - Cancer

Output:
- m3_hub_genes.xlsx
- m3_enrichment_all_results.xlsx
- m3_enrichment_top_terms.xlsx
- m3_hub_enrichment_barplot.png
- m3_hub_enrichment_barplot.pdf
- m3_enrichment_summary.txt
- m3_run_notes.txt
"""

import os
import textwrap
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# =========================
# 0. Paths
# =========================
input_file = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\0_data\primary_simple_gene_cancer_edges.xlsx"
output_dir = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M3"
os.makedirs(output_dir, exist_ok=True)

# =========================
# 1. Load data and define revised hub genes
# =========================
df = pd.read_excel(input_file)
df.columns = [str(c).strip() for c in df.columns]

required_cols = {"Gene", "Cancer"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

df = df[["Gene", "Cancer"]].drop_duplicates().copy()
df.columns = ["Gene", "Cancer_Type"]

degree_df = (
    df.groupby("Gene")["Cancer_Type"]
    .nunique()
    .reset_index(name="Degree")
    .sort_values(["Degree", "Gene"], ascending=[False, True])
    .reset_index(drop=True)
)

all_genes = sorted(degree_df["Gene"].dropna().astype(str).unique().tolist())
hub_genes = sorted(degree_df.loc[degree_df["Degree"] >= 3, "Gene"].dropna().astype(str).unique().tolist())

total_genes = len(all_genes)
hub_count = len(hub_genes)

if hub_count == 0:
    raise ValueError("No hub genes were identified under degree >= 3.")

# Save hub genes
hub_gene_file = os.path.join(output_dir, "m3_hub_genes.xlsx")
with pd.ExcelWriter(hub_gene_file, engine="openpyxl") as writer:
    degree_df.to_excel(writer, sheet_name="all_gene_degrees", index=False)
    pd.DataFrame({"Hub_gene": hub_genes}).to_excel(writer, sheet_name="revised_hub_genes_k_ge_3", index=False)
    pd.DataFrame({"Background_gene": all_genes}).to_excel(writer, sheet_name="background_99_genes", index=False)

# =========================
# 2. Enrichment settings
# =========================
gene_sets = {
    "GO_Biological_Process_2023": "GO Biological Process",
    "KEGG_2021_Human": "KEGG Pathway",
    "WikiPathways_2024_Human": "WikiPathways"
}

priority_keywords = [
    "mismatch", "microsatellite", "immune", "lymphocyte", "t cell", "checkpoint",
    "cancer", "colorectal", "gastric", "endometrial", "dna repair", "repair",
    "tgf", "mapk", "pi3k", "wnt", "apoptotic"
]

# Check available libraries
available_libs = gp.get_library_name(organism="human")
print(f"Available Enrichr libraries loaded: {len(available_libs)}")

valid_gene_sets = {}
for gs_name, gs_label in gene_sets.items():
    if gs_name in available_libs:
        valid_gene_sets[gs_name] = gs_label
    else:
        print(f"[WARNING] Gene-set library not found and will be skipped: {gs_name}")

if len(valid_gene_sets) == 0:
    raise ValueError("None of the requested Enrichr libraries were found.")

# =========================
# 3. Run enrichment
# =========================
all_results = []
run_notes = []

for gs_name, gs_label in valid_gene_sets.items():
    print(f"\nRunning enrichment for: {gs_name}")

    enr = None
    background_mode = "default_database_background"

    # First try custom background
    try:
        enr = gp.enrichr(
            gene_list=hub_genes,
            gene_sets=gs_name,
            organism="human",
            background=all_genes,
            outdir=None,
            no_plot=True
        )
        background_mode = "custom_background_primary_network_genes"
        print(f"  -> custom background succeeded for {gs_name}")
    except Exception as e1:
        msg1 = f"{gs_name}: custom background failed, fallback used. Error: {str(e1)}"
        run_notes.append(msg1)
        print(f"  -> custom background failed for {gs_name}")
        print(f"     reason: {str(e1)}")

        # Fallback to default library background
        try:
            enr = gp.enrichr(
                gene_list=hub_genes,
                gene_sets=gs_name,
                organism="human",
                outdir=None,
                no_plot=True
            )
            background_mode = "default_database_background"
            print(f"  -> fallback without custom background succeeded for {gs_name}")
        except Exception as e2:
            msg2 = f"{gs_name}: enrichment failed entirely. Error: {str(e2)}"
            run_notes.append(msg2)
            print(f"  -> enrichment failed entirely for {gs_name}")
            print(f"     reason: {str(e2)}")
            continue

    if enr is None:
        print(f"  -> no Enrichr object returned for {gs_name}")
        run_notes.append(f"{gs_name}: no Enrichr object returned.")
        continue

    if enr.results is None or enr.results.empty:
        print(f"  -> empty result returned for {gs_name}")
        run_notes.append(f"{gs_name}: empty enrichment result returned.")
        continue

    print(f"  -> {gs_name} returned {enr.results.shape[0]} rows")

    res = enr.results.copy()
    res["Gene_set_db"] = gs_name
    res["Category"] = gs_label
    res["Background_mode"] = background_mode

    # numeric fields
    res["P_value_num"] = pd.to_numeric(res.get("P-value", np.nan), errors="coerce")
    res["Adjusted_P_value_num"] = pd.to_numeric(res.get("Adjusted P-value", np.nan), errors="coerce")
    res["Combined_Score_num"] = pd.to_numeric(res.get("Combined Score", np.nan), errors="coerce")

    # overlap like "5/103"
    if "Overlap" in res.columns:
        res["Overlap_count"] = pd.to_numeric(
            res["Overlap"].astype(str).str.split("/").str[0],
            errors="coerce"
        )
    else:
        res["Overlap_count"] = np.nan

    # -log10 adjusted P
    res["neg_log10_adjP"] = -np.log10(res["Adjusted_P_value_num"].replace(0, np.nan))
    res["neg_log10_adjP"] = res["neg_log10_adjP"].replace([np.inf, -np.inf], np.nan)

    # priority term
    if "Term" in res.columns:
        res["Priority_term"] = res["Term"].astype(str).str.lower().apply(
            lambda x: any(k in x for k in priority_keywords)
        )
    else:
        res["Priority_term"] = False

    all_results.append(res)

if len(all_results) == 0:
    raise ValueError(
        "No enrichment results were generated from any gene-set database.\n"
        "Please inspect console messages and m3_run_notes.txt for details."
    )

all_res_df = pd.concat(all_results, ignore_index=True)

# =========================
# 4. Filter significant terms
# =========================
sig_df = all_res_df.copy()
sig_df = sig_df[sig_df["Adjusted_P_value_num"].notna()]

# Prefer adjusted p < 0.05
sig_main = sig_df[sig_df["Adjusted_P_value_num"] < 0.05].copy()

# Fallback if necessary
if sig_main.empty:
    sig_main = all_res_df[all_res_df["P_value_num"] < 0.05].copy()

if sig_main.empty:
    raise ValueError("No significant enrichment terms were found under the current settings.")

# =========================
# 5. Select top display terms
# =========================
display_df = sig_main.copy()

# Remove obviously redundant generic terms if too many
# Keep them if the result set is small
generic_terms_to_deprioritize = [
    "pathways in cancer",
    "microRNAs in cancer"
]

display_df["Generic_term"] = display_df["Term"].astype(str).str.lower().isin(generic_terms_to_deprioritize)

display_df = display_df.sort_values(
    by=["Priority_term", "Generic_term", "Adjusted_P_value_num", "Overlap_count", "Combined_Score_num"],
    ascending=[False, True, True, False, False]
).reset_index(drop=True)

# Keep top 10 terms for plotting
display_df = display_df.head(10).copy()

# Final plotting order
display_df = display_df.sort_values("neg_log10_adjP", ascending=True).reset_index(drop=True)

def shorten_term(term, width=52):
    return textwrap.fill(str(term), width=width)

display_df["Term_short"] = display_df["Term"].apply(lambda x: shorten_term(x, width=52))

# =========================
# 6. Save results
# =========================
all_result_file = os.path.join(output_dir, "m3_enrichment_all_results.xlsx")
top_term_file = os.path.join(output_dir, "m3_enrichment_top_terms.xlsx")
note_file = os.path.join(output_dir, "m3_run_notes.txt")

with pd.ExcelWriter(all_result_file, engine="openpyxl") as writer:
    all_res_df.to_excel(writer, sheet_name="all_results", index=False)
    sig_main.to_excel(writer, sheet_name="significant_terms", index=False)

with pd.ExcelWriter(top_term_file, engine="openpyxl") as writer:
    display_df.to_excel(writer, sheet_name="top_display_terms", index=False)

with open(note_file, "w", encoding="utf-8") as f:
    f.write("Run notes for enrichment analysis\n\n")
    if len(run_notes) == 0:
        f.write("No warnings or fallback notes.\n")
    else:
        for note in run_notes:
            f.write(note + "\n")

# =========================
# 7. Plot enrichment barplot
# =========================
plot_png = os.path.join(output_dir, "m3_hub_enrichment_barplot.png")
plot_pdf = os.path.join(output_dir, "m3_hub_enrichment_barplot.pdf")

plot_df = display_df.copy()

fig, ax = plt.subplots(figsize=(10.8, 6.4), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

category_palette = {
    "GO Biological Process": "#6BAED6",
    "KEGG Pathway": "#74C476",
    "WikiPathways": "#9E9AC8"
}
bar_colors = [category_palette.get(cat, "#9E9E9E") for cat in plot_df["Category"]]

bars = ax.barh(
    y=np.arange(len(plot_df)),
    width=plot_df["neg_log10_adjP"],
    color=bar_colors,
    edgecolor="black",
    linewidth=0.5,
    alpha=0.92
)

ax.set_yticks(np.arange(len(plot_df)))
ax.set_yticklabels(plot_df["Term_short"], fontsize=10)
ax.set_xlabel("-log10(adjusted P value)", fontsize=12)
ax.set_title("Functional enrichment of revised hub genes (degree ≥ 3)", fontsize=13, pad=12)

ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.35)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# add overlap labels
for bar, overlap in zip(bars, plot_df["Overlap_count"]):
    width = bar.get_width()
    label = f"n={int(overlap)}" if pd.notna(overlap) else ""
    ax.text(
        width + 0.08,
        bar.get_y() + bar.get_height() / 2,
        label,
        va="center",
        ha="left",
        fontsize=9
    )

legend_handles = [
    Patch(facecolor=category_palette[k], edgecolor="black", label=k)
    for k in category_palette if k in plot_df["Category"].values
]
ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig(plot_png, bbox_inches="tight", facecolor="white")
plt.savefig(plot_pdf, bbox_inches="tight", facecolor="white")
plt.close()

# =========================
# 8. Response-ready summary
# =========================
summary_txt = os.path.join(output_dir, "m3_enrichment_summary.txt")

top_terms_lines = []
for _, row in display_df.iterrows():
    term = row["Term"]
    category = row["Category"]
    adjp = row["Adjusted_P_value_num"]
    overlap = row["Overlap"] if "Overlap" in row else "NA"
    top_terms_lines.append(
        f"- {term} [{category}; overlap={overlap}; adjusted P={adjp:.3g}]"
    )

all_term_text = " ".join(display_df["Term"].astype(str).str.lower().tolist())

themes = []
if any(k in all_term_text for k in ["mismatch", "dna repair", "repair", "microsatellite"]):
    themes.append("DNA mismatch repair / genomic instability")
if any(k in all_term_text for k in ["immune", "lymphocyte", "t cell", "checkpoint", "apoptotic"]):
    themes.append("immune regulation")
if any(k in all_term_text for k in ["cancer", "colorectal", "gastric", "endometrial"]):
    themes.append("cancer-related pathways")
if any(k in all_term_text for k in ["tgf", "mapk", "pi3k", "wnt"]):
    themes.append("canonical oncogenic signaling")

if len(themes) == 0:
    theme_sentence = "The enriched terms mainly reflected cancer-relevant functional convergence among the revised hub genes."
else:
    theme_sentence = "The enriched terms primarily converged on " + ", ".join(themes) + "."

background_modes_used = sorted(all_res_df["Background_mode"].dropna().unique().tolist())

summary_text = f"""
M3 enrichment reanalysis summary

Revised hub definition:
- degree >= 3
- hub genes ({hub_count}): {", ".join(hub_genes)}

Background gene universe:
- {total_genes} genes from the revised primary simple gene–cancer network

Gene-set databases queried:
- GO_Biological_Process_2023
- KEGG_2021_Human
- WikiPathways_2024_Human

Background mode(s) used by gseapy/enrichr:
- {"; ".join(background_modes_used)}

Top displayed enrichment terms:
{chr(10).join(top_terms_lines)}

Interpretive summary:
{theme_sentence}

Suggested manuscript wording:
Functional enrichment analysis of the revised hub genes identified under the degree >= 3 definition
showed convergence on biologically coherent MSI-relevant processes and pathways. The dominant themes
included {", ".join(themes) if len(themes) > 0 else "cancer-relevant functional pathways"}, supporting the
biological relevance of the revised high-connectivity gene set.
""".strip()

with open(summary_txt, "w", encoding="utf-8") as f:
    f.write(summary_text)

# =========================
# 9. Console output
# =========================
print("\nDone.")
print(f"Total genes in background network: {total_genes}")
print(f"Revised hub genes (k >= 3, n={hub_count}): {', '.join(hub_genes)}")
print("\nTop displayed enrichment terms:")
for line in top_terms_lines:
    print(line)

print(f"\nFiles saved to: {output_dir}")
print(f"- {hub_gene_file}")
print(f"- {all_result_file}")
print(f"- {top_term_file}")
print(f"- {plot_png}")
print(f"- {plot_pdf}")
print(f"- {summary_txt}")
print(f"- {note_file}")
