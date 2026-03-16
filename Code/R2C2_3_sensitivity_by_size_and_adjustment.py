import os
import numpy as np
import pandas as pd

# =========================
# 0) Paths
# =========================
BASE_DATA_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\23"
OUT_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\23"

STUDY_LEVEL_FILE = os.path.join(BASE_DATA_DIR, "study_level_unique_gene_cancer_pmid.xlsx")
PRIMARY_SIMPLE_FILE = os.path.join(BASE_DATA_DIR, "primary_simple_gene_cancer_edges.xlsx")
QUALITY_FILE = os.path.join(BASE_DATA_DIR, "2026.03.11quality_assessment_template_filled.xlsx")  # 如果不在0_data，改这里

os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# 1) Load data
# =========================
study = pd.read_excel(STUDY_LEVEL_FILE)
primary = pd.read_excel(PRIMARY_SIMPLE_FILE)
quality = pd.read_excel(QUALITY_FILE)

# normalize columns
study.columns = study.columns.str.strip()
primary.columns = primary.columns.str.strip()
quality.columns = quality.columns.str.strip()

# Ensure types
study["PMID"] = study["PMID"].astype(str).str.strip()
quality["PMID"] = quality["PMID"].astype(str).str.strip()

# Keep primary-analysis only at study-level
study_yes = study[study["include_in_primary_analysis"].astype(str).str.lower().eq("yes")].copy()

# =========================
# 2) Build a PMID-level quality table with explicit aggregation rules
#    - sample_size: max
#    - adjustment/study_design/multicenter/validation: max
# =========================
for col in ["sample_size", "adjustment", "study_design", "multicenter", "validation"]:
    if col in quality.columns:
        quality[col] = pd.to_numeric(quality[col], errors="coerce")

agg_dict = {}
if "sample_size" in quality.columns: agg_dict["sample_size"] = "max"
for col in ["adjustment", "study_design", "multicenter", "validation"]:
    if col in quality.columns: agg_dict[col] = "max"

pmid_q = quality.groupby("PMID", as_index=False).agg(agg_dict)

# =========================
# 3) Helper functions
# =========================
def dedup_gene_cancer_edges(df, gene_col="Source", cancer_col="harmonized_cancer_name"):
    out = df[[gene_col, cancer_col]].dropna().copy()
    out[gene_col] = out[gene_col].astype(str).str.strip()
    out[cancer_col] = out[cancer_col].astype(str).str.strip()
    out = out.drop_duplicates()
    return out

def top_hubs_by_degree(edge_df, gene_col="Gene", cancer_col="Cancer", topn=9):
    # degree = number of distinct cancers connected to gene
    deg = (edge_df.groupby(gene_col)[cancer_col]
           .nunique()
           .sort_values(ascending=False))
    top = deg.head(topn)
    return top

def jaccard(a, b):
    a, b = set(a), set(b)
    if len(a | b) == 0:
        return np.nan
    return len(a & b) / len(a | b)

# Main network hubs from primary_simple (already unique gene-cancer)
primary_edges = primary.rename(columns={"Gene": "Gene", "Cancer": "Cancer"})[["Gene", "Cancer"]].drop_duplicates()
main_top = top_hubs_by_degree(primary_edges, gene_col="Gene", cancer_col="Cancer", topn=9)
main_hubs = list(main_top.index)

# =========================
# 4) Merge study-level with quality annotations
# =========================
study_merged = study_yes.merge(pmid_q, on="PMID", how="left")

# For edge building in subsets: use study-level unique (gene,cancer,PMID) then dedup to (gene,cancer)
# study-level gene/cancer cols:
GENE_COL = "Source"
CANCER_COL = "harmonized_cancer_name"

# =========================
# 5) Define thresholds (median and Q3) among studies with sample_size available
# =========================
sample_avail = study_merged["sample_size"].dropna()
median_n = float(sample_avail.median()) if len(sample_avail) else np.nan
q3_n = float(sample_avail.quantile(0.75)) if len(sample_avail) else np.nan

# =========================
# 6) Run sensitivity subsets
# =========================
results = []

def run_subset(name, subset_df):
    edges = dedup_gene_cancer_edges(subset_df, gene_col=GENE_COL, cancer_col=CANCER_COL)
    edges = edges.rename(columns={GENE_COL: "Gene", CANCER_COL: "Cancer"})
    top = top_hubs_by_degree(edges, gene_col="Gene", cancer_col="Cancer", topn=9)
    hubs = list(top.index)
    row = {
        "Subset": name,
        "n_PMIDs": subset_df["PMID"].nunique(),
        "n_study_level_rows": len(subset_df),
        "n_edges_unique_gene_cancer": len(edges),
        "n_genes": edges["Gene"].nunique(),
        "n_cancers": edges["Cancer"].nunique(),
        "Top9_hubs": "; ".join(hubs),
        "Overlap_with_main_hubs_(count/9)": f"{len(set(hubs)&set(main_hubs))}/9",
        "Jaccard_vs_main_top9": jaccard(hubs, main_hubs),
    }
    # also keep overlap list for convenience
    row["Overlap_genes"] = "; ".join(sorted(set(hubs) & set(main_hubs)))
    results.append(row)

# Main reference reconstructed from study-level (should match primary_simple if consistent)
run_subset("Reference (primary_simple)", study_yes.copy())

# sample size subsets
if np.isfinite(median_n):
    run_subset(f"sample_size >= median (>= {median_n:.0f})",
               study_merged[study_merged["sample_size"].notna() & (study_merged["sample_size"] >= median_n)].copy())
if np.isfinite(q3_n):
    run_subset(f"sample_size >= Q3 (>= {q3_n:.0f})",
               study_merged[study_merged["sample_size"].notna() & (study_merged["sample_size"] >= q3_n)].copy())

# adjustment subsets
if "adjustment" in study_merged.columns:
    run_subset("adjustment >= 1",
               study_merged[study_merged["adjustment"].notna() & (study_merged["adjustment"] >= 1)].copy())
    run_subset("adjustment = 2",
               study_merged[study_merged["adjustment"].notna() & (study_merged["adjustment"] == 2)].copy())

# =========================
# 7) Save outputs
# =========================
res_df = pd.DataFrame(results)

out_xlsx = os.path.join(OUT_DIR, "R2C2_3_size_and_adjustment_sensitivity.xlsx")
with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
    res_df.to_excel(writer, index=False, sheet_name="summary")

    # also save thresholds + aggregation rules
    meta = pd.DataFrame({
        "Item": [
            "Study-level input (primary analysis only)",
            "Primary simple network input",
            "Quality table input",
            "PMID-level aggregation rule (sample_size)",
            "PMID-level aggregation rule (adjustment/design/multicenter/validation)",
            "Median sample_size (among available)",
            "Q3 sample_size (among available)",
            "Hub definition",
        ],
        "Value": [
            STUDY_LEVEL_FILE,
            PRIMARY_SIMPLE_FILE,
            QUALITY_FILE,
            "max within PMID",
            "max within PMID",
            median_n,
            q3_n,
            "Top 9 genes by degree (degree = number of distinct cancers per gene)",
        ]
    })
    meta.to_excel(writer, index=False, sheet_name="meta")

print("Saved:", out_xlsx)
print("Main top-9 hubs (primary_simple):", "; ".join(main_hubs))
print("Median sample_size:", median_n, "Q3:", q3_n)
