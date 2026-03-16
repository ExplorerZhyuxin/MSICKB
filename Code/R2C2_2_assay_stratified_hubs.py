import os
import pandas as pd

BASE_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\22"

STUDY_WITH_ASSAY = os.path.join(BASE_DIR, "study_level_with_assay_group.xlsx")
PRIMARY_SIMPLE   = os.path.join(BASE_DIR, "primary_simple_gene_cancer_edges.xlsx")  # from your main pipeline (dedup gene-cancer)

OUT_HUB_TABLE    = os.path.join(BASE_DIR, "assay_stratified_hub_overlap.xlsx")
OUT_DEG_TABLE    = os.path.join(BASE_DIR, "assay_stratified_gene_degree.xlsx")

# ----------------------------
# Load
# ----------------------------
study = pd.read_excel(STUDY_WITH_ASSAY)
primary = pd.read_excel(PRIMARY_SIMPLE)

# Standardize col names
primary.columns = [c.strip() for c in primary.columns]
study.columns = [c.strip() for c in study.columns]

# primary simple has columns: Gene, Cancer
if set(["Gene","Cancer"]) <= set(primary.columns):
    primary_edges = primary[["Gene","Cancer"]].drop_duplicates()
else:
    raise ValueError("primary_simple_gene_cancer_edges.xlsx must have columns: Gene, Cancer")

# study-level has Source (gene) and harmonized_cancer_name (cancer)
gene_col = "Source"
cancer_col = "harmonized_cancer_name"
assay_col = "assay_group"

required = {gene_col, cancer_col, assay_col}
missing = required - set(study.columns)
if missing:
    raise ValueError(f"Missing columns in study_level_with_assay_group.xlsx: {missing}")

# ----------------------------
# Helper: compute degree (cross-cancer breadth)
# ----------------------------
def compute_gene_degree(edge_df, gene="Gene", cancer="Cancer"):
    # degree = number of unique cancers per gene
    deg = edge_df.groupby(gene)[cancer].nunique().reset_index()
    deg = deg.rename(columns={cancer: "degree"})
    return deg.sort_values("degree", ascending=False)

# ----------------------------
# Main network degree + hubs
# ----------------------------
main_deg = compute_gene_degree(primary_edges, "Gene", "Cancer")
# define hubs: match your manuscript; if not sure, use top 9 to mirror "9 hubs"
TOP_K = 9
main_hubs = set(main_deg.head(TOP_K)["Gene"].tolist())

# ----------------------------
# Assay-stratified subnetworks (PMID-level filter -> dedup to gene-cancer)
# ----------------------------
assay_groups = ["PCR-only", "IHC-only", "Mixed/Multiple", "Unclear/Not reported"]
deg_rows = []
overlap_rows = []

for g in assay_groups:
    sub = study[study[assay_col] == g].copy()

    # deduplicate to gene-cancer level (align with main network definition)
    sub_edges = sub[[gene_col, cancer_col]].drop_duplicates()
    sub_edges = sub_edges.rename(columns={gene_col: "Gene", cancer_col: "Cancer"})

    # degree
    sub_deg = compute_gene_degree(sub_edges, "Gene", "Cancer")
    sub_deg["assay_group"] = g
    deg_rows.append(sub_deg)

    # hubs in this subnetwork
    sub_hubs = set(sub_deg.head(TOP_K)["Gene"].tolist()) if len(sub_deg) > 0 else set()

    # overlap with main hubs
    overlap_n = len(main_hubs & sub_hubs)
    jaccard = overlap_n / len(main_hubs | sub_hubs) if len(main_hubs | sub_hubs) > 0 else 0.0

    overlap_rows.append({
        "assay_group": g,
        "n_edges_gene_cancer": len(sub_edges),
        "n_genes": sub_edges["Gene"].nunique(),
        "n_cancers": sub_edges["Cancer"].nunique(),
        f"top{TOP_K}_hubs_in_subnet": len(sub_hubs),
        f"overlap_with_main_top{TOP_K}": overlap_n,
        "jaccard_with_main": round(jaccard, 3),
        "main_top_hubs": ", ".join(sorted(main_hubs)),
        "sub_top_hubs": ", ".join(sub_deg.head(TOP_K)["Gene"].tolist()) if len(sub_deg) else ""
    })

deg_all = pd.concat(deg_rows, ignore_index=True)
overlap_df = pd.DataFrame(overlap_rows)

# Save
with pd.ExcelWriter(OUT_HUB_TABLE) as w:
    overlap_df.to_excel(w, index=False, sheet_name="hub_overlap_summary")

deg_all.to_excel(OUT_DEG_TABLE, index=False)

print("Done.")
print("Saved:")
print(" -", OUT_HUB_TABLE)
print(" -", OUT_DEG_TABLE)
