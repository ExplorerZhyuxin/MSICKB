import os
import numpy as np
import pandas as pd

BASE_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\23"
OUT_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\24"

STUDY_LEVEL_FILE = os.path.join(BASE_DIR, "study_level_unique_gene_cancer_pmid.xlsx")
PRIMARY_SIMPLE_FILE = os.path.join(BASE_DIR, "primary_simple_gene_cancer_edges.xlsx")

study = pd.read_excel(STUDY_LEVEL_FILE)
primary = pd.read_excel(PRIMARY_SIMPLE_FILE)

study.columns = study.columns.str.strip()
primary.columns = primary.columns.str.strip()

# Primary-analysis only
study["PMID"] = study["PMID"].astype(str).str.strip()
study_yes = study[study["include_in_primary_analysis"].astype(str).str.lower().eq("yes")].copy()

# Cancer publication volume = number of unique PMIDs per cancer
CANCER_COL = "harmonized_cancer_name"
GENE_COL = "Source"

cancer_pmid = (study_yes[[CANCER_COL, "PMID"]]
               .dropna()
               .assign(**{CANCER_COL: lambda d: d[CANCER_COL].astype(str).str.strip()})
               .drop_duplicates()
               .groupby(CANCER_COL)["PMID"].nunique()
               .rename("n_PMIDs_cancer")
               .reset_index())

# Primary simple edges: unique (Gene, Cancer)
primary_edges = primary[["Gene", "Cancer"]].dropna().drop_duplicates().copy()
primary_edges["Gene"] = primary_edges["Gene"].astype(str).str.strip()
primary_edges["Cancer"] = primary_edges["Cancer"].astype(str).str.strip()

# Attach cancer volume
edges_w = primary_edges.merge(cancer_pmid, left_on="Cancer", right_on=CANCER_COL, how="left")
edges_w = edges_w.drop(columns=[CANCER_COL])

# Sanity: any cancers missing PMID counts?
missing = edges_w["n_PMIDs_cancer"].isna().sum()

# Normalized degree: sum(1 / n_PMIDs_cancer) across connected cancers
edges_w["inv_cancer_pmid"] = 1.0 / edges_w["n_PMIDs_cancer"].astype(float)

norm_deg = (edges_w.groupby("Gene")["inv_cancer_pmid"]
            .sum()
            .sort_values(ascending=False)
            .rename("normalized_degree")
            .reset_index())

# Standard degree too (for reference)
deg = (primary_edges.groupby("Gene")["Cancer"].nunique()
       .sort_values(ascending=False)
       .rename("degree")
       .reset_index())

summary = deg.merge(norm_deg, on="Gene", how="left")

top9_norm = summary.sort_values(["normalized_degree", "degree"], ascending=False).head(9)
top9_deg = summary.sort_values(["degree", "normalized_degree"], ascending=False).head(9)

main_hubs = list(top9_deg["Gene"])
norm_hubs = list(top9_norm["Gene"])

def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else np.nan

overlap = sorted(set(main_hubs) & set(norm_hubs))

out_xlsx = os.path.join(OUT_DIR, "R2C2_4_cancer_volume_normalized_degree.xlsx")
with pd.ExcelWriter(out_xlsx) as writer:
    cancer_pmid.sort_values("n_PMIDs_cancer", ascending=False).to_excel(writer, index=False, sheet_name="cancer_PMID_volume")
    summary.sort_values("normalized_degree", ascending=False).to_excel(writer, index=False, sheet_name="gene_degree_summary")
    top9_deg.to_excel(writer, index=False, sheet_name="top9_degree")
    top9_norm.to_excel(writer, index=False, sheet_name="top9_normalized_degree")

    meta = pd.DataFrame({
        "Item": [
            "Definition (normalized_degree)",
            "Cancer volume definition",
            "Primary simple network input",
            "Study-level input (primary-analysis only)",
            "Missing cancer PMID counts (should be 0)",
            "Main top9 hubs (degree)",
            "Top9 hubs (normalized_degree)",
            "Overlap (count/9)",
            "Jaccard"
        ],
        "Value": [
            "sum over connected cancers of 1 / n_PMIDs(cancer)",
            "n_PMIDs(cancer) computed from study_level unique cancer–PMID (include_in_primary_analysis = yes)",
            PRIMARY_SIMPLE_FILE,
            STUDY_LEVEL_FILE,
            missing,
            "; ".join(main_hubs),
            "; ".join(norm_hubs),
            f"{len(overlap)}/9",
            jaccard(main_hubs, norm_hubs)
        ]
    })
    meta.to_excel(writer, index=False, sheet_name="meta")

print("Saved:", out_xlsx)
print("Missing cancer PMID counts:", missing)
print("Main hubs:", "; ".join(main_hubs))
print("Normalized hubs:", "; ".join(norm_hubs))
print("Overlap:", f"{len(overlap)}/9", "Jaccard:", jaccard(main_hubs, norm_hubs))
