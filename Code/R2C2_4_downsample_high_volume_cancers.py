import os
import numpy as np
import pandas as pd
from collections import Counter

# =========================
# Paths
# =========================
BASE_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\23"
OUT_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\24"

STUDY_LEVEL_FILE = os.path.join(BASE_DIR, "study_level_unique_gene_cancer_pmid.xlsx")
PRIMARY_SIMPLE_FILE = os.path.join(BASE_DIR, "primary_simple_gene_cancer_edges.xlsx")

# =========================
# Parameters
# =========================
N_ITER = 1000
RANDOM_SEED = 20260313
TOP_N = 9

# =========================
# Load data
# =========================
study = pd.read_excel(STUDY_LEVEL_FILE)
primary = pd.read_excel(PRIMARY_SIMPLE_FILE)

study.columns = study.columns.str.strip()
primary.columns = primary.columns.str.strip()

# Keep only primary-analysis studies
study["PMID"] = study["PMID"].astype(str).str.strip()
study["include_in_primary_analysis"] = study["include_in_primary_analysis"].astype(str).str.strip().str.lower()

study_yes = study[study["include_in_primary_analysis"].eq("yes")].copy()

# Standardize key columns
study_yes["Source"] = study_yes["Source"].astype(str).str.strip()
study_yes["harmonized_cancer_name"] = study_yes["harmonized_cancer_name"].astype(str).str.strip()
study_yes["PMID"] = study_yes["PMID"].astype(str).str.strip()

primary["Gene"] = primary["Gene"].astype(str).str.strip()
primary["Cancer"] = primary["Cancer"].astype(str).str.strip()

# Unique study-level evidence table
study_ucp = study_yes[["Source", "harmonized_cancer_name", "PMID"]].dropna().drop_duplicates().copy()
study_ucp.columns = ["Gene", "Cancer", "PMID"]

# Primary main hubs from degree
primary_edges = primary[["Gene", "Cancer"]].dropna().drop_duplicates().copy()

main_deg = (
    primary_edges.groupby("Gene")["Cancer"]
    .nunique()
    .sort_values(ascending=False)
    .reset_index(name="degree")
)

main_top = main_deg.head(TOP_N)["Gene"].tolist()
main_top_set = set(main_top)

# Cancer-specific PMID volume
cancer_pmid = (
    study_ucp.groupby("Cancer")["PMID"]
    .nunique()
    .sort_values(ascending=False)
    .reset_index(name="n_PMIDs")
)

# Q3 cap
cap_value = float(cancer_pmid["n_PMIDs"].quantile(0.75))
cap_int = int(np.floor(cap_value))  # conservative integer cap

# Cancer -> PMID list
cancer_to_pmids = (
    study_ucp[["Cancer", "PMID"]]
    .drop_duplicates()
    .groupby("Cancer")["PMID"]
    .apply(list)
    .to_dict()
)

# PMID -> rows in study_ucp
# (We'll filter by PMID, then collapse back to unique Gene-Cancer)
rng = np.random.default_rng(RANDOM_SEED)

def top_hubs_from_sampled_pmids(sampled_pmid_set, study_ucp_df, top_n=9):
    sub = study_ucp_df[study_ucp_df["PMID"].isin(sampled_pmid_set)].copy()
    if sub.empty:
        return [], pd.DataFrame(columns=["Gene", "degree"])
    edges = sub[["Gene", "Cancer"]].drop_duplicates()
    deg = (
        edges.groupby("Gene")["Cancer"]
        .nunique()
        .sort_values(ascending=False)
        .reset_index(name="degree")
    )
    # Stable tie-break: sort by degree desc, then gene asc
    deg = deg.sort_values(["degree", "Gene"], ascending=[False, True]).reset_index(drop=True)
    top = deg.head(top_n)["Gene"].tolist()
    return top, deg

def jaccard(a, b):
    a, b = set(a), set(b)
    union = a | b
    return len(a & b) / len(union) if union else np.nan

# Identify high-volume cancers
cancer_pmid["is_high_volume"] = cancer_pmid["n_PMIDs"] > cap_int

# =========================
# Downsampling iterations
# =========================
iter_rows = []
hub_counter = Counter()
main_hub_retention = Counter()
top9_signature_counter = Counter()

for i in range(1, N_ITER + 1):
    sampled_pmids = []

    for _, row in cancer_pmid.iterrows():
        cancer = row["Cancer"]
        n_pmids = int(row["n_PMIDs"])
        pmid_list = cancer_to_pmids[cancer]

        if n_pmids > cap_int:
            chosen = rng.choice(pmid_list, size=cap_int, replace=False).tolist()
        else:
            chosen = pmid_list

        sampled_pmids.extend(chosen)

    sampled_pmids = sorted(set(sampled_pmids))

    top_hubs, deg_tbl = top_hubs_from_sampled_pmids(sampled_pmids, study_ucp, top_n=TOP_N)
    top_set = set(top_hubs)
    overlap_n = len(main_top_set & top_set)
    jac = jaccard(main_top, top_hubs)

    for g in top_hubs:
        hub_counter[g] += 1

    for g in main_top:
        if g in top_set:
            main_hub_retention[g] += 1

    sig = "; ".join(top_hubs)
    top9_signature_counter[sig] += 1

    iter_rows.append({
        "iteration": i,
        "n_sampled_PMIDs_total": len(sampled_pmids),
        "top9_hubs": sig,
        "overlap_with_main_hubs_(count/9)": f"{overlap_n}/9",
        "overlap_n": overlap_n,
        "jaccard_vs_main_top9": jac
    })

iter_df = pd.DataFrame(iter_rows)

# =========================
# Summaries
# =========================
summary_stats = pd.DataFrame({
    "Metric": [
        "Iterations",
        "Random seed",
        "Top N hubs",
        "Cancer PMID cap rule",
        "Cancer PMID cap (Q3, raw)",
        "Cancer PMID cap (integer used)",
        "Number of high-volume cancers",
        "High-volume cancers",
        "Main top9 hubs",
        "Mean sampled PMIDs per iteration",
        "Median sampled PMIDs per iteration",
        "Mean overlap_n",
        "Median overlap_n",
        "Min overlap_n",
        "Max overlap_n",
        "Mean Jaccard",
        "Median Jaccard",
        "Min Jaccard",
        "Max Jaccard"
    ],
    "Value": [
        N_ITER,
        RANDOM_SEED,
        TOP_N,
        "Downsample cancers with n_PMIDs > Q3 to Q3 (integer floor)",
        cap_value,
        cap_int,
        int(cancer_pmid["is_high_volume"].sum()),
        "; ".join(cancer_pmid.loc[cancer_pmid["is_high_volume"], "Cancer"].tolist()),
        "; ".join(main_top),
        iter_df["n_sampled_PMIDs_total"].mean(),
        iter_df["n_sampled_PMIDs_total"].median(),
        iter_df["overlap_n"].mean(),
        iter_df["overlap_n"].median(),
        iter_df["overlap_n"].min(),
        iter_df["overlap_n"].max(),
        iter_df["jaccard_vs_main_top9"].mean(),
        iter_df["jaccard_vs_main_top9"].median(),
        iter_df["jaccard_vs_main_top9"].min(),
        iter_df["jaccard_vs_main_top9"].max(),
    ]
})

main_hub_retention_df = pd.DataFrame({
    "Gene": main_top,
    "retained_in_top9_count": [main_hub_retention[g] for g in main_top],
    "retained_in_top9_freq": [main_hub_retention[g] / N_ITER for g in main_top]
}).sort_values(["retained_in_top9_freq", "Gene"], ascending=[False, True]).reset_index(drop=True)

all_hub_freq_df = pd.DataFrame({
    "Gene": list(hub_counter.keys()),
    "top9_count_across_iterations": list(hub_counter.values()),
})
all_hub_freq_df["top9_freq_across_iterations"] = all_hub_freq_df["top9_count_across_iterations"] / N_ITER
all_hub_freq_df = all_hub_freq_df.sort_values(
    ["top9_freq_across_iterations", "Gene"], ascending=[False, True]
).reset_index(drop=True)

signature_df = pd.DataFrame({
    "top9_signature": list(top9_signature_counter.keys()),
    "count": list(top9_signature_counter.values())
})
signature_df["freq"] = signature_df["count"] / N_ITER
signature_df = signature_df.sort_values(["count", "top9_signature"], ascending=[False, True]).reset_index(drop=True)

# Per-cancer downsampling plan
cancer_plan = cancer_pmid.copy()
cancer_plan["cap_applied"] = np.where(cancer_plan["n_PMIDs"] > cap_int, cap_int, cancer_plan["n_PMIDs"])
cancer_plan["n_removed_per_iteration"] = cancer_plan["n_PMIDs"] - cancer_plan["cap_applied"]

# =========================
# Save
# =========================
out_xlsx = os.path.join(OUT_DIR, "R2C2_4_downsample_high_volume_cancers.xlsx")
with pd.ExcelWriter(out_xlsx) as writer:
    summary_stats.to_excel(writer, index=False, sheet_name="summary")
    cancer_plan.to_excel(writer, index=False, sheet_name="cancer_downsample_plan")
    main_hub_retention_df.to_excel(writer, index=False, sheet_name="main_hub_retention")
    all_hub_freq_df.to_excel(writer, index=False, sheet_name="all_hub_frequency")
    signature_df.to_excel(writer, index=False, sheet_name="top9_signatures")
    iter_df.to_excel(writer, index=False, sheet_name="iteration_results")

print("Saved:", out_xlsx)
print("Main top9 hubs:", "; ".join(main_top))
print("Cancer PMID Q3 cap (raw):", cap_value, " integer cap used:", cap_int)
print("High-volume cancers:", "; ".join(cancer_pmid.loc[cancer_pmid["is_high_volume"], "Cancer"].tolist()))
print("Mean overlap_n:", round(iter_df["overlap_n"].mean(), 4))
print("Median overlap_n:", round(iter_df["overlap_n"].median(), 4))
print("Min overlap_n:", int(iter_df["overlap_n"].min()), "Max overlap_n:", int(iter_df["overlap_n"].max()))
print("Mean Jaccard:", round(iter_df["jaccard_vs_main_top9"].mean(), 4))
print("Median Jaccard:", round(iter_df["jaccard_vs_main_top9"].median(), 4))
print("Min Jaccard:", round(iter_df["jaccard_vs_main_top9"].min(), 4), "Max Jaccard:", round(iter_df["jaccard_vs_main_top9"].max(), 4))

print("\nMain hub retention frequencies:")
for _, r in main_hub_retention_df.iterrows():
    print(f"  {r['Gene']}: {r['retained_in_top9_count']}/{N_ITER} ({r['retained_in_top9_freq']:.3f})")
