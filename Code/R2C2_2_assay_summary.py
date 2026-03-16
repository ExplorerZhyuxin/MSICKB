import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 0) Path config (EDIT IF NEEDED)
# ============================================================
BASE_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\22"

# Input files
SAMPLE_XLSX = os.path.join(BASE_DIR, "sample.xlsx")
STUDY_LEVEL_XLSX = os.path.join(BASE_DIR, "study_level_unique_gene_cancer_pmid.xlsx")

# Output files
OUT_SUMMARY_XLSX = os.path.join(BASE_DIR, "pmid_assay_summary.xlsx")
OUT_COUNTS_XLSX  = os.path.join(BASE_DIR, "assay_group_counts.xlsx")
OUT_BAR_PNG      = os.path.join(BASE_DIR, "assay_group_barplot.png")
OUT_STUDY_MERGED = os.path.join(BASE_DIR, "study_level_with_assay_group.xlsx")


# ============================================================
# 1) Helper functions
# ============================================================
def normalize_text(x) -> str:
    """Lowercase + normalize separators; safe for NaN."""
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    # unify separators
    s = s.replace("&", " and ")
    s = s.replace("/", " or ")
    s = re.sub(r"[;,|]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s

def detect_modalities(text: str) -> dict:
    t = normalize_text(text)

    # IHC: require explicit IHC/immunohistochem/staining/protein loss
    has_ihc = bool(re.search(r"\bihc\b|immunohistochem|stain|staining|protein loss|loss of .*mmr", t))

    # PCR: require explicit PCR/panel/marker keywords
    has_pcr = bool(re.search(
        r"\bpcr\b|bethesda|nci panel|promega|multiplex pcr|fragment analysis|"
        r"bat-?25|bat-?26|nr-?21|nr-?24|mono-?27|d2s123|d5s346|d17s250",
        t
    ))

    # NGS: require explicit NGS/WES/msisensor etc.
    has_ngs = bool(re.search(
        r"\bngs\b|next generation sequencing|\bwes\b|whole exome|whole-exome|"
        r"msisensor|mantis|msi[_\s-]?sensor|msi score",
        t
    ))

    return {"ihc": has_ihc, "pcr": has_pcr, "ngs": has_ngs}


def assign_group(mods: dict) -> str:
    """Assign assay_group given modalities flags."""
    n_true = sum(bool(v) for v in mods.values())
    if n_true == 0:
        return "Unclear/Not reported"
    if n_true >= 2:
        return "Mixed/Multiple"
    if mods["pcr"]:
        return "PCR-only"
    if mods["ihc"]:
        return "IHC-only"
    if mods["ngs"]:
        return "NGS-only"
    return "Unclear/Not reported"


# ============================================================
# 2) Load data
# ============================================================
if not os.path.exists(SAMPLE_XLSX):
    raise FileNotFoundError(f"Cannot find sample file: {SAMPLE_XLSX}")
if not os.path.exists(STUDY_LEVEL_XLSX):
    raise FileNotFoundError(f"Cannot find study-level file: {STUDY_LEVEL_XLSX}")

sample_df = pd.read_excel(SAMPLE_XLSX)
study_df  = pd.read_excel(STUDY_LEVEL_XLSX)

# Basic checks
if "PMID" not in sample_df.columns:
    raise ValueError("sample.xlsx must contain a 'PMID' column.")
if "test_method" not in sample_df.columns:
    raise ValueError("sample.xlsx must contain a 'test_method' column.")
if "PMID" not in study_df.columns:
    raise ValueError("study_level_unique_gene_cancer_pmid.xlsx must contain a 'PMID' column.")

# Standardize PMID as string (strip .0 etc.)
def pmid_to_str(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    # if Excel reads as float like 38691939.0
    s = re.sub(r"\.0$", "", s)
    return s

sample_df["PMID"] = sample_df["PMID"].apply(pmid_to_str)
study_df["PMID"]  = study_df["PMID"].apply(pmid_to_str)

# Only keep PMIDs that appear in the gene-network study-level file
pmids_in_network = set(study_df["PMID"].dropna().astype(str))
sample_sub = sample_df[sample_df["PMID"].isin(pmids_in_network)].copy()

# ============================================================
# 3) Build PMID-level assay summary (union across rows)
# ============================================================
# Combine multiple columns for detection: test_method + panel/platform + other_information (if exist)
cols_for_text = ["test_method"]  # ONLY use test_method for assay classification (robust & consistent)


def row_text_concat(row) -> str:
    parts = []
    for c in cols_for_text:
        parts.append(normalize_text(row.get(c, "")))
    return " | ".join([p for p in parts if p])

sample_sub["_detect_text"] = sample_sub.apply(row_text_concat, axis=1)

# detect modalities per row
mods_row = sample_sub["_detect_text"].apply(detect_modalities).apply(pd.Series)
sample_sub = pd.concat([sample_sub, mods_row], axis=1)

# PMID-level union
pmid_union = (
    sample_sub.groupby("PMID")[["ihc", "pcr", "ngs"]]
    .max()
    .reset_index()
)

pmid_union["assay_group"] = pmid_union.apply(lambda r: assign_group({"ihc": r["ihc"], "pcr": r["pcr"], "ngs": r["ngs"]}), axis=1)

# For transparency: keep raw test_method strings aggregated
raw_methods = (
    sample_sub.groupby("PMID")["test_method"]
    .apply(lambda x: "; ".join(sorted(set([str(v) for v in x.dropna().tolist()]))))
    .reset_index()
    .rename(columns={"test_method": "test_method_raw_values"})
)

pmid_summary = pmid_union.merge(raw_methods, on="PMID", how="left")

# Coverage check: PMIDs in network but not found in sample.xlsx subset
pmid_summary_set = set(pmid_summary["PMID"])
missing_pmids = sorted(list(pmids_in_network - pmid_summary_set))
missing_df = pd.DataFrame({"PMID": missing_pmids})
missing_df["assay_group"] = "Unclear/Not reported"
missing_df["test_method_raw_values"] = ""

# append missing as unclear
pmid_summary_full = pd.concat([pmid_summary, missing_df], ignore_index=True)

# ensure one row per PMID
pmid_summary_full = pmid_summary_full.drop_duplicates(subset=["PMID"]).reset_index(drop=True)

# ============================================================
# 4) Output summary + counts + plot
# ============================================================
pmid_summary_full = pmid_summary_full.sort_values(by=["assay_group", "PMID"]).reset_index(drop=True)
pmid_summary_full.to_excel(OUT_SUMMARY_XLSX, index=False)

counts = (
    pmid_summary_full["assay_group"]
    .value_counts()
    .rename_axis("assay_group")
    .reset_index(name="pmid_count")
)

# add percentages
total = counts["pmid_count"].sum()
counts["pmid_percent"] = (counts["pmid_count"] / total * 100).round(1)

# order
order = ["PCR-only", "IHC-only", "NGS-only", "Mixed/Multiple", "Unclear/Not reported"]
counts["assay_group"] = pd.Categorical(counts["assay_group"], categories=order, ordered=True)
counts = counts.sort_values("assay_group").reset_index(drop=True)
counts.to_excel(OUT_COUNTS_XLSX, index=False)

# Plot
plt.figure(figsize=(8, 4.5))
plt.bar(counts["assay_group"].astype(str), counts["pmid_count"], color="#4C78A8")
plt.xticks(rotation=25, ha="right")
plt.ylabel("Number of PMIDs")
plt.title("MSI ascertainment methods among gene-related studies (PMID-level)")
for i, (c, p) in enumerate(zip(counts["pmid_count"], counts["pmid_percent"])):
    plt.text(i, c + max(counts["pmid_count"])*0.01, f"{c} ({p}%)", ha="center", va="bottom", fontsize=9)
plt.tight_layout()
plt.savefig(OUT_BAR_PNG, dpi=300)
plt.close()

# ============================================================
# 5) Merge assay_group back to study-level file
# ============================================================
study_merged = study_df.merge(
    pmid_summary_full[["PMID", "assay_group"]],
    on="PMID",
    how="left"
)
study_merged["assay_group"] = study_merged["assay_group"].fillna("Unclear/Not reported")
study_merged.to_excel(OUT_STUDY_MERGED, index=False)

print("Done.")
print("Outputs:")
print(" -", OUT_SUMMARY_XLSX)
print(" -", OUT_COUNTS_XLSX)
print(" -", OUT_BAR_PNG)
print(" -", OUT_STUDY_MERGED)
