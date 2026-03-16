import os
import re
import pandas as pd

# ============================================================
# Config
# ============================================================
BASE_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\2_R2\22"

STUDY_WITH_ASSAY = os.path.join(BASE_DIR, "study_level_with_assay_group.xlsx")
PRIMARY_SIMPLE   = os.path.join(BASE_DIR, "primary_simple_gene_cancer_edges.xlsx")
SAMPLE_CSV       = os.path.join(BASE_DIR, "sample.csv")

OUT_XLSX         = os.path.join(BASE_DIR, "threshold_classic_PCR_subset_results.xlsx")

TOP_K = 9
MAIN_HUBS = ["BRAF","CD274","KRAS","MLH1","MSH2","PTEN","RNF43","TGFBR2","TP53"]

# ============================================================
# Helpers
# ============================================================
def pmid_to_str(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = re.sub(r"\.0$", "", s)
    return s

def compute_gene_degree(edge_df, gene="Gene", cancer="Cancer"):
    deg = edge_df.groupby(gene)[cancer].nunique().reset_index(name="degree")
    return deg.sort_values("degree", ascending=False)

def overlap_and_jaccard(sub_hubs, main_hubs):
    main = set(main_hubs)
    sub = set(sub_hubs)
    overlap = len(main & sub)
    jacc = overlap / len(main | sub) if (main | sub) else 0.0
    return overlap, round(jacc, 3)

def normalize_text(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    s = s.replace("≥", "greater than or equal to")
    s = s.replace(" >=", " greater than or equal to ")
    s = re.sub(r"\s+", " ", s)
    return s

def is_classic_pcr_definition(def_text: str) -> bool:
    """
    Classic MSI-H definition for 5-marker panels: >=2 of 5 unstable markers.
    We match flexible phrasings: '>=2 of 5', '2 of 5', 'two or more of the five', etc.
    """
    t = normalize_text(def_text)
    # strong patterns
    patterns = [
        r"greater than or equal to\s*2\s*of\s*5",
        r"\b2\s*/\s*5\b",
        r"\b2\s*of\s*5\b",
        r"two\s*or\s*more\s*of\s*the\s*five",
        r"two\s*or\s*more\s*markers.*\bfive\b",
        r"shifts?\s*in\s*two\s*or\s*more\s*markers"  # often used for Promega 5 mononucleotide markers
    ]
    return any(re.search(p, t) for p in patterns)

def has_pcr_method(test_method: str) -> bool:
    t = normalize_text(test_method)
    return bool(re.search(r"\bpcr\b|bethesda|nci|promega|fragment analysis|bat-?25|bat-?26|nr-?21|nr-?24|mono-?27", t))

# ============================================================
# Load data
# ============================================================
study = pd.read_excel(STUDY_WITH_ASSAY)
primary = pd.read_excel(PRIMARY_SIMPLE)
sample = pd.read_csv(SAMPLE_CSV)

# standardize PMIDs
study["PMID"] = study["PMID"].apply(pmid_to_str)
sample["PMID"] = sample["PMID"].apply(pmid_to_str)

# primary edges
primary.columns = [c.strip() for c in primary.columns]
if not set(["Gene", "Cancer"]) <= set(primary.columns):
    raise ValueError("primary_simple_gene_cancer_edges.xlsx must have columns: Gene, Cancer")
primary_edges = primary[["Gene","Cancer"]].drop_duplicates()

# main hubs from primary network (top9 by degree) and also keep the manuscript hubs list
main_deg = compute_gene_degree(primary_edges, "Gene", "Cancer")
main_top9_by_degree = main_deg.head(TOP_K)["Gene"].tolist()

# ============================================================
# Build classic-PCR PMID subset (restricted to PMIDs contributing to network)
# ============================================================
pmids_in_network = set(study["PMID"].dropna().astype(str))

sample_net = sample[sample["PMID"].isin(pmids_in_network)].copy()

# robustly require PCR in test_method AND classic definition
sample_net["test_method_txt"] = sample_net.get("test_method", "").apply(normalize_text)
sample_net["def_txt"] = sample_net.get("Microsatellte_instability_definition", "").apply(normalize_text)

sample_net["is_pcr"] = sample_net["test_method"].apply(has_pcr_method)
sample_net["is_classic_def"] = sample_net["Microsatellte_instability_definition"].apply(is_classic_pcr_definition)

classic_rows = sample_net[(sample_net["is_pcr"]) & (sample_net["is_classic_def"])].copy()

classic_pmids = sorted(classic_rows["PMID"].dropna().unique().tolist())

# ============================================================
# Filter study-level edges by classic PMIDs -> dedup to gene-cancer -> degree/hubs
# ============================================================
gene_col = "Source"
cancer_col = "harmonized_cancer_name"

sub_study = study[study["PMID"].isin(classic_pmids)].copy()

sub_edges = sub_study[[gene_col, cancer_col]].drop_duplicates().rename(
    columns={gene_col:"Gene", cancer_col:"Cancer"}
)

sub_deg = compute_gene_degree(sub_edges, "Gene", "Cancer")
sub_top9 = sub_deg.head(TOP_K)["Gene"].tolist()

# overlaps
overlap_manuscript, jacc_manuscript = overlap_and_jaccard(sub_top9, MAIN_HUBS)
overlap_mainDegree, jacc_mainDegree = overlap_and_jaccard(sub_top9, main_top9_by_degree)

# ============================================================
# Summaries for Supplement
# ============================================================
summary = pd.DataFrame([{
    "subset": "Classic PCR definition (>=2 of 5 unstable markers)",
    "n_pmids_in_subset": len(classic_pmids),
    "n_sample_rows_matched": len(classic_rows),
    "n_edges_gene_cancer": len(sub_edges),
    "n_genes": sub_edges["Gene"].nunique(),
    "n_cancers": sub_edges["Cancer"].nunique(),
    "top9_hubs": ", ".join(sub_top9),
    "overlap_with_manuscript_main_hubs(9)": overlap_manuscript,
    "jaccard_with_manuscript_main_hubs": jacc_manuscript,
    "main_top9_by_degree": ", ".join(main_top9_by_degree),
    "overlap_with_main_top9_by_degree": overlap_mainDegree,
    "jaccard_with_main_top9_by_degree": jacc_mainDegree,
    "manuscript_main_hubs": ", ".join(MAIN_HUBS),
}])

# save
with pd.ExcelWriter(OUT_XLSX) as w:
    summary.to_excel(w, index=False, sheet_name="summary")
    pd.DataFrame({"PMID": classic_pmids}).to_excel(w, index=False, sheet_name="classic_pmids")
    classic_rows[["PMID","test_method","panel_or_platform","Microsatellite_markers","Microsatellte_instability_definition"]].drop_duplicates() \
        .to_excel(w, index=False, sheet_name="matched_sample_rows")
    sub_edges.to_excel(w, index=False, sheet_name="sub_edges_gene_cancer")
    sub_deg.to_excel(w, index=False, sheet_name="sub_gene_degree")

print("Done.")
print("Saved:", OUT_XLSX)
print("Classic PMIDs:", len(classic_pmids))
print("Subset edges (gene-cancer):", len(sub_edges))
print("Top9 hubs:", sub_top9)
print("Overlap w/ manuscript hubs:", overlap_manuscript, "Jaccard:", jacc_manuscript)
