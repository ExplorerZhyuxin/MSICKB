# -*- coding: utf-8 -*-
"""
07_M36_TCGA_mutation_validation_revised.py

External molecular support of revised hub genes in TCGA MSI-associated cohorts:
Mutation prevalence in MSI-H vs non–MSI-H tumors

Outputs:
- per-cancer mutation tables
- pooled mutation summary across selected TCGA cohorts
- response-ready summary txt
"""

import os
import time
import requests
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# =========================
# 0. Config
# =========================
OUTPUT_DIR = r"G:\其它\博士\分1\CSBJ\CSBJ\大修\Major revision\1_M\M36_TCGA_validation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

BASE_URL = "https://www.cbioportal.org/api"
HEADERS_JSON = {"Accept": "application/json"}
HEADERS_POST = {"Accept": "application/json", "Content-Type": "application/json"}

REVISED_HUB_GENES = ["BRAF", "CD274", "KRAS", "MLH1", "MSH2", "PTEN", "RNF43", "TGFBR2", "TP53"]

CANCER_CONFIGS = {
    "UCEC": {
        "study_id": "ucec_tcga_pan_can_atlas_2018",
        "display_name": "Endometrial Cancer",
        "msi_threshold": 3.5
    },
    "COADREAD": {
        "study_id": "coadread_tcga_pan_can_atlas_2018",
        "display_name": "Colorectal Cancer",
        "msi_threshold": 3.5
    },
    "STAD": {
        "study_id": "stad_tcga_pan_can_atlas_2018",
        "display_name": "Gastric Cancer",
        "msi_threshold": 3.5
    }
}

# =========================
# 1. Helper functions
# =========================
def safe_get(url, params=None):
    r = requests.get(url, params=params, headers=HEADERS_JSON, timeout=60)
    r.raise_for_status()
    return r.json()

def safe_post(url, payload, params=None, max_retries=5):
    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.post(url, json=payload, params=params, headers=HEADERS_POST, timeout=180)
            r.raise_for_status()
            data = r.json()
            # basic sanity: must be list for fetch endpoints
            if data is None:
                raise ValueError("Empty JSON")
            return data
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            print(f"  Warning: POST attempt {attempt}/{max_retries} failed: {e}. Retrying in {wait}s...")
            time.sleep(wait)
    raise RuntimeError(f"POST failed after {max_retries} attempts: {last_err}")


def get_gene_entrez_ids(genes):
    gene_entrez = {}
    for gene in genes:
        try:
            resp = safe_get(f"{BASE_URL}/genes/{gene}")
            gene_entrez[gene] = resp["entrezGeneId"]
            time.sleep(0.05)
        except Exception as e:
            print(f"  Warning: failed to get Entrez ID for {gene}: {e}")
    return gene_entrez

def get_sample_clinical_data(study_id):
    return safe_get(
        f"{BASE_URL}/studies/{study_id}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"}
    )

def choose_msisensor_attr_id(clinical_data):
    """
    Choose a single MSIsensor-related clinical attribute per study.
    Excludes MANTIS.
    """
    attr_ids = sorted({str(item.get("clinicalAttributeId", "")) for item in clinical_data})

    candidates = []
    for attr_id in attr_ids:
        attr_upper = attr_id.upper().replace("_", "").replace(" ", "")
        if "MSISENSOR" in attr_upper and "MANTIS" not in attr_upper:
            candidates.append(attr_id)

    # Priority rules
    preferred_exact = [
        "MSISENSOR_SCORE",
        "MSISENSOR",
        "MSI_SENSOR_SCORE",
        "MSI_SENSOR"
    ]

    for pref in preferred_exact:
        for attr_id in candidates:
            if attr_id.upper() == pref:
                return attr_id, candidates

    # Secondary preference: contains both MSI and SCORE and SENSOR
    for attr_id in candidates:
        au = attr_id.upper().replace("_", "").replace(" ", "")
        if "MSI" in au and "SENSOR" in au and "SCORE" in au:
            return attr_id, candidates

    return (candidates[0] if candidates else None), candidates

def extract_msi_groups(clinical_data, threshold=3.5):
    """
    Identify MSI-H and non–MSI-H groups from sample-level clinical data
    using one chosen MSIsensor-related attribute.
    """
    chosen_attr, candidate_attrs = choose_msisensor_attr_id(clinical_data)

    if chosen_attr is None:
        return set(), set(), candidate_attrs, None

    msi_scores = {}
    for item in clinical_data:
        attr_id = str(item.get("clinicalAttributeId", ""))
        if attr_id != chosen_attr:
            continue

        sample_id = item.get("sampleId")
        value = item.get("value")
        try:
            msi_scores[sample_id] = float(value)
        except (TypeError, ValueError):
            continue

    msih = {s for s, score in msi_scores.items() if score >= threshold}
    non_msih = {s for s, score in msi_scores.items() if score < threshold}

    return msih, non_msih, candidate_attrs, chosen_attr

def get_molecular_profiles(study_id):
    return safe_get(f"{BASE_URL}/studies/{study_id}/molecular-profiles")

def choose_mutation_profile(study_id):
    profiles = get_molecular_profiles(study_id)
    profile_ids = [p.get("molecularProfileId", "") for p in profiles]

    preferred = [
        f"{study_id}_mutations",
        f"{study_id}_mutations_unfiltered"
    ]

    for pid in preferred:
        if pid in profile_ids:
            return pid

    for p in profiles:
        pid = p.get("molecularProfileId", "")
        datatype = str(p.get("molecularAlterationType", "")).upper()
        if datatype == "MUTATION_EXTENDED":
            return pid

    raise ValueError(f"No mutation profile found for study {study_id}. Available profiles: {profile_ids}")

def get_mutations(study_id, sample_ids, gene_entrez):
    mutation_profile = choose_mutation_profile(study_id)
    entrez_to_gene = {v: k for k, v in gene_entrez.items()}
    

    payload = {
        "sampleIds": list(sample_ids),
        "entrezGeneIds": list(gene_entrez.values())
    }

    data = safe_post(
        f"{BASE_URL}/molecular-profiles/{mutation_profile}/mutations/fetch",
        payload,
        params={"projection": "SUMMARY"}
    )
    print(f"  DEBUG: mutation fetch returned {len(data)} records for {study_id} ({mutation_profile})")
    returned_entrez = {m.get("entrezGeneId") for m in data if m.get("entrezGeneId") is not None}
    missing_entrez = set(gene_entrez.values()) - returned_entrez
    if missing_entrez:
        print(f"  Warning: missing entrez IDs in mutation response: {missing_entrez}. Will retry once...")

    gene_mutations = {gene: set() for gene in gene_entrez.keys()}
    for mut in data:
        entrez_id = mut.get("entrezGeneId")
        sample_id = mut.get("sampleId")
        gene = entrez_to_gene.get(entrez_id)
        if gene and sample_id:
            gene_mutations[gene].add(sample_id)

    return gene_mutations, mutation_profile
    

def fisher_with_haldane(msih_mut, msih_total, non_msih_mut, non_msih_total):
    msih_wt = msih_total - msih_mut
    non_msih_wt = non_msih_total - non_msih_mut

    table = [[msih_mut, msih_wt], [non_msih_mut, non_msih_wt]]
    _, p = fisher_exact(table)

    a = msih_mut + 0.5
    b = msih_wt + 0.5
    c = non_msih_mut + 0.5
    d = non_msih_wt + 0.5

    or_val = (a * d) / (b * c)
    log_or = np.log(or_val)
    se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    ci_low = np.exp(log_or - 1.96 * se)
    ci_high = np.exp(log_or + 1.96 * se)

    return table, or_val, p, ci_low, ci_high

def note_from_counts(msih_mut, non_msih_mut):
    if msih_mut == 0 and non_msih_mut == 0:
        return "No mutations in either group"
    elif non_msih_mut == 0:
        return "No mutations in non-MSI-H"
    elif msih_mut == 0:
        return "No mutations in MSI-H"
    return ""

def sigmark_from_fdr(fdr):
    if fdr < 0.001:
        return "***"
    elif fdr < 0.01:
        return "**"
    elif fdr < 0.05:
        return "*"
    return ""

# =========================
# 2. Per-cancer analysis
# =========================
def analyze_cancer_mutation(cancer_key, config, gene_entrez):
    study_id = config["study_id"]
    display_name = config["display_name"]
    threshold = config["msi_threshold"]

    print(f"\n{'='*70}")
    print(f"Mutation validation: {cancer_key} ({display_name})")
    print(f"{'='*70}")

    clinical_data = get_sample_clinical_data(study_id)
    msih_samples, non_msih_samples, candidate_attrs, chosen_attr = extract_msi_groups(
        clinical_data, threshold=threshold
    )

    print(f"  MSIsensor candidate attributes: {candidate_attrs if candidate_attrs else 'None'}")
    print(f"  Chosen MSI attribute: {chosen_attr if chosen_attr else 'None'}")
    print(f"  MSI-H samples: {len(msih_samples)}")
    print(f"  non-MSI-H samples: {len(non_msih_samples)}")

    if chosen_attr is None:
        print("  Warning: no MSIsensor attribute found, skipping.")
        return None

    if len(msih_samples) < 10 or len(non_msih_samples) < 10:
        print("  Warning: insufficient sample size, skipping.")
        return None

    all_samples = msih_samples | non_msih_samples
    gene_mutations, mutation_profile = get_mutations(study_id, all_samples, gene_entrez)

    results = []
    print("DEBUG genes:", REVISED_HUB_GENES)

    for gene in REVISED_HUB_GENES:
        print("DEBUG processing gene:", gene)

        mutated_samples = gene_mutations.get(gene, set())
        msih_mut = len(mutated_samples & msih_samples)
        non_msih_mut = len(mutated_samples & non_msih_samples)

        table, or_val, p, ci_low, ci_high = fisher_with_haldane(
            msih_mut, len(msih_samples), non_msih_mut, len(non_msih_samples)
        )

        msih_rate = msih_mut / len(msih_samples) * 100
        non_msih_rate = non_msih_mut / len(non_msih_samples) * 100

        results.append({
            "Cancer_key": cancer_key,
            "Cancer": display_name,
            "Study_id": study_id,
            "Mutation_profile_used": mutation_profile,
            "Gene": gene,
            "MSI_attribute_used": chosen_attr,
            "MSI_candidate_attributes": "; ".join(candidate_attrs),
            "MSI_threshold": threshold,
            "MSIH_mut": msih_mut,
            "MSIH_total": len(msih_samples),
            "MSIH_rate": round(msih_rate, 2),
            "Non_MSIH_mut": non_msih_mut,
            "Non_MSIH_total": len(non_msih_samples),
            "Non_MSIH_rate": round(non_msih_rate, 2),
            "Difference": round(msih_rate - non_msih_rate, 2),
            "OR_haldane": round(or_val, 4),
            "CI_lower": round(ci_low, 4),
            "CI_upper": round(ci_high, 4),
            "P_value": p,
            "Note": note_from_counts(msih_mut, non_msih_mut)
        })

    df = pd.DataFrame(results)
    df["FDR"] = multipletests(df["P_value"], method="fdr_bh")[1]
    df["Significance"] = df["FDR"].apply(sigmark_from_fdr)
    df["P_value"] = df["P_value"].round(6)
    df["FDR"] = df["FDR"].round(6)

    out_csv = os.path.join(OUTPUT_DIR, f"{cancer_key.lower()}_revised_hub_mutation_validation.csv")
    df.to_csv(out_csv, index=False)

    print(f"  Saved: {out_csv}")
    print(df[["Gene", "MSIH_rate", "Non_MSIH_rate", "OR_haldane", "FDR", "Significance"]].to_string(index=False))

    return df

# =========================
# 3. Pooled summary
# =========================
def pooled_mutation_summary(all_results):
    pooled_rows = []
    for gene in REVISED_HUB_GENES:
        sub = []
        for cancer_key, df in all_results.items():
            row = df[df["Gene"] == gene].iloc[0]
            sub.append(row)

        msih_mut = sum(int(r["MSIH_mut"]) for r in sub)
        msih_total = sum(int(r["MSIH_total"]) for r in sub)
        non_msih_mut = sum(int(r["Non_MSIH_mut"]) for r in sub)
        non_msih_total = sum(int(r["Non_MSIH_total"]) for r in sub)

        table, or_val, p, ci_low, ci_high = fisher_with_haldane(
            msih_mut, msih_total, non_msih_mut, non_msih_total
        )

        pooled_rows.append({
            "Cancer_key": "PAN_MSI_TCGA",
            "Cancer": "Pooled selected MSI-associated TCGA cohorts",
            "Study_id": "UCEC + COADREAD + STAD",
            "Gene": gene,
            "MSIH_mut": msih_mut,
            "MSIH_total": msih_total,
            "MSIH_rate": round(msih_mut / msih_total * 100, 2),
            "Non_MSIH_mut": non_msih_mut,
            "Non_MSIH_total": non_msih_total,
            "Non_MSIH_rate": round(non_msih_mut / non_msih_total * 100, 2),
            "Difference": round(msih_mut / msih_total * 100 - non_msih_mut / non_msih_total * 100, 2),
            "OR_haldane": round(or_val, 4),
            "CI_lower": round(ci_low, 4),
            "CI_upper": round(ci_high, 4),
            "P_value": p,
            "Note": note_from_counts(msih_mut, non_msih_mut)
        })

    pooled_df = pd.DataFrame(pooled_rows)
    pooled_df["FDR"] = multipletests(pooled_df["P_value"], method="fdr_bh")[1]
    pooled_df["Significance"] = pooled_df["FDR"].apply(sigmark_from_fdr)
    pooled_df["P_value"] = pooled_df["P_value"].round(6)
    pooled_df["FDR"] = pooled_df["FDR"].round(6)

    out_csv = os.path.join(OUTPUT_DIR, "pan_tcga_revised_hub_mutation_validation.csv")
    pooled_df.to_csv(out_csv, index=False)

    print(f"\nSaved pooled mutation summary: {out_csv}")
    print(pooled_df[["Gene", "MSIH_rate", "Non_MSIH_rate", "OR_haldane", "FDR", "Significance"]].to_string(index=False))

    return pooled_df

# =========================
# 4. Response-ready text
# =========================
def write_summary_txt(all_results, pooled_df):
    txt_path = os.path.join(OUTPUT_DIR, "m36_tcga_mutation_validation_summary.txt")

    lines = []
    lines.append("TCGA mutation validation summary")
    lines.append("")
    lines.append("Revised hub genes:")
    lines.append("- " + ", ".join(REVISED_HUB_GENES))
    lines.append("")
    lines.append("TCGA cohorts analyzed:")
    for k, v in CANCER_CONFIGS.items():
        lines.append(f"- {k}: {v['study_id']} ({v['display_name']})")
    lines.append("")

    for cancer_key, df in all_results.items():
        meta = df.iloc[0]
        lines.append(f"{cancer_key}:")
        lines.append(f"- MSI attribute used: {meta['MSI_attribute_used']}")
        lines.append(f"- MSI candidate attributes: {meta['MSI_candidate_attributes']}")
        lines.append(f"- Threshold for MSI-H: {meta['MSI_threshold']}")
        sig = df[df["FDR"] < 0.05]["Gene"].tolist()
        lines.append(f"- Significant genes at FDR < 0.05: {len(sig)}/{len(REVISED_HUB_GENES)}")
        if sig:
            lines.append("- " + ", ".join(sig))
        lines.append("")

    sig_pooled = pooled_df[pooled_df["FDR"] < 0.05]["Gene"].tolist()
    lines.append(f"Pooled selected MSI-associated TCGA cohorts: significant genes at FDR < 0.05 = {len(sig_pooled)}/{len(REVISED_HUB_GENES)}")
    if sig_pooled:
        lines.append("- " + ", ".join(sig_pooled))
    lines.append("")
    lines.append("Top pooled results:")
    lines.append(
        pooled_df.sort_values(["FDR", "P_value", "Gene"])[
            ["Gene", "MSIH_rate", "Non_MSIH_rate", "OR_haldane", "P_value", "FDR", "Significance"]
        ].to_string(index=False)
    )

    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"\nSaved summary text: {txt_path}")

# =========================
# 5. Main
# =========================
if __name__ == "__main__":
    print("Starting TCGA mutation validation for revised hub genes...")

    gene_entrez = get_gene_entrez_ids(REVISED_HUB_GENES)
    print(f"Resolved Entrez IDs for {len(gene_entrez)}/{len(REVISED_HUB_GENES)} genes.")

    all_results = {}
    for cancer_key, config in CANCER_CONFIGS.items():
        try:
            df = analyze_cancer_mutation(cancer_key, config, gene_entrez)
            if df is not None:
                all_results[cancer_key] = df
        except Exception as e:
            print(f"Error in {cancer_key}: {e}")

    if all_results:
        pooled_df = pooled_mutation_summary(all_results)
        write_summary_txt(all_results, pooled_df)

        excel_path = os.path.join(OUTPUT_DIR, "m36_tcga_mutation_validation_all.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            for cancer_key, df in all_results.items():
                df.to_excel(writer, sheet_name=cancer_key[:31], index=False)
            pooled_df.to_excel(writer, sheet_name="PAN_MSI_TCGA", index=False)

        print(f"Saved Excel workbook: {excel_path}")

    print("\nDone.")
