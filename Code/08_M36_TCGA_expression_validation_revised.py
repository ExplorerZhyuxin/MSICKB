# -*- coding: utf-8 -*-
"""
08_M36_TCGA_expression_validation_revised.py

External molecular support of revised hub genes in TCGA MSI-associated cohorts:
Expression difference in MSI-H vs non–MSI-H tumors

Outputs:
- per-cancer expression tables
- pooled combined p-value summary
- response-ready summary txt
"""

import os
import time
import requests
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, chi2
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

def safe_post(url, payload, params=None):
    r = requests.post(url, json=payload, params=params, headers=HEADERS_POST, timeout=120)
    r.raise_for_status()
    return r.json()

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

def choose_expression_profile(study_id):
    profiles = get_molecular_profiles(study_id)
    profile_ids = [p.get("molecularProfileId", "") for p in profiles]

    preferred = [
        f"{study_id}_rna_seq_v2_mrna",
        f"{study_id}_rna_seq_mrna",
        f"{study_id}_rna_seq_v2_mrna_median_Zscores",
        f"{study_id}_mrna"
    ]

    for pid in preferred:
        if pid in profile_ids:
            return pid

    for p in profiles:
        pid = p.get("molecularProfileId", "")
        datatype = str(p.get("molecularAlterationType", "")).upper()
        if datatype == "MRNA_EXPRESSION":
            return pid

    raise ValueError(f"No expression profile found for study {study_id}. Available profiles: {profile_ids}")

def get_expression_data(study_id, sample_ids, gene_entrez):
    mrna_profile = choose_expression_profile(study_id)
    entrez_to_gene = {v: k for k, v in gene_entrez.items()}

    payload = {
        "sampleIds": list(sample_ids),
        "entrezGeneIds": list(gene_entrez.values())
    }

    data = safe_post(
        f"{BASE_URL}/molecular-profiles/{mrna_profile}/molecular-data/fetch",
        payload,
        params={"projection": "SUMMARY"}
    )

    records = []
    for item in data:
        gene = entrez_to_gene.get(item.get("entrezGeneId"))
        sample_id = item.get("sampleId")
        value = item.get("value")
        if gene is not None and sample_id is not None:
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue

            records.append({
                "sample_id": sample_id,
                "gene": gene,
                "expression": value
            })

    return pd.DataFrame(records), mrna_profile

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
def analyze_cancer_expression(cancer_key, config, gene_entrez):
    study_id = config["study_id"]
    display_name = config["display_name"]
    threshold = config["msi_threshold"]

    print(f"\n{'='*70}")
    print(f"Expression validation: {cancer_key} ({display_name})")
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
    expr_df, expression_profile = get_expression_data(study_id, all_samples, gene_entrez)

    if expr_df is None or expr_df.empty:
        print("  Warning: no expression data returned.")
        return None

    results = []
    for gene in REVISED_HUB_GENES:
        sub = expr_df[expr_df["gene"] == gene].copy()
        msih_expr = sub[sub["sample_id"].isin(msih_samples)]["expression"].dropna()
        non_msih_expr = sub[sub["sample_id"].isin(non_msih_samples)]["expression"].dropna()

        if len(msih_expr) < 5 or len(non_msih_expr) < 5:
            continue

        _, p = mannwhitneyu(msih_expr, non_msih_expr, alternative="two-sided")
        msih_med = float(np.median(msih_expr))
        non_msih_med = float(np.median(non_msih_expr))

        if msih_med > 0 and non_msih_med > 0:
            log2fc = np.log2(msih_med / non_msih_med)
        else:
            log2fc = np.nan

        results.append({
            "Cancer_key": cancer_key,
            "Cancer": display_name,
            "Study_id": study_id,
            "Expression_profile_used": expression_profile,
            "Gene": gene,
            "MSI_attribute_used": chosen_attr,
            "MSI_candidate_attributes": "; ".join(candidate_attrs),
            "MSI_threshold": threshold,
            "MSIH_n": len(msih_expr),
            "Non_MSIH_n": len(non_msih_expr),
            "MSIH_median": round(msih_med, 4),
            "Non_MSIH_median": round(non_msih_med, 4),
            "Log2FC_median_ratio": round(log2fc, 4) if pd.notna(log2fc) else np.nan,
            "P_value": p
        })

    if not results:
        print("  Warning: no analyzable genes.")
        return None

    df = pd.DataFrame(results)
    df["FDR"] = multipletests(df["P_value"], method="fdr_bh")[1]
    df["Significance"] = df["FDR"].apply(sigmark_from_fdr)
    df["P_value"] = df["P_value"].round(6)
    df["FDR"] = df["FDR"].round(6)

    out_csv = os.path.join(OUTPUT_DIR, f"{cancer_key.lower()}_revised_hub_expression_validation.csv")
    df.to_csv(out_csv, index=False)

    print(f"  Saved: {out_csv}")
    print(df[["Gene", "MSIH_median", "Non_MSIH_median", "Log2FC_median_ratio", "FDR", "Significance"]].to_string(index=False))

    return df

# =========================
# 3. Pooled summary
# =========================
def pan_expression_summary(all_results):
    rows = []
    for gene in REVISED_HUB_GENES:
        gene_rows = []
        for cancer_key, df in all_results.items():
            sub = df[df["Gene"] == gene]
            if len(sub) == 1:
                gene_rows.append(sub.iloc[0])

        if len(gene_rows) == 0:
            continue

        pvals = [max(float(r["P_value"]), 1e-300) for r in gene_rows]
        chi2_stat = -2 * np.sum(np.log(pvals))
        combined_p = 1 - chi2.cdf(chi2_stat, df=2 * len(pvals))

        valid_effect_rows = [r for r in gene_rows if pd.notna(r["Log2FC_median_ratio"])]
        if valid_effect_rows:
            weights = [int(r["MSIH_n"]) for r in valid_effect_rows]
            effect_vals = [float(r["Log2FC_median_ratio"]) for r in valid_effect_rows]
            weighted_log2fc = np.average(effect_vals, weights=weights)
        else:
            weighted_log2fc = np.nan

        rows.append({
            "Cancer_key": "PAN_MSI_TCGA",
            "Cancer": "Pooled selected MSI-associated TCGA cohorts",
            "Study_id": "UCEC + COADREAD + STAD",
            "Gene": gene,
            "N_cancers_contributed": len(gene_rows),
            "MSIH_n_total": sum(int(r["MSIH_n"]) for r in gene_rows),
            "Non_MSIH_n_total": sum(int(r["Non_MSIH_n"]) for r in gene_rows),
            "Weighted_median_log2FC_summary": round(weighted_log2fc, 4) if pd.notna(weighted_log2fc) else np.nan,
            "Combined_P_value": combined_p
        })

    pooled_df = pd.DataFrame(rows)
    if pooled_df.empty:
        return pooled_df

    pooled_df["FDR"] = multipletests(pooled_df["Combined_P_value"], method="fdr_bh")[1]
    pooled_df["Significance"] = pooled_df["FDR"].apply(sigmark_from_fdr)
    pooled_df["Combined_P_value"] = pooled_df["Combined_P_value"].round(6)
    pooled_df["FDR"] = pooled_df["FDR"].round(6)

    out_csv = os.path.join(OUTPUT_DIR, "pan_tcga_revised_hub_expression_validation.csv")
    pooled_df.to_csv(out_csv, index=False)

    print(f"\nSaved pooled expression summary: {out_csv}")
    print(pooled_df[["Gene", "Weighted_median_log2FC_summary", "Combined_P_value", "FDR", "Significance"]].to_string(index=False))

    return pooled_df

# =========================
# 4. Response-ready text
# =========================
def write_summary_txt(all_results, pooled_df):
    txt_path = os.path.join(OUTPUT_DIR, "m36_tcga_expression_validation_summary.txt")

    lines = []
    lines.append("TCGA expression validation summary")
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
        lines.append(f"- Significant genes at FDR < 0.05: {len(sig)}/{len(df)}")
        if sig:
            lines.append("- " + ", ".join(sig))
        lines.append("")

    sig_pooled = pooled_df[pooled_df["FDR"] < 0.05]["Gene"].tolist() if not pooled_df.empty else []
    lines.append(f"Pooled selected MSI-associated TCGA cohorts: significant genes at FDR < 0.05 = {len(sig_pooled)}/{len(pooled_df) if not pooled_df.empty else 0}")
    if sig_pooled:
        lines.append("- " + ", ".join(sig_pooled))
    lines.append("")

    if not pooled_df.empty:
        lines.append("Top pooled results:")
        lines.append(
            pooled_df[["Gene", "Weighted_median_log2FC_summary", "Combined_P_value", "FDR", "Significance"]]
            .sort_values(["FDR", "Combined_P_value", "Gene"])
            .to_string(index=False)
        )

    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"\nSaved summary text: {txt_path}")

# =========================
# 5. Main
# =========================
if __name__ == "__main__":
    print("Starting TCGA expression validation for revised hub genes...")

    gene_entrez = get_gene_entrez_ids(REVISED_HUB_GENES)
    print(f"Resolved Entrez IDs for {len(gene_entrez)}/{len(REVISED_HUB_GENES)} genes.")

    all_results = {}
    for cancer_key, config in CANCER_CONFIGS.items():
        try:
            df = analyze_cancer_expression(cancer_key, config, gene_entrez)
            if df is not None:
                all_results[cancer_key] = df
        except Exception as e:
            print(f"Error in {cancer_key}: {e}")

    if all_results:
        pooled_df = pan_expression_summary(all_results)
        write_summary_txt(all_results, pooled_df)

        excel_path = os.path.join(OUTPUT_DIR, "m36_tcga_expression_validation_all.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            for cancer_key, df in all_results.items():
                df.to_excel(writer, sheet_name=cancer_key[:31], index=False)
            if pooled_df is not None and not pooled_df.empty:
                pooled_df.to_excel(writer, sheet_name="PAN_MSI_TCGA", index=False)

        print(f"Saved Excel workbook: {excel_path}")

    print("\nDone.")
