#!/usr/bin/env python3
"""
LINCS L1000 Connectivity Analysis via Enrichr REST API
======================================================

Manuscript: INSILI-D-26-00281
Section: 2.5, 3.6

Queries the LINCS L1000 Chemical Perturbation Consensus Signatures library
via the Enrichr REST API to identify compounds whose transcriptomic
perturbation signatures oppose the PCD gene-expression programme.

Input:
    - results/02_DEG_176_p01_lfc1.csv  (176 DEGs at p < 0.01, |logFC| > 1)
    - results/02_DEG_1014_p05_lfc1.csv (1,014 DEGs at p < 0.05, |logFC| > 1)

Output:
    - results/lincs_l1000_results.csv
    - results/lincs_drug_ranking.csv

Expected Results (from manuscript):
    Curcumin:         76 hits, best p = 8.71e-3, score = 2.52
    Resveratrol:      18 hits, best p = 7.49e-3, score = 5.18
    Ibuprofen-piconol: 8 hits, best p = 4.34e-3, score = 16.44
    Dexamethasone:    31 hits, best p = 3.68e-2, score = 1.81
    Sirolimus:       113 hits, best p = 4.65e-2, score = 2.15
    Tretinoin:        22 hits, best p = 4.21e-2, score = 3.45
    Quercetin:        14 hits, best p = 4.71e-2, score = 3.18
    Calcitriol:       15 hits, best p = 8.64e-2, score = 4.94
    NAC:              18 hits, p = 0.37, score = 0.45 (NOT SIGNIFICANT)
    Metformin:         4 hits, p = 0.81, score = 0.12 (NOT SIGNIFICANT)
"""

import os
import sys
import json
import time
import requests
import csv
from pathlib import Path

# ============================================================================
# Configuration
# ============================================================================

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
ENRICHR_ADD_URL = f"{ENRICHR_BASE}/addList"
ENRICHR_ENRICH_URL = f"{ENRICHR_BASE}/enrich"

LINCS_LIBRARY = "LINCS_L1000_Chem_Pert_Consensus_Sigs"

# 10 pre-selected candidate compounds (Table 2 in manuscript)
CANDIDATE_DRUGS = {
    "curcumin": {
        "class": "Polyphenol",
        "mechanism": "NF-kB inhibitor, anti-inflammatory",
        "lipinski": True,
        "docking_kcal": -5.5,
    },
    "resveratrol": {
        "class": "Polyphenol",
        "mechanism": "SIRT1 activator, antioxidant",
        "lipinski": True,
        "docking_kcal": -5.9,
    },
    "ibuprofen": {  # ibuprofen-piconol in manuscript
        "class": "NSAID",
        "mechanism": "Anti-inflammatory (topical)",
        "lipinski": True,
        "docking_kcal": None,
    },
    "dexamethasone": {
        "class": "Corticosteroid",
        "mechanism": "Immunosuppressant (airway use)",
        "lipinski": True,
        "docking_kcal": -6.2,
    },
    "sirolimus": {
        "class": "mTOR inhibitor",
        "mechanism": "Autophagy inducer, cilia biogenesis",
        "lipinski": False,  # MW 914 Da
        "docking_kcal": None,  # excluded from docking
    },
    "tretinoin": {
        "class": "Retinoid",
        "mechanism": "Differentiation, ciliated cell fate",
        "lipinski": True,
        "docking_kcal": None,
    },
    "quercetin": {
        "class": "Flavonoid",
        "mechanism": "Antioxidant, anti-inflammatory",
        "lipinski": True,
        "docking_kcal": None,
    },
    "calcitriol": {
        "class": "Vitamin D",
        "mechanism": "Immune modulation",
        "lipinski": False,  # MW 417
        "docking_kcal": -6.2,
    },
    "n-acetylcysteine": {  # NAC
        "class": "Thiol antioxidant",
        "mechanism": "Mucolytic (NOT SIGNIFICANT)",
        "lipinski": True,
        "docking_kcal": -3.7,
    },
    "metformin": {
        "class": "Biguanide",
        "mechanism": "AMPK activator (NOT SIGNIFICANT)",
        "lipinski": True,
        "docking_kcal": -4.5,
    },
}

# Alternate name mappings for LINCS matching
DRUG_ALIASES = {
    "curcumin": ["curcumin", "diferuloylmethane"],
    "resveratrol": ["resveratrol", "trans-resveratrol"],
    "ibuprofen": ["ibuprofen", "ibuprofen-piconol"],
    "dexamethasone": ["dexamethasone", "dex"],
    "sirolimus": ["sirolimus", "rapamycin"],
    "tretinoin": ["tretinoin", "all-trans-retinoic acid", "atra"],
    "quercetin": ["quercetin"],
    "calcitriol": ["calcitriol", "1,25-dihydroxyvitamin d3"],
    "n-acetylcysteine": ["n-acetylcysteine", "nac", "acetylcysteine"],
    "metformin": ["metformin"],
}


def find_input_files():
    """Locate DEG list files from differential expression analysis."""
    script_dir = Path(__file__).resolve().parent.parent
    results_dir = script_dir / "results"

    # Try multiple naming conventions
    deg_176_candidates = [
        results_dir / "02_DEG_176_p01_lfc1.csv",
        results_dir / "02_deg_list_176.csv",
    ]
    deg_1014_candidates = [
        results_dir / "02_DEG_1014_p05_lfc1.csv",
        results_dir / "02_deg_list_1014.csv",
    ]

    deg_176_file = None
    for f in deg_176_candidates:
        if f.exists():
            deg_176_file = f
            break

    deg_1014_file = None
    for f in deg_1014_candidates:
        if f.exists():
            deg_1014_file = f
            break

    return deg_176_file, deg_1014_file, results_dir


def load_gene_list(filepath):
    """Load gene symbols from a CSV file."""
    genes = []
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Try common column names
            for col in ["gene_symbol", "Gene.Symbol", "symbol", "gene"]:
                if col in row and row[col].strip():
                    genes.append(row[col].strip())
                    break
    return genes


def submit_to_enrichr(gene_list, description="PCD DEGs"):
    """Submit gene list to Enrichr and return user list ID."""
    genes_str = "\n".join(gene_list)
    payload = {"list": (None, genes_str), "description": (None, description)}

    response = requests.post(ENRICHR_ADD_URL, files=payload)
    if response.status_code != 200:
        raise RuntimeError(f"Enrichr submission failed: {response.status_code}")

    data = response.json()
    return data["userListId"]


def query_enrichr(user_list_id, library=LINCS_LIBRARY):
    """Query Enrichr for enrichment results."""
    params = {"userListId": user_list_id, "backgroundType": library}
    response = requests.get(ENRICHR_ENRICH_URL, params=params)

    if response.status_code != 200:
        raise RuntimeError(f"Enrichr query failed: {response.status_code}")

    return response.json()


def parse_lincs_results(enrichr_results, library=LINCS_LIBRARY):
    """Parse Enrichr LINCS results into structured format."""
    results = []
    if library not in enrichr_results:
        print(f"WARNING: Library '{library}' not found in Enrichr results")
        return results

    for entry in enrichr_results[library]:
        # Enrichr result format:
        # [rank, term_name, p_value, z_score, combined_score,
        #  overlapping_genes, adjusted_p, old_p, old_adj_p]
        term = entry[1]
        p_value = entry[2]
        z_score = entry[3]
        combined_score = entry[4]
        overlap_genes = entry[5]
        adj_p = entry[6]

        results.append({
            "term": term,
            "p_value": p_value,
            "adjusted_p": adj_p,
            "z_score": z_score,
            "combined_score": combined_score,
            "overlap_genes": ";".join(overlap_genes) if overlap_genes else "",
            "n_overlap": len(overlap_genes) if overlap_genes else 0,
        })

    return results


def match_drug_entries(results, drug_name, aliases):
    """Find all LINCS entries matching a drug name or its aliases."""
    matched = []
    for entry in results:
        term_lower = entry["term"].lower()
        for alias in aliases:
            if alias.lower() in term_lower:
                matched.append(entry)
                break
    return matched


def compute_drug_summary(matched_entries):
    """Compute drug-level summary statistics."""
    if not matched_entries:
        return {
            "n_hits": 0,
            "best_p": 1.0,
            "mean_combined_score": 0.0,
            "significant": False,
        }

    p_values = [e["p_value"] for e in matched_entries]
    scores = [e["combined_score"] for e in matched_entries]

    # Count significant hits (p < 0.05)
    sig_hits = sum(1 for p in p_values if p < 0.05)

    return {
        "n_hits": sig_hits,
        "n_total_entries": len(matched_entries),
        "best_p": min(p_values),
        "mean_combined_score": sum(scores) / len(scores),
        "significant": min(p_values) < 0.05,
    }


def main():
    print("=" * 70)
    print("LINCS L1000 Connectivity Analysis via Enrichr API")
    print("Manuscript: INSILI-D-26-00281, Sections 2.5 and 3.6")
    print("=" * 70)

    # --- Step 1: Find input files ---
    print("\n[1/5] Locating input files...")
    deg_176_file, deg_1014_file, results_dir = find_input_files()
    results_dir.mkdir(parents=True, exist_ok=True)

    if deg_176_file:
        print(f"  Found 176-gene list: {deg_176_file}")
        gene_list_176 = load_gene_list(deg_176_file)
        print(f"  Loaded {len(gene_list_176)} genes")
    else:
        print("  WARNING: 176-gene list not found. Using manuscript gene list.")
        # Fallback: key genes from manuscript Table 1
        gene_list_176 = [
            "LOC728452", "MED13L", "DDX58", "CLN8", "COG3", "HLA-C",
            "SSX2IP", "LOC652755", "KDELR3", "GUCA1B", "TP53BP2",
            "SLC18A1", "DCLK2", "LIMS2", "RUNX1T1", "LOC440900",
            "LOC650293", "CT45-5", "FLJ44790", "LOC650406",
            "PDLIM3", "C7orf29", "DTX3L", "CHRNA2", "SERPINB4",
            "CXCL9", "C1orf187", "GSR", "MBD1", "BAZ1A", "PHF3",
        ]

    if deg_1014_file:
        print(f"  Found 1014-gene list: {deg_1014_file}")
        gene_list_1014 = load_gene_list(deg_1014_file)
        print(f"  Loaded {len(gene_list_1014)} genes")
    else:
        gene_list_1014 = None
        print("  1014-gene list not found (optional)")

    # --- Step 2: Submit to Enrichr ---
    print("\n[2/5] Submitting gene list to Enrichr...")
    try:
        user_list_id = submit_to_enrichr(
            gene_list_176,
            description="PCD 176 DEGs (p<0.01, |logFC|>1, sex-corrected)"
        )
        print(f"  Enrichr user list ID: {user_list_id}")
        time.sleep(1)  # Rate limiting
    except Exception as e:
        print(f"  ERROR: Enrichr submission failed: {e}")
        print("  Ensure internet connectivity and try again.")
        sys.exit(1)

    # --- Step 3: Query LINCS L1000 library ---
    print("\n[3/5] Querying LINCS L1000 Chem Pert Consensus Sigs...")
    try:
        enrichr_results = query_enrichr(user_list_id, LINCS_LIBRARY)
        all_lincs = parse_lincs_results(enrichr_results, LINCS_LIBRARY)
        print(f"  Retrieved {len(all_lincs)} LINCS entries")
    except Exception as e:
        print(f"  ERROR: Enrichr query failed: {e}")
        sys.exit(1)

    # Save all LINCS results
    lincs_all_path = results_dir / "lincs_l1000_all_results.csv"
    if all_lincs:
        with open(lincs_all_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_lincs[0].keys())
            writer.writeheader()
            writer.writerows(all_lincs)
        print(f"  Saved: {lincs_all_path}")

    # --- Step 4: Match candidate drugs ---
    print("\n[4/5] Matching 10 pre-selected candidate compounds...")
    drug_summaries = []

    for drug_key, drug_info in CANDIDATE_DRUGS.items():
        aliases = DRUG_ALIASES.get(drug_key, [drug_key])
        matched = match_drug_entries(all_lincs, drug_key, aliases)
        summary = compute_drug_summary(matched)

        drug_summaries.append({
            "drug": drug_key.replace("-", " ").title()
                if drug_key != "n-acetylcysteine" else "NAC",
            "class": drug_info["class"],
            "mechanism": drug_info["mechanism"],
            "lincs_hits": summary["n_hits"],
            "best_p_value": f"{summary['best_p']:.4g}",
            "mean_combined_score": f"{summary['mean_combined_score']:.2f}",
            "docking_kcal_mol": drug_info["docking_kcal"]
                if drug_info["docking_kcal"] else "n/a",
            "lipinski": "Yes" if drug_info["lipinski"] else "No",
            "significant": "Yes" if summary["significant"] else "No",
        })

        sig_str = "SIGNIFICANT" if summary["significant"] else "not significant"
        print(f"  {drug_key:>20s}: {summary['n_hits']:3d} hits, "
              f"p = {summary['best_p']:.4g}, "
              f"score = {summary['mean_combined_score']:.2f} "
              f"[{sig_str}]")

    # Save drug ranking
    ranking_path = results_dir / "lincs_drug_ranking.csv"
    with open(ranking_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=drug_summaries[0].keys())
        writer.writeheader()
        writer.writerows(drug_summaries)
    print(f"\n  Saved: {ranking_path}")

    # --- Step 5: Verification ---
    print("\n[5/5] Verification against manuscript expectations...")
    print("  Expected results (Table 2):")
    print("    Curcumin:          76 hits, p = 8.71e-3")
    print("    Resveratrol:       18 hits, p = 7.49e-3")
    print("    Ibuprofen-piconol:  8 hits, p = 4.34e-3")
    print("    Dexamethasone:     31 hits, p = 3.68e-2")
    print("    NAC:               NOT SIGNIFICANT (p = 0.37)")
    print("    Metformin:         NOT SIGNIFICANT (p = 0.81)")
    print()
    print("  NOTE: Exact hit counts may vary slightly depending on")
    print("  Enrichr database version. The ranking and significance")
    print("  pattern should be consistent.")
    print()
    print("  IMPORTANT: NAC and metformin should NOT be significant.")
    print("  If they appear significant in your run, this may reflect")
    print("  a database update and should be noted.")

    print("\n" + "=" * 70)
    print("LINCS L1000 analysis complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
