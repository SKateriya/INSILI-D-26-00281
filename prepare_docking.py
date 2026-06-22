#!/usr/bin/env python3
"""
Molecular Docking Preparation Script
=====================================

Manuscript: INSILI-D-26-00281, Section 2.5 and 3.7

Downloads the MED13L AlphaFold structure and prepares the receptor
for docking with AutoDock Vina 1.2.5.

Requirements:
    - Open Babel 3.1.1 (for format conversion)
    - AutoDock Vina 1.2.5 (for docking)
    - Internet access (for AlphaFold download)

Usage:
    python prepare_docking.py
"""

import os
import sys
import subprocess

ALPHAFOLD_URL = (
    "https://alphafold.ebi.ac.uk/files/"
    "AF-Q71F56-F1-model_v6.pdb"
)

COMPOUNDS = {
    "dexamethasone": {"cid": "5743", "mw": 392.5},
    "calcitriol":    {"cid": "5280453", "mw": 416.6},
    "resveratrol":   {"cid": "445154", "mw": 228.2},
    "curcumin":      {"cid": "969516", "mw": 368.4},
    "metformin":     {"cid": "4091", "mw": 129.2},
    "epa":           {"cid": "446284", "mw": 302.5},
    "nac":           {"cid": "12035", "mw": 163.2},
}


def download_alphafold_structure():
    """Download MED13L AlphaFold structure."""
    output = "AF-Q71F56-F1-model_v6.pdb"
    if os.path.exists(output):
        print(f"  AlphaFold structure already exists: {output}")
        return output

    print(f"  Downloading from AlphaFold: {ALPHAFOLD_URL}")
    try:
        import urllib.request
        urllib.request.urlretrieve(ALPHAFOLD_URL, output)
        print(f"  Saved: {output}")
    except Exception as e:
        print(f"  ERROR: Download failed: {e}")
        print("  Please download manually from:")
        print(f"    {ALPHAFOLD_URL}")
        sys.exit(1)
    return output


def extract_domain(pdb_file, start=100, end=510):
    """Extract CDK8-interaction domain (residues 100-510)."""
    output = "MED13L_CDK8_domain.pdb"
    print(f"  Extracting residues {start}-{end} from {pdb_file}...")

    with open(pdb_file, "r") as fin, open(output, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    resnum = int(line[22:26].strip())
                    if start <= resnum <= end:
                        fout.write(line)
                except ValueError:
                    continue
            elif line.startswith("END"):
                fout.write(line)

    print(f"  Saved domain: {output}")
    return output


def compute_centroid(pdb_file):
    """Compute centroid of the domain for grid box placement."""
    coords = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))

    if not coords:
        print("  ERROR: No ATOM records found")
        return None

    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)

    print(f"  Domain centroid: ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    print(f"  Total atoms: {len(coords)}")
    return cx, cy, cz


def update_vina_config(cx, cy, cz):
    """Update vina_config.txt with computed centroid coordinates."""
    config_lines = [
        "# AutoDock Vina 1.2.5 Configuration",
        "# Target: MED13L CDK8-interaction domain (residues 100-510)",
        "# Structure: AlphaFold AF-Q71F56-F1, model v6",
        "",
        "receptor = MED13L_CDK8_domain.pdbqt",
        "",
        f"center_x = {cx:.3f}",
        f"center_y = {cy:.3f}",
        f"center_z = {cz:.3f}",
        "",
        "size_x = 40",
        "size_y = 40",
        "size_z = 40",
        "",
        "exhaustiveness = 32",
        "num_modes = 10",
        "energy_range = 4",
        "",
        "seed = 42",
    ]

    with open("vina_config.txt", "w") as f:
        f.write("\n".join(config_lines) + "\n")
    print("  Updated vina_config.txt with centroid coordinates")


def main():
    print("=" * 60)
    print("Molecular Docking Preparation")
    print("Target: MED13L (AF-Q71F56-F1, residues 100-510)")
    print("Software: AutoDock Vina 1.2.5")
    print("=" * 60)

    # Step 1: Download AlphaFold structure
    print("\n[1/4] Downloading MED13L AlphaFold structure...")
    pdb_file = download_alphafold_structure()

    # Step 2: Extract CDK8-interaction domain
    print("\n[2/4] Extracting CDK8-interaction domain (residues 100-510)...")
    domain_pdb = extract_domain(pdb_file, start=100, end=510)

    # Step 3: Compute centroid for grid box
    print("\n[3/4] Computing domain centroid for grid placement...")
    centroid = compute_centroid(domain_pdb)
    if centroid:
        update_vina_config(*centroid)

    # Step 4: Instructions
    print("\n[4/4] Next steps...")
    print("  1. Convert receptor to PDBQT:")
    print("     obabel MED13L_CDK8_domain.pdb -O MED13L_CDK8_domain.pdbqt -xr")
    print()
    print("  2. Download ligand SDF files from PubChem:")
    for name, info in COMPOUNDS.items():
        print(f"     {name:20s} CID {info['cid']:>10s}  (MW {info['mw']:.1f} Da)")
    print()
    print("  3. Convert ligands to PDBQT:")
    print("     obabel ligand.sdf -O ligand.pdbqt --gen3d -h")
    print()
    print("  4. Run docking for each compound:")
    print("     vina --config vina_config.txt --ligand compound.pdbqt \\")
    print("          --out compound_docked.pdbqt --log compound_log.txt")
    print()
    print("  NOTE: Sirolimus (MW 914 Da) is excluded from docking.")
    print()
    print("  Expected binding energies (kcal/mol):")
    expected = [
        ("Dexamethasone", -6.2), ("Calcitriol", -6.2),
        ("Resveratrol", -5.9), ("Curcumin", -5.5),
        ("Metformin", -4.5), ("EPA", -4.5), ("NAC", -3.7),
    ]
    for name, energy in expected:
        print(f"    {name:20s} {energy:5.1f}")

    print("\n" + "=" * 60)
    print("Preparation complete. See README_docking.md for full protocol.")
    print("=" * 60)


if __name__ == "__main__":
    main()
