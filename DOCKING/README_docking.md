# Molecular Docking Protocol

## Software
- **AutoDock Vina 1.2.5** (Trott & Olson, J Comput Chem 2010;31:455-461)
- **Open Babel 3.1.1** (for ligand preparation)

## Target Protein
- **MED13L** (Mediator complex subunit 13L)
- **AlphaFold structure:** AF-Q71F56-F1, model v6
- **Download:** https://alphafold.ebi.ac.uk/entry/Q71F56
- **UniProt ID:** Q71F56
- **Docking domain:** CDK8-interaction domain, residues 100–510

## Receptor Preparation
1. Download the AlphaFold PDB file (AF-Q71F56-F1-model_v6.pdb)
2. Extract residues 100–510 using PyMOL or similar:
   ```
   select domain, resi 100-510
   save MED13L_CDK8_domain.pdb, domain
   ```
3. Convert to PDBQT format:
   ```
   prepare_receptor -r MED13L_CDK8_domain.pdb -o MED13L_CDK8_domain.pdbqt
   ```
   Or using Open Babel:
   ```
   obabel MED13L_CDK8_domain.pdb -O MED13L_CDK8_domain.pdbqt -xr
   ```

## Ligand Preparation
1. Download 3D structures from PubChem or DrugBank in SDF format
2. Add hydrogens and minimise energy
3. Convert to PDBQT:
   ```
   obabel ligand.sdf -O ligand.pdbqt --gen3d -h
   ```

## Seven Compounds Docked
| Compound        | PubChem CID | MW (Da) | Notes                    |
|-----------------|-------------|---------|--------------------------|
| Dexamethasone   | 5743        | 392.5   | Corticosteroid           |
| Calcitriol      | 5280453     | 416.6   | Active vitamin D         |
| Resveratrol     | 445154      | 228.2   | Polyphenol               |
| Curcumin        | 969516      | 368.4   | Polyphenol               |
| Metformin       | 4091        | 129.2   | Biguanide                |
| EPA             | 446284      | 302.5   | Omega-3 fatty acid       |
| NAC             | 12035       | 163.2   | Thiol antioxidant        |

**Sirolimus** (MW 914 Da) was **excluded** — exceeds reliable size range for AutoDock Vina.

## Docking Parameters (vina_config.txt)
```
center_x = [calculated from domain centroid]
center_y = [calculated from domain centroid]
center_z = [calculated from domain centroid]
size_x = 40
size_y = 40
size_z = 40
spacing = 0.375
exhaustiveness = 32
num_modes = 10
energy_range = 4
```

## Running Docking
```bash
vina --receptor MED13L_CDK8_domain.pdbqt \
     --ligand dexamethasone.pdbqt \
     --config vina_config.txt \
     --out dexamethasone_docked.pdbqt \
     --log dexamethasone_log.txt
```

## Expected Results (from manuscript Table S2)

| Compound        | Binding Energy (kcal/mol) | Key Contact Residues                          |
|-----------------|---------------------------|-----------------------------------------------|
| Dexamethasone   | −6.2                      | GLN289, TYR150, ASN290, HIS162, PRO370        |
| Calcitriol      | −6.2                      | GLU154, LEU163, PHE139, ARG148, TYR150        |
| Resveratrol     | −5.9                      | ARG148, HIS162, GLU161, CYS165                |
| Curcumin        | −5.5                      | LEU163, ARG148, TYR150, GLU161, HIS162        |
| Metformin       | −4.5                      | —                                             |
| EPA             | −4.5                      | —                                             |
| NAC             | −3.7                      | (weakest affinity)                             |

## Shared Binding Pocket
Residues engaged by multiple compounds:
**LEU163, ARG148, GLU161, HIS162, PHE139, TYR150**

## Important Notes
- The MED13L structure is an **AlphaFold prediction**, not experimental.
- All docking results are **hypothesis-generating** and require experimental
  validation (co-crystallography or cryo-EM).
- GSTT2B, previously used as docking target, has **no probe on GPL6947**
  and therefore cannot be a valid target derived from GSE25186 data.
