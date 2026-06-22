# Expected Numerical Outputs

All values below match the manuscript (INSILI-D-26-00281) exactly.

## Differential Expression (Script 02)
| Metric                   | Value                     |
|--------------------------|---------------------------|
| DEGs (p<0.01, |logFC|>1) | 176                       |
| DEGs (p<0.05, |logFC|>1) | 1,014                     |
| Genes at FDR < 0.05      | 0                         |
| LOC728452 logFC           | −3.22 (p = 1.08 × 10⁻³)  |
| MED13L logFC             | −3.19                     |
| MED13L p-value           | 1.63 × 10⁻³              |
| DDX58 logFC              | −3.15                     |
| HLA-C logFC              | −3.04                     |
| CLN8 logFC               | −3.10                     |
| TP53BP2 logFC            | +3.20                     |
| SLC18A1 logFC            | +3.18                     |
| XIST logFC (after corr.) | +0.55 (p = 0.41)          |
| RPS4Y1 logFC (after corr.)| −0.15 (p = 0.77)         |
| EIF1AY logFC (after corr.)| −0.06 (p = 0.93)         |

## PCA (Script 01)
| Metric                   | Value          |
|--------------------------|----------------|
| PC1 variance explained   | 30.6%          |
| PC2 variance explained   | 23.3%          |
| PC1 vs disease ANOVA p   | 2.1 × 10⁻⁵    |

## WGCNA (Script 03)
| Metric                    | Value                    |
|---------------------------|--------------------------|
| Soft-threshold power      | β = 6                    |
| Scale-free R²             | > 0.85                   |
| Mean connectivity         | 45.2                     |
| Top hub gene              | MED13L                   |
| Module-trait significance | None (all perm. p > 0.05)|
| GSTT2B on GPL6947         | ABSENT                   |

## Nested LOOCV (Script 07)
| Metric                   | Value                    |
|--------------------------|--------------------------|
| Leaky AUC                | 0.991 ± 0.006            |
| Nested AUC               | 0.750                    |
| Nested 95% CI            | [0.43, 1.00]             |
| Permutation p            | 0.062                    |
| Leakage gap              | 0.24                     |
| Stable signature (≥8/15) | 8 genes                  |

### 8-Gene Stable Signature
PDLIM3, LOC650406, C7orf29, DTX3L, CHRNA2, SERPINB4, CXCL9, C1orf187

## Cross-Cohort Validation (Script 08)
| Metric                   | Value                    |
|--------------------------|--------------------------|
| Transfer AUC (→GSE272189)| 0.30 (below chance)      |
| Reverse AUC (→GSE25186)  | 0.49 (at chance)         |
| GSTA1 logFC (GSE272189)  | +0.46, p = 0.14          |
| GSTA2 logFC (GSE272189)  | +0.47, p = 0.35          |

## LINCS L1000 (Python script)
| Drug               | Hits | Best p-value | Score |
|--------------------|------|--------------|-------|
| Curcumin           | 76   | 8.71×10⁻³    | 2.52  |
| Resveratrol        | 18   | 7.49×10⁻³    | 5.18  |
| Ibuprofen-piconol  | 8    | 4.34×10⁻³    | 16.44 |
| Dexamethasone      | 31   | 3.68×10⁻²    | 1.81  |
| Sirolimus          | 113  | 4.65×10⁻²    | 2.15  |
| Tretinoin          | 22   | 4.21×10⁻²    | 3.45  |
| Quercetin          | 14   | 4.71×10⁻²    | 3.18  |
| Calcitriol         | 15   | 8.64×10⁻²    | 4.94  |
| NAC                | 18   | 0.37         | 0.45  |
| Metformin          | 4    | 0.81         | 0.12  |

## Molecular Docking (docking/)
| Compound        | Binding Energy (kcal/mol) |
|-----------------|---------------------------|
| Dexamethasone   | −6.2                      |
| Calcitriol      | −6.2                      |
| Resveratrol     | −5.9                      |
| Curcumin        | −5.5                      |
| Metformin       | −4.5                      |
| EPA             | −4.5                      |
| NAC             | −3.7                      |

## Experimental Validation (Fig. 8)
| Gene    | Fold Reduction | p-value    | Reference |
|---------|---------------|------------|-----------|
| CCDC40  | 2.6           | < 0.001    | HPRT      |
| DNAI1   | 7.1           | < 0.001    | HPRT      |
