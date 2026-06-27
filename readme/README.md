# PCD Differential Expression Analysis — GSE25186

## Overview

Differential expression analysis of Primary Ciliary Dyskinesia (PCD) transcriptome data from GEO dataset GSE25186 (Illumina HumanHT-12 v4 microarray). The script performs two complementary analyses to test robustness of findings.

## Analyses

### Analysis 1 — Primary (6 PCD vs 9 Control)

All 15 samples included. Sex-corrected limma empirical Bayes. This is the primary analysis used for DEG discovery and all headline statistics.

### Analysis 2 — Sensitivity (5 PCD vs 8 Control)

Two samples excluded based on ciliary marker expression quality control. Sex-corrected limma empirical Bayes. Tests whether key gene effect sizes are robust to removal of biologically ambiguous samples.

**Excluded samples:**

**GSM1289627 (PCD)** — Expression profile indistinguishable from healthy controls across key PCD-associated genes. Six of seven disease-relevant genes show control-range expression:

| Gene   | GSM1289627 | Other PCD mean | Control mean | Profile |
|--------|-----------|----------------|--------------|---------|
| DNAH5  | 6.94      | 2.05           | 4.93         | Control |
| DNAH11 | 6.26      | 2.55           | 5.05         | Control |
| MED13L | 4.47      | 1.98           | 4.82         | Control |
| DDX58  | 4.94      | 2.31           | 5.47         | Control |
| HLA-C  | 3.81      | 2.16           | 4.95         | Control |
| CXCL9  | 4.89      | 2.59           | 5.69         | Control |
| FRMPD3 | 3.57      | 4.41           | 2.84         | PCD     |

Possible explanations: very mild phenotype, sample collection error, or clinical misdiagnosis.

**GSM1289637 (Control)** — Near-absent expression of ciliary dynein genes atypical for a healthy donor:

| Gene   | GSM1289637 | Other Ctrl mean | PCD mean | Profile  |
|--------|-----------|-----------------|----------|----------|
| DNAH5  | 0.50      | 5.44            | 2.87     | PCD-like |
| DNAH11 | 1.17      | 5.65            | 3.16     | PCD-like |

Near-zero ciliary dynein in a "healthy" control suggests subclinical PCD, carrier status, or technical artifact.

With a total cohort of 15 samples, exclusions were restricted to the two most biologically unambiguous outliers to preserve statistical power.

## Files

```
data/
  GSE25186_gene_expression.csv       Gene-level expression (25,159 genes x 15 samples)
  sample_metadata_primary.csv        All 15 samples — Group + Sex
  sample_metadata_sensitivity.csv    13 QC-filtered samples — Group + Sex
run_analysis.R                       Runs both analyses + comparison
README.md                           This file
```

## How to Run

```bash
Rscript run_analysis.R
```

Requires R with Bioconductor `limma`. The script installs limma automatically if not present.

## Output

The script prints three sections:
1. **Primary analysis** — key gene results, DEG counts, top 20 DEGs
2. **Sensitivity analysis** — same output on QC-filtered samples
3. **Comparison** — side-by-side logFC with agreement percentage and direction consistency

Two CSV files are saved:
- `results_primary_6v9.csv` — full DE results from all 15 samples
- `results_sensitivity_5v8.csv` — full DE results from 13 QC-filtered samples

## Method

- Expression: log2-transformed (floor at 1)
- DE engine: limma empirical Bayes (Smyth 2004)
- Design: ~ Sex + Group (Control as reference)
- Multiple testing: Benjamini-Hochberg FDR
- Sex included as covariate to reduce unexplained variance
- Sex determined from expression of XIST (>3 = Female) and RPS4Y1 (>5 = Male)

## Gene Name Notes

- **GSTT2B** is annotated as **GSTT2** on the Illumina HumanHT-12 v4 platform
- **METTL26** is annotated as **C1orf194** on the Illumina HumanHT-12 v4 platform