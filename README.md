# INSILI-D-26-00281

## Integrative Transcriptomic Analysis Reveals Molecular Signatures and Candidate Therapeutic Targets in Primary Ciliary Dyskinesia

**Manuscript:** INSILI-D-26-00281 (In Silico Research in Biomedicine)

**Authors:** Jitender, Md Wamique Hossain, Shilpa Mohanty, Suneel Kateriya

**Affiliation:** Laboratory of Optobiotechnology, School of Biotechnology, Jawaharlal Nehru University, New Delhi, India

---

## Overview

This repository contains all analysis code for the manuscript. The pipeline re-analyses the GSE25186 microarray dataset (6 PCD vs. 9 controls, Illumina HumanHT-12 V3.0, GPL6947) with three critical corrections applied to the original analysis:

1. **Sex-covariate correction** — PCD group: 4F/2M; Control group: 2F/7M. Sex is included as a covariate in the linear model design matrix.
2. **Proper multiple-testing correction** — The original 1,249 DEGs at "FDR < 0.05" reflected nominal p-values mislabelled as FDR-adjusted. Corrected count: 176 DEGs at nominal p < 0.01, |log₂FC| > 1.
3. **Nested cross-validation** — Fully nested LOOCV eliminates data leakage. Corrected AUC = 0.750 (vs. 0.991 from leaky pipeline).

## Repository Structure

```
INSILI-D-26-00281/
├── README.md
├── LICENSE
├── data/
│   ├── download_gse25186.R        # Download and preprocess GSE25186
│   └── sample_metadata.csv        # Sample annotations (sex, group)
├── R/
│   ├── 00_install_packages.R      # Install all required R packages
│   ├── 01_preprocessing_qc.R      # Data loading, normalisation, QC
│   ├── 02_differential_expression.R  # limma DE with sex covariate
│   ├── 03_wgcna_exploratory.R     # WGCNA (exploratory, n=15 < 20)
│   ├── 04_network_hub_analysis.R  # Hub gene identification
│   ├── 05_pathway_enrichment.R    # GO, KEGG, Reactome via Enrichr
│   ├── 06_tf_enrichment.R         # TF enrichment (ChEA, TRRUST)
│   ├── 07_nested_loocv.R          # Nested LOOCV (RF+LASSO+RF)
│   ├── 08_cross_cohort_validation.R  # GSE272189 transfer test
│   └── 09_figure_generation.R     # Generate all manuscript figures
├── python/
│   ├── requirements.txt
│   └── lincs_l1000_connectivity.py  # LINCS L1000 via Enrichr API
├── docking/
│   ├── README_docking.md          # Docking protocol and parameters
│   ├── vina_config.txt            # AutoDock Vina configuration
│   └── prepare_docking.py         # Receptor/ligand preparation
├── results/
│   └── expected_outputs.md        # Expected numerical results
└── figures/
    └── .gitkeep
```

## Software Requirements

### R (v4.3.1 or later)

| Package       | Version  | Purpose                            |
|---------------|----------|------------------------------------|
| limma         | 3.56.2   | Differential expression            |
| WGCNA         | 1.72-5   | Co-expression network analysis     |
| glmnet        | 4.1-8    | LASSO-regularised logistic regression |
| randomForest  | 4.7-1.1  | Random Forest classification       |
| pROC          | 1.18.5   | ROC curve analysis                 |
| GEOquery      | 2.68.0   | Download GEO datasets              |
| enrichR       | 3.2      | Enrichr API access                 |
| ggplot2       | 3.4+     | Visualisation                      |
| pheatmap      | 1.0.12   | Heatmaps                          |
| org.Hs.eg.db  | 3.17.0   | Gene annotation                   |

### Python (3.9+)

| Package  | Purpose                     |
|----------|-----------------------------|
| requests | Enrichr REST API calls      |
| pandas   | Data manipulation           |
| json     | API response parsing        |

### External Software

| Tool             | Version | Purpose                |
|------------------|---------|------------------------|
| AutoDock Vina    | 1.2.5   | Molecular docking      |
| Open Babel       | 3.1.1   | Ligand format conversion |

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/<username>/INSILI-D-26-00281.git
cd INSILI-D-26-00281

# 2. Install R packages
Rscript R/00_install_packages.R

# 3. Download and preprocess GSE25186 data
Rscript data/download_gse25186.R

# 4. Run the complete analysis pipeline (in order)
Rscript R/01_preprocessing_qc.R
Rscript R/02_differential_expression.R
Rscript R/03_wgcna_exploratory.R
Rscript R/04_network_hub_analysis.R
Rscript R/05_pathway_enrichment.R
Rscript R/06_tf_enrichment.R
Rscript R/07_nested_loocv.R
Rscript R/08_cross_cohort_validation.R
Rscript R/09_figure_generation.R

# 5. LINCS L1000 connectivity analysis (Python)
pip install -r python/requirements.txt
python python/lincs_l1000_connectivity.py

# 6. Molecular docking (requires AutoDock Vina)
# See docking/README_docking.md for detailed instructions
```

## Key Parameters

All parameters are set to match the manuscript exactly:

| Analysis                | Parameter              | Value                |
|-------------------------|------------------------|----------------------|
| Differential expression | Method                 | limma eBayes (moderated t-test) |
|                         | Design matrix          | ~ group + sex        |
|                         | FDR correction         | Benjamini-Hochberg   |
|                         | DEG threshold          | nominal p < 0.01, \|log₂FC\| > 1 |
| WGCNA                   | Soft-threshold power   | β = 6                |
|                         | Top genes              | 3,000 by MAD         |
|                         | Correlation            | Pearson              |
|                         | Min module size        | 30                   |
|                         | Hub metric             | kTotal × Gene Significance |
| Nested LOOCV            | Outer folds            | 15 (LOOCV)           |
|                         | Feature pre-filter     | RF importance (500 trees, Gini) |
|                         | Top candidates         | 50                   |
|                         | Feature selection      | LASSO (α=1, 5-fold inner CV) |
|                         | Classifier             | Random Forest (500 trees) |
|                         | Permutations           | 1,000                |
|                         | Bootstrap resamples    | 1,000                |
| Molecular docking       | Software               | AutoDock Vina 1.2.5  |
|                         | Target                 | MED13L (AF-Q71F56-F1, v6) |
|                         | Domain                 | CDK8-interaction (residues 100–510) |
|                         | Grid                   | 40 × 40 × 40 Å      |
|                         | Spacing                | 0.375 Å              |
|                         | Exhaustiveness         | 32                   |
|                         | Modes                  | 10                   |
|                         | Energy range           | 4 kcal/mol           |
| LINCS L1000             | Library                | Chem Pert Consensus Sigs |
|                         | API                    | Enrichr REST API     |

## Expected Results

See `results/expected_outputs.md` for a complete table of expected numerical outputs matching the manuscript.

**Key headline numbers:**
- 176 DEGs at nominal p < 0.01, |log₂FC| > 1
- 0 genes at FDR < 0.05 (after BH correction)
- MED13L: log₂FC = −3.19, p = 1.63 × 10⁻³
- Nested LOOCV AUC = 0.750 (95% CI 0.43–1.00, permutation p = 0.062)
- Leaky pipeline AUC = 0.991 (leakage gap = 0.24)
- 8-gene stable signature: PDLIM3, LOC650406, C7orf29, DTX3L, CHRNA2, SERPINB4, CXCL9, C1orf187
- Curcumin: 76 LINCS hits, p = 8.71 × 10⁻³
- Resveratrol: 18 LINCS hits, p = 7.49 × 10⁻³
- NAC: NOT significant (p = 0.37)
- Dexamethasone docking: −6.2 kcal/mol; Resveratrol: −5.9; NAC: −3.7

## Data Sources

| Dataset    | Source     | Description                              |
|------------|-----------|------------------------------------------|
| GSE25186   | NCBI GEO  | 6 PCD + 9 control nasal epithelial brushings (Geremek et al. 2014) |
| GSE272189  | NCBI GEO  | 71,396 cells, 13 donors, DNAH5-mutant PCD scRNA-seq (Koenitzer et al. 2024) |
| AF-Q71F56  | AlphaFold | MED13L predicted structure (model v6)     |

## Citation

If you use this code, please cite:

> Jitender, Hossain MW, Mohanty S, Kateriya S. Integrative Transcriptomic Analysis Reveals Molecular Signatures and Candidate Therapeutic Targets in Primary Ciliary Dyskinesia. *In Silico Research in Biomedicine*. 2026. (INSILI-D-26-00281)

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Suneel Kateriya — Laboratory of Optobiotechnology, School of Biotechnology, Jawaharlal Nehru University, New Delhi, India.
