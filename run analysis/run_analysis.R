# ============================================================
# PCD Differential Expression Analysis — GSE25186
# Platform: Illumina HumanHT-12 v4 microarray
#
# This script performs TWO analyses:
#   1. PRIMARY   — All 15 samples (6 PCD vs 9 Control),
#                  sex-corrected limma
#   2. SENSITIVITY — 13 samples (5 PCD vs 8 Control),
#                    2 QC outliers excluded, sex-corrected limma
#
# The sensitivity analysis tests robustness of key findings
# after removing samples with ambiguous cilia marker profiles.
#
# RUN:  Rscript run_analysis.R
# ============================================================

# --- Install limma if needed --------------------------------
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("limma", ask = FALSE)
}
library(limma)

# --- Load expression data -----------------------------------
expr <- read.csv("data/GSE25186_gene_expression.csv",
                 row.names = 1, check.names = FALSE)
expr_log <- log2(pmax(as.matrix(expr), 1))

# --- Key genes of interest ----------------------------------
# Note: GSTT2B is annotated as GSTT2 on HumanHT-12 v4
# Note: METTL26 is annotated as C1orf194 on HumanHT-12 v4
key_genes <- c("MED13L", "DDX58", "HLA-C", "CXCL9",
               "DNAH5", "DNAH11", "FRMPD3",
               "GSTT2", "C1orf194")


# =============================================================
# FUNCTION: Run DE analysis
# =============================================================
run_de <- function(meta, expr_log, label) {

  cat("\n")
  cat("############################################################\n")
  cat(sprintf("  %s\n", label))
  cat("############################################################\n\n")

  cat(sprintf("Samples: %d\n", nrow(meta)))
  cat(sprintf("  PCD     (n=%d): %s\n",
      sum(meta$Group == "PCD"),
      paste(meta$GSM_ID[meta$Group == "PCD"], collapse = ", ")))
  cat(sprintf("  Control (n=%d): %s\n",
      sum(meta$Group == "Control"),
      paste(meta$GSM_ID[meta$Group == "Control"], collapse = ", ")))
  cat(sprintf("  Sex: %d F, %d M\n",
      sum(meta$Sex == "F"), sum(meta$Sex == "M")))

  # Design matrix with sex covariate
  grp <- factor(meta$Group, levels = c("Control", "PCD"))
  sex <- factor(meta$Sex, levels = c("F", "M"))
  design <- model.matrix(~ sex + grp)

  # Fit limma
  fit <- eBayes(lmFit(expr_log[, meta$GSM_ID], design))
  res <- topTable(fit, coef = "grpPCD", number = Inf, sort.by = "none")

  # Key gene results
  cat("\n------------------------------------------------------------\n")
  cat("KEY GENE RESULTS\n")
  cat("------------------------------------------------------------\n\n")
  cat(sprintf("%-12s %10s %12s %12s %6s\n",
              "Gene", "logFC", "P.Value", "adj.P.Val", "Sig"))
  cat(paste(rep("-", 56), collapse = ""), "\n")

  for (g in key_genes) {
    if (g %in% rownames(res)) {
      r   <- res[g, ]
      sig <- ifelse(r$adj.P.Val < 0.05, "***",
             ifelse(r$adj.P.Val < 0.10, "**",
             ifelse(r$P.Value  < 0.05, "*", "")))
      cat(sprintf("%-12s %10.3f %12.6f %12.6f %6s\n",
                  g, r$logFC, r$P.Value, r$adj.P.Val, sig))
    }
  }

  # DEG counts
  n_fdr    <- sum(res$adj.P.Val < 0.05)
  n_fdr_fc <- sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1)
  n_nom    <- sum(res$P.Value < 0.05)
  n_nom_fc <- sum(res$P.Value < 0.05 & abs(res$logFC) > 1)

  cat("\n------------------------------------------------------------\n")
  cat("DIFFERENTIALLY EXPRESSED GENES\n")
  cat("------------------------------------------------------------\n\n")
  cat(sprintf("  FDR < 0.05                     : %d\n", n_fdr))
  cat(sprintf("  FDR < 0.05 & |log2FC| > 1      : %d\n", n_fdr_fc))
  cat(sprintf("  Nominal p < 0.05               : %d\n", n_nom))
  cat(sprintf("  Nominal p < 0.05 & |log2FC| > 1: %d\n", n_nom_fc))

  # Top 20 by adj.P.Val
  cat("\n------------------------------------------------------------\n")
  cat("TOP 20 DEGs (by adjusted p-value)\n")
  cat("------------------------------------------------------------\n\n")
  top20 <- head(res[order(res$adj.P.Val), ], 20)
  cat(sprintf("%-20s %10s %12s %12s\n",
              "Gene", "logFC", "P.Value", "adj.P.Val"))
  cat(paste(rep("-", 58), collapse = ""), "\n")
  for (i in seq_len(nrow(top20))) {
    cat(sprintf("%-20s %10.3f %12.6f %12.6f\n",
                rownames(top20)[i],
                top20$logFC[i], top20$P.Value[i], top20$adj.P.Val[i]))
  }

  # Top 20 by |logFC| among nominal significant
  cat("\n------------------------------------------------------------\n")
  cat("TOP 20 by |logFC| (nominal p < 0.05)\n")
  cat("------------------------------------------------------------\n\n")
  sig_nom <- res[res$P.Value < 0.05, ]
  top_fc  <- head(sig_nom[order(-abs(sig_nom$logFC)), ], 20)
  cat(sprintf("%-20s %10s %12s %12s\n",
              "Gene", "logFC", "P.Value", "adj.P.Val"))
  cat(paste(rep("-", 58), collapse = ""), "\n")
  for (i in seq_len(nrow(top_fc))) {
    cat(sprintf("%-20s %10.3f %12.6f %12.6f\n",
                rownames(top_fc)[i],
                top_fc$logFC[i], top_fc$P.Value[i], top_fc$adj.P.Val[i]))
  }

  return(res)
}


# =============================================================
# ANALYSIS 1: PRIMARY (6 PCD vs 9 Control — all 15 samples)
# =============================================================
meta_primary <- read.csv("data/sample_metadata_primary.csv")
res_primary  <- run_de(meta_primary, expr_log,
                       "ANALYSIS 1: PRIMARY (6 PCD vs 9 Control, sex-corrected)")

# Save primary results
res_primary$Gene <- rownames(res_primary)
write.csv(res_primary[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")],
          "results_primary_6v9.csv", row.names = FALSE)


# =============================================================
# ANALYSIS 2: SENSITIVITY (5 PCD vs 8 Control — QC-filtered)
# =============================================================
meta_sens <- read.csv("data/sample_metadata_sensitivity.csv")
res_sens  <- run_de(meta_sens, expr_log,
                    "ANALYSIS 2: SENSITIVITY (5 PCD vs 8 Control, QC-filtered, sex-corrected)")

# Save sensitivity results
res_sens$Gene <- rownames(res_sens)
write.csv(res_sens[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")],
          "results_sensitivity_5v8.csv", row.names = FALSE)


# =============================================================
# COMPARISON: Primary vs Sensitivity
# =============================================================
cat("\n\n")
cat("############################################################\n")
cat("  COMPARISON: Primary (6v9) vs Sensitivity (5v8)\n")
cat("############################################################\n\n")

cat(sprintf("%-12s %10s %10s %10s %8s\n",
            "Gene", "Primary", "Sensitiv.", "Diff", "Agree%"))
cat(paste(rep("-", 54), collapse = ""), "\n")

total_agree <- 0
n_genes     <- 0

for (g in key_genes) {
  if (g %in% rownames(res_primary) && g %in% rownames(res_sens)) {
    lfc_p <- res_primary[g, "logFC"]
    lfc_s <- res_sens[g, "logFC"]
    diff  <- lfc_s - lfc_p
    agree <- (1 - abs(diff) / abs(lfc_p)) * 100
    total_agree <- total_agree + agree
    n_genes     <- n_genes + 1
    cat(sprintf("%-12s %10.3f %10.3f %+10.3f %7.1f%%\n",
                g, lfc_p, lfc_s, diff, agree))
  }
}

cat(paste(rep("-", 54), collapse = ""), "\n")
cat(sprintf("%-12s %10s %10s %10s %7.1f%%\n",
            "AVERAGE", "", "", "", total_agree / n_genes))

# DEG comparison
n_fdr_p <- sum(res_primary$adj.P.Val < 0.05 & abs(res_primary$logFC) > 1)
n_fdr_s <- sum(res_sens$adj.P.Val < 0.05 & abs(res_sens$logFC) > 1)
cat(sprintf("\nDEGs (FDR<0.05, |logFC|>1):  Primary = %d,  Sensitivity = %d\n",
            n_fdr_p, n_fdr_s))

# Direction consistency
common <- intersect(rownames(res_primary), rownames(res_sens))
same_dir <- sum(sign(res_primary[common, "logFC"]) == sign(res_sens[common, "logFC"]))
cat(sprintf("Direction consistency:        %d / %d genes (%.1f%%)\n",
            same_dir, length(common),
            100 * same_dir / length(common)))

cat("\n############################################################\n")
cat("  Analysis complete. Output files:\n")
cat("    results_primary_6v9.csv\n")
cat("    results_sensitivity_5v8.csv\n")
cat("############################################################\n")