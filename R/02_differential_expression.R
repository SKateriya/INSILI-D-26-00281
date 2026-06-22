##############################################################################
## Script:  02_differential_expression.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Perform differential expression analysis on the normalised,
##          gene-level expression matrix produced by 01_preprocessing_qc.R.
##          Uses empirical Bayes moderated t-test (limma) with sex as a
##          covariate in the design matrix and BH FDR correction.
##
## Dataset: GSE25186 — 6 PCD vs 9 controls
##          Platform: Illumina HumanHT-12 V3.0, GPL6947
##          Sex imbalance: PCD 4F/2M, Control 2F/7M
##
## Inputs:
##   results/01_normalised_expression.RData — from 01_preprocessing_qc.R
##       Contains: expr_gene (gene x sample matrix), sample_meta (metadata)
##
## Outputs:
##   results/02_DE_full_results.csv         — all genes, full topTable output
##   results/02_DEG_176_p01_lfc1.csv        — 176 DEGs (p < 0.01, |log2FC| > 1)
##   results/02_DEG_1014_p05_lfc1.csv       — 1,014 DEGs (p < 0.05, |log2FC| > 1)
##   results/02_volcano_data.csv            — data for volcano plot
##   figures/02_volcano_plot.pdf             — volcano plot
##   results/02_DE_results.RData            — all DE objects for downstream scripts
##
## Expected results (from manuscript):
##   176 DEGs at nominal p < 0.01 with |log2FC| > 1
##   1,014 DEGs at nominal p < 0.05 with |log2FC| > 1
##   NO gene survives FDR < 0.05 (all adj.P.Val = 1.00)
##   Top gene: MED13L logFC = -3.19, p = 1.63e-3
##
## Usage:   Rscript 02_differential_expression.R
##############################################################################

cat("=== 02_differential_expression.R ===\n")
cat("Starting differential expression analysis for GSE25186 ...\n\n")

set.seed(42)

## ---- load libraries --------------------------------------------------------
suppressPackageStartupMessages({
    library(limma)
    library(ggplot2)
    library(RColorBrewer)
})

## ---- create output directories ---------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load normalised data -----------------------------------------------
cat("[1/7] Loading normalised expression data ...\n")

rdata_file <- "results/01_normalised_expression.RData"
if (!file.exists(rdata_file)) {
    stop("Normalised data not found: ", rdata_file,
         "\n  Please run 01_preprocessing_qc.R first.")
}
load(rdata_file)  # loads expr_gene, sample_meta, pca_result, var_explained, pval_pc1

cat(sprintf("  Expression matrix: %d genes x %d samples\n",
            nrow(expr_gene), ncol(expr_gene)))
cat(sprintf("  Groups: PCD=%d, Control=%d\n",
            sum(sample_meta$group == "PCD"),
            sum(sample_meta$group == "Control")))
cat(sprintf("  Sex:    F=%d, M=%d\n",
            sum(sample_meta$sex == "F"),
            sum(sample_meta$sex == "M")))

## ---- 2. Build design matrix -----------------------------------------------
cat("\n[2/7] Building design matrix: ~ group + sex ...\n")

# Set factor levels — Control is the reference
sample_meta$group <- factor(sample_meta$group, levels = c("Control", "PCD"))
sample_meta$sex   <- factor(sample_meta$sex,   levels = c("M", "F"))

design <- model.matrix(~ group + sex, data = sample_meta)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- gsub("sex",   "", colnames(design))

cat("  Design matrix:\n")
print(design)
cat(sprintf("\n  Coefficients: %s\n", paste(colnames(design), collapse = ", ")))
cat("  The 'PCD' coefficient captures the PCD vs Control effect, adjusted for sex.\n")

## ---- 3. Fit linear model ---------------------------------------------------
cat("\n[3/7] Fitting linear model (limma::lmFit) ...\n")

fit <- lmFit(expr_gene, design)
cat("  Linear model fitted.\n")

## ---- 4. Empirical Bayes moderation -----------------------------------------
cat("\n[4/7] Applying empirical Bayes moderation (limma::eBayes) ...\n")

fit_eb <- eBayes(fit)
cat(sprintf("  Prior df:       %.2f\n", fit_eb$df.prior))
cat(sprintf("  Prior variance: %.4f\n", fit_eb$s2.prior))

## ---- 5. Extract results ----------------------------------------------------
cat("\n[5/7] Extracting DE results (BH FDR correction) ...\n")

# Extract for the PCD coefficient
de_results <- topTable(fit_eb, coef = "PCD", adjust.method = "BH",
                       number = Inf, sort.by = "none")

# Add gene name column
de_results$gene <- rownames(de_results)

# Reorder by p-value
de_results <- de_results[order(de_results$P.Value), ]

cat(sprintf("  Total genes tested: %d\n", nrow(de_results)))
cat(sprintf("  Min nominal p-value: %.2e\n", min(de_results$P.Value)))
cat(sprintf("  Min adjusted p-value: %.2f\n", min(de_results$adj.P.Val)))

## ---- 6. Define DEG sets ----------------------------------------------------
cat("\n[6/7] Defining DEG sets ...\n")

# Strict threshold: nominal p < 0.01 AND |log2FC| > 1
deg_strict <- de_results[de_results$P.Value < 0.01 &
                         abs(de_results$logFC) > 1, ]
cat(sprintf("  DEGs (p < 0.01, |log2FC| > 1): %d (manuscript expects 176)\n",
            nrow(deg_strict)))

# Broader threshold: nominal p < 0.05 AND |log2FC| > 1
deg_broad <- de_results[de_results$P.Value < 0.05 &
                        abs(de_results$logFC) > 1, ]
cat(sprintf("  DEGs (p < 0.05, |log2FC| > 1): %d (manuscript expects 1,014)\n",
            nrow(deg_broad)))

# FDR check
n_fdr <- sum(de_results$adj.P.Val < 0.05)
cat(sprintf("  DEGs at FDR < 0.05: %d (manuscript expects 0)\n", n_fdr))

if (n_fdr == 0) {
    cat("  CONFIRMED: No gene survives BH FDR < 0.05 correction.\n")
    cat("  All adjusted p-values = 1.00 as expected.\n")
} else {
    cat("  NOTE: Some genes survive FDR — differs from manuscript expectation.\n")
}

# Direction breakdown
n_up_strict   <- sum(deg_strict$logFC > 0)
n_down_strict <- sum(deg_strict$logFC < 0)
cat(sprintf("\n  Strict DEGs: %d up-regulated, %d down-regulated in PCD\n",
            n_up_strict, n_down_strict))

n_up_broad   <- sum(deg_broad$logFC > 0)
n_down_broad <- sum(deg_broad$logFC < 0)
cat(sprintf("  Broad DEGs:  %d up-regulated, %d down-regulated in PCD\n",
            n_up_broad, n_down_broad))

## ---- Top genes report ------------------------------------------------------
cat("\n--- Top genes (by nominal p-value) ---\n")
top_n <- min(20, nrow(de_results))
top_genes <- de_results[1:top_n, c("gene", "logFC", "AveExpr", "t",
                                    "P.Value", "adj.P.Val")]
print(top_genes, digits = 4, row.names = FALSE)

## ---- Manuscript gene checks ------------------------------------------------
cat("\n--- Manuscript key gene verification ---\n")

check_gene <- function(gene_name, expected_lfc, expected_p) {
    if (gene_name %in% de_results$gene) {
        row <- de_results[de_results$gene == gene_name, ]
        cat(sprintf("  %-10s  logFC = %+.2f (expect %+.2f)  p = %.2e (expect %.2e)  adj.P = %.2f\n",
                    gene_name, row$logFC, expected_lfc, row$P.Value, expected_p, row$adj.P.Val))
    } else {
        cat(sprintf("  %-10s  NOT FOUND in gene-level matrix\n", gene_name))
    }
}

cat("  Top DEGs from manuscript:\n")
check_gene("MED13L",  -3.19, 1.63e-3)
check_gene("LOC728452", -3.22, 1.08e-3)
check_gene("DDX58",   -3.15, NA)
check_gene("HLA-C",   -3.04, NA)
check_gene("CLN8",    -3.10, NA)
check_gene("TP53BP2", +3.20, NA)
check_gene("SLC18A1", +3.18, NA)

cat("\n  Sex-linked gene check (should be non-significant after adjustment):\n")
check_gene("XIST",    +0.55, 0.41)
check_gene("RPS4Y1",  -0.15, 0.77)
check_gene("EIF1AY",  -0.06, 0.93)

# Verify sex adjustment worked
if ("XIST" %in% de_results$gene) {
    xist_row <- de_results[de_results$gene == "XIST", ]
    if (xist_row$P.Value > 0.05) {
        cat("  CONFIRMED: XIST is non-significant (p > 0.05) — sex covariate is working.\n")
    } else {
        cat("  WARNING: XIST is nominally significant — sex adjustment may be incomplete.\n")
    }
}
if ("RPS4Y1" %in% de_results$gene) {
    rps4y1_row <- de_results[de_results$gene == "RPS4Y1", ]
    if (rps4y1_row$P.Value > 0.05) {
        cat("  CONFIRMED: RPS4Y1 is non-significant (p > 0.05) — sex covariate is working.\n")
    } else {
        cat("  WARNING: RPS4Y1 is nominally significant — sex adjustment may be incomplete.\n")
    }
}

## ---- 7. Save outputs -------------------------------------------------------
cat("\n[7/7] Saving results ...\n")

# Full DE results
write.csv(de_results, "results/02_DE_full_results.csv", row.names = FALSE)
cat("  Saved: results/02_DE_full_results.csv\n")

# 176-gene strict list
write.csv(deg_strict, "results/02_DEG_176_p01_lfc1.csv", row.names = FALSE)
cat(sprintf("  Saved: results/02_DEG_176_p01_lfc1.csv  (%d genes)\n", nrow(deg_strict)))

# 1,014-gene broad list
write.csv(deg_broad, "results/02_DEG_1014_p05_lfc1.csv", row.names = FALSE)
cat(sprintf("  Saved: results/02_DEG_1014_p05_lfc1.csv  (%d genes)\n", nrow(deg_broad)))

# Volcano plot data
volcano_data <- data.frame(
    gene      = de_results$gene,
    logFC     = de_results$logFC,
    neg_log10_p = -log10(de_results$P.Value),
    P.Value   = de_results$P.Value,
    adj.P.Val = de_results$adj.P.Val,
    significant_strict = de_results$P.Value < 0.01 & abs(de_results$logFC) > 1,
    significant_broad  = de_results$P.Value < 0.05 & abs(de_results$logFC) > 1,
    stringsAsFactors = FALSE
)
write.csv(volcano_data, "results/02_volcano_data.csv", row.names = FALSE)
cat("  Saved: results/02_volcano_data.csv\n")

# RData for downstream
save(de_results, deg_strict, deg_broad, fit, fit_eb, design,
     file = "results/02_DE_results.RData")
cat("  Saved: results/02_DE_results.RData\n")

## ---- Volcano plot ----------------------------------------------------------
cat("\nGenerating volcano plot ...\n")

# Assign categories for colouring
volcano_data$category <- "Not significant"
volcano_data$category[volcano_data$significant_broad &
                      !volcano_data$significant_strict] <- "p < 0.05 & |log2FC| > 1"
volcano_data$category[volcano_data$significant_strict] <- "p < 0.01 & |log2FC| > 1"

volcano_data$category <- factor(volcano_data$category,
    levels = c("Not significant",
               "p < 0.05 & |log2FC| > 1",
               "p < 0.01 & |log2FC| > 1"))

# Label top genes
top_labels <- c("MED13L", "DDX58", "HLA-C", "CLN8", "TP53BP2", "SLC18A1")
volcano_data$label <- ifelse(volcano_data$gene %in% top_labels,
                             volcano_data$gene, NA)

volcano_plot <- ggplot(volcano_data,
                       aes(x = logFC, y = neg_log10_p, colour = category)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_colour_manual(
        values = c("Not significant"            = "grey70",
                   "p < 0.05 & |log2FC| > 1"    = "#FDB863",
                   "p < 0.01 & |log2FC| > 1"    = "#D73027"),
        name = "Significance"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey60") +
    labs(
        title = "GSE25186 — Differential Expression: PCD vs Control",
        subtitle = sprintf(
            "Design: ~ group + sex | %d DEGs (p<0.01, |lfc|>1) | No gene at FDR<0.05",
            nrow(deg_strict)),
        x = expression(log[2]~fold~change~(PCD~vs~Control)),
        y = expression(-log[10]~(nominal~italic(p)))
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        legend.position = "bottom"
    )

# Add labels for key genes (if ggrepel available, use it; otherwise annotate)
if (requireNamespace("ggrepel", quietly = TRUE)) {
    volcano_plot <- volcano_plot +
        ggrepel::geom_text_repel(
            aes(label = label),
            size = 3, max.overlaps = 20,
            segment.color = "grey50",
            na.rm = TRUE, show.legend = FALSE
        )
} else {
    label_data <- volcano_data[!is.na(volcano_data$label), ]
    if (nrow(label_data) > 0) {
        volcano_plot <- volcano_plot +
            geom_text(
                data = label_data,
                aes(label = label),
                size = 3, vjust = -0.8, hjust = 0.5,
                show.legend = FALSE
            )
    }
}

ggsave("figures/02_volcano_plot.pdf", volcano_plot, width = 9, height = 7)
cat("  Saved: figures/02_volcano_plot.pdf\n")

## ---- Summary ---------------------------------------------------------------
cat("\n--- Differential Expression Summary ---\n")
cat(sprintf("  Design:              ~ group + sex\n"))
cat(sprintf("  Method:              empirical Bayes moderated t-test (limma)\n"))
cat(sprintf("  Multiple testing:    BH FDR correction\n"))
cat(sprintf("  Total genes tested:  %d\n", nrow(de_results)))
cat(sprintf("  DEGs (p<0.01, |lfc|>1):  %d  (manuscript: 176)\n", nrow(deg_strict)))
cat(sprintf("  DEGs (p<0.05, |lfc|>1):  %d  (manuscript: 1,014)\n", nrow(deg_broad)))
cat(sprintf("  DEGs at FDR < 0.05:      %d  (manuscript: 0)\n", n_fdr))
if ("MED13L" %in% de_results$gene) {
    med_row <- de_results[de_results$gene == "MED13L", ]
    cat(sprintf("  Top gene MED13L:     logFC = %+.2f (manuscript: -3.19)\n", med_row$logFC))
    cat(sprintf("                       p     = %.2e (manuscript: 1.63e-3)\n", med_row$P.Value))
}

cat("\n=== 02_differential_expression.R complete ===\n")
