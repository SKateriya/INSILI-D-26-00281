##############################################################################
## Script:  01_preprocessing_qc.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Load raw expression data from GSE25186 (Illumina HumanHT-12 V3.0,
##          GPL6947), perform background correction, quantile normalisation,
##          log2 transformation, probe-to-gene collapse, sex verification,
##          PCA, and QC plots.
##
## Dataset: GSE25186 — 6 PCD samples vs 9 healthy controls
##          Platform: Illumina HumanHT-12 V3.0 (GPL6947)
##          Sex imbalance: PCD 4F/2M, Control 2F/7M
##
## Inputs:
##   data/GSE25186_expression_matrix.csv   — probe-level expression (probes x samples)
##   data/sample_metadata.csv              — columns: sample_id, group (PCD/Control), sex (F/M)
##   data/GPL6947_probe_annotation.csv     — columns: probe_id, gene_symbol, entrez_id, ...
##
## Outputs:
##   results/01_normalised_expression.RData — normalised gene-level matrix + metadata
##   results/01_normalised_expression.csv   — same, as flat CSV
##   results/01_sex_verification.csv        — XIST / RPS4Y1 expression vs metadata sex
##   figures/01_boxplot_distributions.pdf   — per-sample expression boxplots
##   figures/01_pca_plot.pdf                — PCA scatter (PC1 vs PC2)
##
## Expected PCA statistics (from manuscript):
##   PC1 explains 30.6% of variance
##   PC2 explains 23.3% of variance
##   ANOVA PC1 ~ disease: p = 2.1e-5
##
## Usage:   Rscript 01_preprocessing_qc.R
##############################################################################

cat("=== 01_preprocessing_qc.R ===\n")
cat("Starting preprocessing and QC for GSE25186 ...\n\n")

## ---- load libraries --------------------------------------------------------
suppressPackageStartupMessages({
    library(limma)
    library(ggplot2)
    library(RColorBrewer)
})

## ---- create output directories ---------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load data ----------------------------------------------------------
cat("[1/8] Loading input data ...\n")

# Expression matrix
expr_file <- "data/GSE25186_expression_matrix.csv"
if (!file.exists(expr_file)) {
    stop("Expression matrix not found: ", expr_file,
         "\n  Please place the GSE25186 expression data in data/")
}
expr_raw <- read.csv(expr_file, row.names = 1, check.names = FALSE)
cat(sprintf("  Expression matrix: %d probes x %d samples\n",
            nrow(expr_raw), ncol(expr_raw)))

# Sample metadata
meta_file <- "data/sample_metadata.csv"
if (!file.exists(meta_file)) {
    stop("Sample metadata not found: ", meta_file)
}
sample_meta <- read.csv(meta_file, stringsAsFactors = FALSE)
cat(sprintf("  Sample metadata: %d samples\n", nrow(sample_meta)))
cat(sprintf("  Groups: %s\n",
            paste(names(table(sample_meta$group)),
                  table(sample_meta$group), sep = "=", collapse = ", ")))
cat(sprintf("  Sex: %s\n",
            paste(names(table(sample_meta$sex)),
                  table(sample_meta$sex), sep = "=", collapse = ", ")))

# Probe annotation
annot_file <- "data/GPL6947_probe_annotation.csv"
if (!file.exists(annot_file)) {
    stop("Probe annotation not found: ", annot_file)
}
probe_annot <- read.csv(annot_file, stringsAsFactors = FALSE)
cat(sprintf("  Probe annotation: %d probes\n", nrow(probe_annot)))

## ---- 2. Align samples -----------------------------------------------------
cat("\n[2/8] Aligning sample order between expression matrix and metadata ...\n")

# Ensure metadata sample IDs match column names
common_samples <- intersect(colnames(expr_raw), sample_meta$sample_id)
if (length(common_samples) == 0) {
    stop("No common sample IDs between expression matrix columns and metadata$sample_id")
}
if (length(common_samples) < ncol(expr_raw)) {
    warning(sprintf("Only %d of %d samples matched — subsetting.",
                    length(common_samples), ncol(expr_raw)))
}

expr_raw   <- expr_raw[, common_samples, drop = FALSE]
sample_meta <- sample_meta[match(common_samples, sample_meta$sample_id), ]
rownames(sample_meta) <- sample_meta$sample_id

cat(sprintf("  Aligned: %d samples (PCD=%d, Control=%d)\n",
            length(common_samples),
            sum(sample_meta$group == "PCD"),
            sum(sample_meta$group == "Control")))

# Verify expected sample counts
stopifnot(sum(sample_meta$group == "PCD") == 6)
stopifnot(sum(sample_meta$group == "Control") == 9)

## ---- 3. Background correction ---------------------------------------------
cat("\n[3/8] Background correction (normexp method) ...\n")
expr_mat <- as.matrix(expr_raw)

# limma::backgroundCorrect expects an EListRaw or matrix
# For Illumina data, normexp offset is standard
expr_bg <- backgroundCorrect(expr_mat, method = "normexp")
# backgroundCorrect returns an EList when given a matrix; extract $E
if (is(expr_bg, "EList") || is(expr_bg, "EListRaw")) {
    expr_bg <- expr_bg$E
}
cat("  Background correction complete.\n")

## ---- 4. Log2 transformation -----------------------------------------------
cat("\n[4/8] Log2 transformation (if needed) ...\n")

# Check if data is already on log scale by examining range
data_range <- range(expr_bg, na.rm = TRUE)
cat(sprintf("  Data range after background correction: [%.2f, %.2f]\n",
            data_range[1], data_range[2]))

if (data_range[2] > 100) {
    cat("  Data appears to be on raw intensity scale — applying log2.\n")
    # Ensure all values are positive before log
    expr_bg[expr_bg <= 0] <- min(expr_bg[expr_bg > 0], na.rm = TRUE) / 2
    expr_log <- log2(expr_bg)
} else {
    cat("  Data appears to be already log-transformed — skipping.\n")
    expr_log <- expr_bg
}

## ---- 5. Quantile normalisation ---------------------------------------------
cat("\n[5/8] Quantile normalisation (limma::normalizeBetweenArrays) ...\n")
expr_norm <- normalizeBetweenArrays(expr_log, method = "quantile")
cat(sprintf("  Normalised matrix: %d probes x %d samples\n",
            nrow(expr_norm), ncol(expr_norm)))

## ---- 6. Probe-to-gene mapping ----------------------------------------------
cat("\n[6/8] Probe-to-gene collapse (keep probe with highest mean expression) ...\n")

# Match probe IDs
common_probes <- intersect(rownames(expr_norm), probe_annot$probe_id)
if (length(common_probes) == 0) {
    # Try matching without explicit probe_id column — first column might be ID
    if ("ID" %in% colnames(probe_annot)) {
        probe_annot$probe_id <- probe_annot$ID
        common_probes <- intersect(rownames(expr_norm), probe_annot$probe_id)
    }
}
cat(sprintf("  Probes with annotation: %d of %d\n",
            length(common_probes), nrow(expr_norm)))

if (length(common_probes) == 0) {
    stop("No probe IDs matched between expression matrix and annotation file.")
}

# Subset to annotated probes
expr_annot <- expr_norm[common_probes, , drop = FALSE]
annot_sub  <- probe_annot[match(common_probes, probe_annot$probe_id), ]

# Identify gene_symbol column (handle variations)
gene_col <- intersect(c("gene_symbol", "Gene.Symbol", "Symbol", "GENE_SYMBOL",
                         "Gene_Symbol", "gene_name"), colnames(annot_sub))
if (length(gene_col) == 0) {
    stop("Cannot find gene symbol column in probe annotation.")
}
gene_col <- gene_col[1]

# Remove probes with missing or empty gene symbols
valid <- !is.na(annot_sub[[gene_col]]) & annot_sub[[gene_col]] != "" &
         annot_sub[[gene_col]] != "---"
expr_annot <- expr_annot[valid, , drop = FALSE]
annot_sub  <- annot_sub[valid, ]
cat(sprintf("  Probes with valid gene symbol: %d\n", nrow(expr_annot)))

# For each gene, keep the probe with the highest mean expression
probe_means <- rowMeans(expr_annot, na.rm = TRUE)
gene_symbols <- annot_sub[[gene_col]]

# Handle multi-mapped probes (e.g., "GENE1 /// GENE2") — take first symbol
gene_symbols <- gsub(" ///.*", "", gene_symbols)

# Split by gene, pick max-mean probe
gene_probe_df <- data.frame(
    probe    = rownames(expr_annot),
    gene     = gene_symbols,
    mean_exp = probe_means,
    stringsAsFactors = FALSE
)

best_probes <- tapply(seq_len(nrow(gene_probe_df)), gene_probe_df$gene,
                      function(idx) {
                          idx[which.max(gene_probe_df$mean_exp[idx])]
                      })
best_probes <- unlist(best_probes)

expr_gene <- expr_annot[gene_probe_df$probe[best_probes], , drop = FALSE]
rownames(expr_gene) <- gene_probe_df$gene[best_probes]

cat(sprintf("  Gene-level matrix: %d genes x %d samples\n",
            nrow(expr_gene), ncol(expr_gene)))

## ---- 7. Sex verification ---------------------------------------------------
cat("\n[7/8] Sex verification using XIST and RPS4Y1 expression ...\n")

sex_genes <- c("XIST", "RPS4Y1")
sex_present <- sex_genes[sex_genes %in% rownames(expr_gene)]

if (length(sex_present) > 0) {
    sex_expr <- data.frame(
        sample  = colnames(expr_gene),
        group   = sample_meta$group,
        sex_meta = sample_meta$sex,
        stringsAsFactors = FALSE
    )
    for (g in sex_present) {
        sex_expr[[g]] <- as.numeric(expr_gene[g, ])
    }

    # Infer sex from expression
    if ("XIST" %in% sex_present && "RPS4Y1" %in% sex_present) {
        # Female: high XIST, low RPS4Y1; Male: low XIST, high RPS4Y1
        xist_med  <- median(sex_expr$XIST)
        rps4y1_med <- median(sex_expr$RPS4Y1)
        sex_expr$sex_inferred <- ifelse(
            sex_expr$XIST > xist_med & sex_expr$RPS4Y1 < rps4y1_med, "F",
            ifelse(sex_expr$XIST < xist_med & sex_expr$RPS4Y1 > rps4y1_med, "M",
                   "ambiguous"))
    } else if ("XIST" %in% sex_present) {
        xist_med <- median(sex_expr$XIST)
        sex_expr$sex_inferred <- ifelse(sex_expr$XIST > xist_med, "F", "M")
    } else {
        rps4y1_med <- median(sex_expr$RPS4Y1)
        sex_expr$sex_inferred <- ifelse(sex_expr$RPS4Y1 > rps4y1_med, "M", "F")
    }

    # Compare
    concordance <- sum(sex_expr$sex_meta == sex_expr$sex_inferred, na.rm = TRUE)
    cat(sprintf("  Sex concordance (metadata vs inferred): %d / %d samples\n",
                concordance, nrow(sex_expr)))
    if (concordance < nrow(sex_expr)) {
        mismatches <- sex_expr[sex_expr$sex_meta != sex_expr$sex_inferred, ]
        cat("  WARNING: sex mismatches detected:\n")
        print(mismatches)
    }

    # Verify sex imbalance: PCD 4F/2M, Control 2F/7M
    pcd_sex  <- table(sample_meta$sex[sample_meta$group == "PCD"])
    ctrl_sex <- table(sample_meta$sex[sample_meta$group == "Control"])
    cat(sprintf("  PCD sex distribution:     F=%d, M=%d (expected 4F/2M)\n",
                pcd_sex["F"], pcd_sex["M"]))
    cat(sprintf("  Control sex distribution: F=%d, M=%d (expected 2F/7M)\n",
                ctrl_sex["F"], ctrl_sex["M"]))

    # Save sex verification
    write.csv(sex_expr, "results/01_sex_verification.csv", row.names = FALSE)
    cat("  Saved: results/01_sex_verification.csv\n")
} else {
    cat("  WARNING: Neither XIST nor RPS4Y1 found in gene-level matrix.\n")
    cat("  Sex verification cannot be performed.\n")
}

## ---- 8. PCA ----------------------------------------------------------------
cat("\n[8/8] PCA analysis ...\n")

# Centre and scale genes (rows = genes, columns = samples → transpose for prcomp)
pca_result <- prcomp(t(expr_gene), center = TRUE, scale. = FALSE)

# Variance explained
var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
cat(sprintf("  PC1 explains %.1f%% of variance (manuscript: 30.6%%)\n", var_explained[1]))
cat(sprintf("  PC2 explains %.1f%% of variance (manuscript: 23.3%%)\n", var_explained[2]))

# ANOVA: PC1 ~ disease group
pc_scores <- data.frame(
    PC1   = pca_result$x[, 1],
    PC2   = pca_result$x[, 2],
    group = sample_meta$group,
    sex   = sample_meta$sex,
    sample = sample_meta$sample_id,
    stringsAsFactors = FALSE
)
anova_pc1 <- anova(lm(PC1 ~ group, data = pc_scores))
pval_pc1  <- anova_pc1$`Pr(>F)`[1]
cat(sprintf("  ANOVA PC1 ~ group: p = %.2e (manuscript: 2.1e-5)\n", pval_pc1))

## ---- QC plots --------------------------------------------------------------
cat("\nGenerating QC plots ...\n")

# 1. Boxplot of expression distributions
pdf("figures/01_boxplot_distributions.pdf", width = 10, height = 6)
par(mar = c(8, 4, 3, 1))
group_cols <- ifelse(sample_meta$group == "PCD",
                     brewer.pal(3, "Set1")[1],   # red
                     brewer.pal(3, "Set1")[2])    # blue
boxplot(expr_gene,
        las = 2,
        col = group_cols,
        main = "GSE25186 — Normalised Expression Distributions",
        ylab = "log2 expression",
        cex.axis = 0.7)
legend("topright",
       legend = c("PCD (n=6)", "Control (n=9)"),
       fill = brewer.pal(3, "Set1")[1:2],
       border = NA, bty = "n")
dev.off()
cat("  Saved: figures/01_boxplot_distributions.pdf\n")

# 2. PCA plot
pca_df <- data.frame(
    PC1   = pca_result$x[, 1],
    PC2   = pca_result$x[, 2],
    Group = sample_meta$group,
    Sex   = sample_meta$sex
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group, shape = Sex)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_colour_manual(values = c("PCD" = "#E41A1C", "Control" = "#377EB8")) +
    labs(
        title = "GSE25186 — PCA of Normalised Expression",
        x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
        y = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
        caption = sprintf("ANOVA PC1 ~ group: p = %.2e", pval_pc1)
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
    )

ggsave("figures/01_pca_plot.pdf", pca_plot, width = 8, height = 6)
cat("  Saved: figures/01_pca_plot.pdf\n")

## ---- Save outputs ----------------------------------------------------------
cat("\nSaving normalised expression data ...\n")

# RData
save(expr_gene, sample_meta, pca_result, var_explained, pval_pc1,
     file = "results/01_normalised_expression.RData")
cat("  Saved: results/01_normalised_expression.RData\n")

# CSV
expr_out <- data.frame(gene = rownames(expr_gene), expr_gene,
                       check.names = FALSE, stringsAsFactors = FALSE)
write.csv(expr_out, "results/01_normalised_expression.csv", row.names = FALSE)
cat("  Saved: results/01_normalised_expression.csv\n")

## ---- Summary ---------------------------------------------------------------
cat("\n--- Preprocessing Summary ---\n")
cat(sprintf("  Dataset:           GSE25186\n"))
cat(sprintf("  Platform:          Illumina HumanHT-12 V3.0 (GPL6947)\n"))
cat(sprintf("  Samples:           %d (PCD=%d, Control=%d)\n",
            ncol(expr_gene),
            sum(sample_meta$group == "PCD"),
            sum(sample_meta$group == "Control")))
cat(sprintf("  Genes:             %d\n", nrow(expr_gene)))
cat(sprintf("  Normalisation:     quantile (limma)\n"))
cat(sprintf("  Background corr:   normexp (limma)\n"))
cat(sprintf("  PC1:               %.1f%% variance\n", var_explained[1]))
cat(sprintf("  PC2:               %.1f%% variance\n", var_explained[2]))
cat(sprintf("  PC1 ~ group:       p = %.2e\n", pval_pc1))

cat("\n=== 01_preprocessing_qc.R complete ===\n")
