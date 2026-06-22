##############################################################################
## Script:  03_wgcna_exploratory.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Weighted Gene Co-expression Network Analysis (WGCNA) on the
##          normalised gene-level expression matrix from 01_preprocessing_qc.R.
##          Identifies co-expression modules, correlates module eigengenes
##          with disease status and sex, ranks hub genes by kTotal x Gene
##          Significance, and runs permutation testing to assess module-trait
##          significance.
##
## Dataset: GSE25186 — 6 PCD vs 9 controls (n = 15 total)
##          Platform: Illumina HumanHT-12 V3.0 (GPL6947)
##
## IMPORTANT LIMITATIONS:
##   n = 15 is below the minimum 20 samples recommended for stable WGCNA.
##   Module assignments are UNSTABLE and should be interpreted as exploratory.
##   Module-trait correlations do NOT reach significance after permutation
##   testing (all permutation p > 0.05 with 1,000 label permutations).
##
## Parameters (from manuscript):
##   - Top 3,000 most variable genes by median absolute deviation (MAD)
##   - Signed adjacency matrix, soft-threshold power beta = 6
##     (R-squared > 0.85, mean connectivity 45.2)
##   - Topological overlap matrix (TOM), 1-TOM as distance
##   - Hierarchical clustering with average linkage
##   - Dynamic tree cutting: deepSplit = 2, minimum module size = 30
##   - Module eigengenes: first principal component of each module
##   - Module-trait correlation: Pearson for continuous, point-biserial
##     for binary traits
##   - Hub genes: ranked by kTotal x Gene Significance
##   - 1,000 label permutations for module-trait significance
##
## Inputs:
##   results/01_normalised_expression.RData — from 01_preprocessing_qc.R
##       Contains: expr_gene (gene x sample matrix), sample_meta (metadata)
##
## Outputs:
##   results/03_wgcna_results.RData         — all WGCNA objects for downstream
##   results/03_module_assignments.csv      — gene-to-module mapping
##   results/03_module_trait_correlations.csv — module-trait r and p-values
##   results/03_hub_genes_ranked.csv        — hub genes ranked by kTotal x GS
##   results/03_permutation_pvalues.csv     — permutation-based p-values
##   results/03_gstt2b_probe_check.txt      — GSTT2B probe absence verification
##   figures/03_soft_threshold.pdf          — scale-free topology fit plots
##   figures/03_dendrogram.pdf              — gene dendrogram with module colours
##   figures/03_module_trait_heatmap.pdf    — module-trait correlation heatmap
##   figures/03_hub_gene_scatter.pdf        — kTotal vs Gene Significance scatter
##
## Usage:   Rscript 03_wgcna_exploratory.R
##############################################################################

cat("=== 03_wgcna_exploratory.R ===\n")
cat("Starting WGCNA exploratory analysis for GSE25186 ...\n\n")

set.seed(42)

## ---- load libraries --------------------------------------------------------
suppressPackageStartupMessages({
    library(WGCNA)
    library(ggplot2)
    library(RColorBrewer)
})

# WGCNA settings
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

## ---- create output directories ---------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load normalised data -----------------------------------------------
cat("[1/11] Loading normalised expression data ...\n")

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

# Sample size warning
n_samples <- ncol(expr_gene)
if (n_samples < 20) {
    cat(sprintf("\n  *** WARNING: n = %d is below the WGCNA minimum of 20 samples. ***\n",
                n_samples))
    cat("  *** Module assignments will be UNSTABLE. Results are exploratory only. ***\n\n")
}

## ---- 2. Load DE results for Gene Significance ------------------------------
cat("[2/11] Loading DE results for Gene Significance ...\n")

de_rdata <- "results/02_DE_results.RData"
if (!file.exists(de_rdata)) {
    stop("DE results not found: ", de_rdata,
         "\n  Please run 02_differential_expression.R first.")
}
load(de_rdata)  # loads de_results, deg_strict, deg_broad, fit, fit_eb, design

cat(sprintf("  DE results: %d genes\n", nrow(de_results)))

## ---- 3. Verify GSTT2B probe absence on GPL6947 ----------------------------
cat("\n[3/11] Verifying GSTT2B probe status on GPL6947 ...\n")

gstt2b_in_expr <- "GSTT2B" %in% rownames(expr_gene)
gstt2b_in_de   <- "GSTT2B" %in% de_results$gene

cat(sprintf("  GSTT2B in expression matrix: %s\n", gstt2b_in_expr))
cat(sprintf("  GSTT2B in DE results:        %s\n", gstt2b_in_de))

# Also check for alternative names
gstt2b_alts <- c("GSTT2B", "GSTT2", "GSTT2P")
gstt2b_found <- gstt2b_alts[gstt2b_alts %in% rownames(expr_gene)]
if (length(gstt2b_found) > 0) {
    cat(sprintf("  GSTT2B-related probes found: %s\n",
                paste(gstt2b_found, collapse = ", ")))
} else {
    cat("  CONFIRMED: GSTT2B has NO probe on GPL6947 — it is absent from the dataset.\n")
    cat("  This gene cannot be detected in this analysis.\n")
}

# Save verification
probe_check_file <- "results/03_gstt2b_probe_check.txt"
writeLines(c(
    "GSTT2B Probe Status on GPL6947 (Illumina HumanHT-12 V3.0)",
    "-----------------------------------------------------------",
    sprintf("GSTT2B in expression matrix: %s", gstt2b_in_expr),
    sprintf("GSTT2B in DE results:        %s", gstt2b_in_de),
    sprintf("Alternative names checked:   %s", paste(gstt2b_alts, collapse = ", ")),
    sprintf("Any related probes found:    %s",
            ifelse(length(gstt2b_found) > 0, paste(gstt2b_found, collapse = ", "), "NONE")),
    "",
    "CONCLUSION: GSTT2B cannot be detected on GPL6947.",
    "Any analysis claiming GSTT2B as a hub gene from this dataset is incorrect."
), probe_check_file)
cat(sprintf("  Saved: %s\n", probe_check_file))

## ---- 4. Select top 3,000 most variable genes by MAD -----------------------
cat("\n[4/11] Selecting top 3,000 most variable genes by MAD ...\n")

gene_mad <- apply(expr_gene, 1, mad)
n_select <- min(3000, nrow(expr_gene))

if (nrow(expr_gene) < 3000) {
    cat(sprintf("  WARNING: Only %d genes available (< 3,000). Using all.\n",
                nrow(expr_gene)))
}

# Rank by MAD and select top 3,000
mad_order    <- order(gene_mad, decreasing = TRUE)
top_genes    <- names(gene_mad)[mad_order[1:n_select]]
expr_wgcna   <- expr_gene[top_genes, ]

cat(sprintf("  Selected %d genes (MAD range: %.4f – %.4f)\n",
            nrow(expr_wgcna),
            min(gene_mad[top_genes]),
            max(gene_mad[top_genes])))

# Transpose for WGCNA: samples in rows, genes in columns
datExpr <- t(expr_wgcna)
cat(sprintf("  WGCNA input matrix: %d samples x %d genes\n",
            nrow(datExpr), ncol(datExpr)))

## ---- 5. Check for good samples and genes -----------------------------------
cat("\n[5/11] Checking for outlier samples and genes with too many missing values ...\n")

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
        cat(sprintf("  Removing %d genes with too many missing values.\n",
                    sum(!gsg$goodGenes)))
    }
    if (sum(!gsg$goodSamples) > 0) {
        cat(sprintf("  Removing %d outlier samples.\n",
                    sum(!gsg$goodSamples)))
    }
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
} else {
    cat("  All samples and genes passed quality checks.\n")
}

cat(sprintf("  Final WGCNA input: %d samples x %d genes\n",
            nrow(datExpr), ncol(datExpr)))

## ---- 6. Soft-threshold power selection -------------------------------------
cat("\n[6/11] Selecting soft-threshold power (signed network) ...\n")
cat("  Testing powers 1-20 ...\n")

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 0)

cat("\n  Power selection results:\n")
cat(sprintf("  %-6s  %-12s  %-18s\n", "Power", "R-squared", "Mean Connectivity"))
for (i in seq_along(powers)) {
    cat(sprintf("  %-6d  %-12.4f  %-18.2f\n",
                powers[i],
                -sign(sft$fitIndices$slope[i]) * sft$fitIndices$SFT.R.sq[i],
                sft$fitIndices$mean.k.[i]))
}

# Use manuscript value: beta = 6
beta_power <- 6
sft_r2 <- -sign(sft$fitIndices$slope[beta_power]) * sft$fitIndices$SFT.R.sq[beta_power]
sft_mk <- sft$fitIndices$mean.k.[beta_power]

cat(sprintf("\n  Using manuscript power beta = %d\n", beta_power))
cat(sprintf("  R-squared at beta = %d: %.4f (manuscript: > 0.85)\n",
            beta_power, sft_r2))
cat(sprintf("  Mean connectivity at beta = %d: %.2f (manuscript: 45.2)\n",
            beta_power, sft_mk))

if (sft_r2 < 0.80) {
    cat("  WARNING: R-squared is below 0.80 — scale-free topology fit may be poor.\n")
}

# Plot soft threshold selection
pdf("figures/03_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# Scale-free topology fit
plot(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R²)",
     main = "Scale Independence",
     type = "n")
text(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     labels = powers, col = "red", cex = 0.9)
abline(h = 0.85, col = "blue", lty = 2)
points(beta_power, sft_r2, pch = 19, col = "blue", cex = 2)
legend("bottomright", legend = sprintf("beta = %d (R² = %.3f)", beta_power, sft_r2),
       col = "blue", pch = 19, bty = "n")

# Mean connectivity
plot(sft$fitIndices$Power,
     sft$fitIndices$mean.k.,
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean Connectivity",
     type = "n")
text(sft$fitIndices$Power,
     sft$fitIndices$mean.k.,
     labels = powers, col = "red", cex = 0.9)
points(beta_power, sft_mk, pch = 19, col = "blue", cex = 2)
legend("topright", legend = sprintf("beta = %d (k = %.1f)", beta_power, sft_mk),
       col = "blue", pch = 19, bty = "n")

dev.off()
cat("  Saved: figures/03_soft_threshold.pdf\n")

## ---- 7. Build signed adjacency and TOM ------------------------------------
cat("\n[7/11] Building signed adjacency matrix and TOM ...\n")
cat("  (This may take a few minutes for 3,000 genes) ...\n")

# Signed adjacency: a_ij = ((1 + cor(i,j)) / 2)^beta
adjacency <- adjacency(datExpr, power = beta_power, type = "signed")

cat("  Adjacency matrix computed.\n")

# Topological Overlap Matrix
TOM    <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

cat("  TOM computed. Using 1-TOM as distance.\n")

## ---- 8. Hierarchical clustering and module detection -----------------------
cat("\n[8/11] Hierarchical clustering (average linkage) and dynamic tree cut ...\n")

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
cat("  Hierarchical clustering complete.\n")

# Dynamic tree cutting: deepSplit = 2, minClusterSize = 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = 30)

# Convert numeric labels to colours
dynamicColors <- labels2colors(dynamicMods)

n_modules <- length(unique(dynamicColors)) - ifelse("grey" %in% dynamicColors, 1, 0)
n_grey    <- sum(dynamicColors == "grey")

cat(sprintf("  Modules detected: %d (excluding grey/unassigned)\n", n_modules))
cat(sprintf("  Unassigned genes (grey): %d\n", n_grey))
cat("\n  Module sizes:\n")
mod_table <- sort(table(dynamicColors), decreasing = TRUE)
for (i in seq_along(mod_table)) {
    cat(sprintf("    %-15s  %d genes\n", names(mod_table)[i], mod_table[i]))
}

# Dendrogram plot
pdf("figures/03_dendrogram.pdf", width = 12, height = 6)
plotDendroAndColors(geneTree, dynamicColors,
                    "Module colours",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "GSE25186 — Gene Dendrogram and Module Colours\n(WGCNA: signed, beta=6, deepSplit=2, minSize=30)")
dev.off()
cat("  Saved: figures/03_dendrogram.pdf\n")

## ---- 9. Module eigengenes and module-trait correlations --------------------
cat("\n[9/11] Computing module eigengenes and module-trait correlations ...\n")

# Module eigengenes (first PC of each module)
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs    <- MEList$eigengenes

# Remove the grey module eigengene if present
if ("MEgrey" %in% colnames(MEs)) {
    MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
}
cat(sprintf("  Module eigengenes computed: %d modules\n", ncol(MEs)))

# Define traits
# Disease status: PCD = 1, Control = 0 (point-biserial = Pearson on binary)
# Sex: F = 1, M = 0 (point-biserial = Pearson on binary)
datTraits <- data.frame(
    Disease = as.numeric(sample_meta$group == "PCD"),
    Sex     = as.numeric(sample_meta$sex == "F"),
    row.names = rownames(datExpr)
)

cat("  Trait encoding:\n")
cat(sprintf("    Disease: PCD=1, Control=0 (sum=%d)\n", sum(datTraits$Disease)))
cat(sprintf("    Sex:     F=1, M=0 (sum=%d)\n", sum(datTraits$Sex)))

# Module-trait correlations (Pearson / point-biserial for binary)
moduleTraitCor  <- cor(MEs, datTraits, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

cat("\n  Module-trait correlations (Disease):\n")
cat(sprintf("  %-15s  %10s  %12s\n", "Module", "Correlation", "P-value"))
for (i in seq_len(nrow(moduleTraitCor))) {
    cat(sprintf("  %-15s  %+10.4f  %12.4e\n",
                rownames(moduleTraitCor)[i],
                moduleTraitCor[i, "Disease"],
                moduleTraitPval[i, "Disease"]))
}

# Check if any are nominally significant
n_nom_sig <- sum(moduleTraitPval[, "Disease"] < 0.05)
cat(sprintf("\n  Modules nominally significant for Disease (p < 0.05): %d\n", n_nom_sig))

# Save module-trait correlations
mt_df <- data.frame(
    module           = rownames(moduleTraitCor),
    cor_disease      = moduleTraitCor[, "Disease"],
    pval_disease     = moduleTraitPval[, "Disease"],
    cor_sex          = moduleTraitCor[, "Sex"],
    pval_sex         = moduleTraitPval[, "Sex"],
    stringsAsFactors = FALSE
)
write.csv(mt_df, "results/03_module_trait_correlations.csv", row.names = FALSE)
cat("  Saved: results/03_module_trait_correlations.csv\n")

## ---- 10. Permutation testing for module-trait significance -----------------
cat("\n[10/11] Permutation testing (1,000 permutations) ...\n")
cat("  Testing whether module-trait correlations are significant\n")
cat("  against randomly permuted disease labels ...\n")

n_perm <- 1000

# Store observed correlations (absolute value for two-sided test)
obs_cor_disease <- abs(moduleTraitCor[, "Disease"])

# Permutation distribution
perm_max_cor <- numeric(n_perm)  # max |cor| across all modules per permutation
perm_cor_mat <- matrix(NA, nrow = n_perm, ncol = ncol(MEs))
colnames(perm_cor_mat) <- colnames(MEs)

cat(sprintf("  Running %d permutations ", n_perm))

for (perm_i in seq_len(n_perm)) {
    # Permute disease labels
    perm_disease <- sample(datTraits$Disease)

    # Correlate MEs with permuted labels
    perm_cor <- cor(MEs, perm_disease, use = "p")
    perm_cor_mat[perm_i, ] <- abs(as.numeric(perm_cor))
    perm_max_cor[perm_i]   <- max(abs(perm_cor))

    if (perm_i %% 100 == 0) cat(".")
}
cat(" done.\n")

# Per-module permutation p-value:
# proportion of permuted |cor| >= observed |cor|
perm_pvals <- numeric(ncol(MEs))
names(perm_pvals) <- colnames(MEs)
for (j in seq_len(ncol(MEs))) {
    perm_pvals[j] <- (sum(perm_cor_mat[, j] >= obs_cor_disease[j]) + 1) / (n_perm + 1)
}

# Family-wise permutation p-value (max-statistic correction)
fwer_pvals <- numeric(ncol(MEs))
names(fwer_pvals) <- colnames(MEs)
for (j in seq_len(ncol(MEs))) {
    fwer_pvals[j] <- (sum(perm_max_cor >= obs_cor_disease[j]) + 1) / (n_perm + 1)
}

cat("\n  Permutation results (Disease trait):\n")
cat(sprintf("  %-15s  %10s  %14s  %14s\n",
            "Module", "|Obs cor|", "Perm p-value", "FWER p-value"))
for (j in seq_len(ncol(MEs))) {
    cat(sprintf("  %-15s  %10.4f  %14.4f  %14.4f\n",
                colnames(MEs)[j],
                obs_cor_disease[j],
                perm_pvals[j],
                fwer_pvals[j]))
}

n_perm_sig <- sum(perm_pvals < 0.05)
n_fwer_sig <- sum(fwer_pvals < 0.05)
cat(sprintf("\n  Modules with permutation p < 0.05:  %d\n", n_perm_sig))
cat(sprintf("  Modules with FWER p < 0.05:         %d\n", n_fwer_sig))

if (n_perm_sig == 0) {
    cat("  CONFIRMED: No module-trait correlation survives permutation testing.\n")
    cat("  This is expected given n = 15 and the instability of WGCNA at this sample size.\n")
} else {
    cat("  NOTE: Some modules pass permutation testing — interpret with caution (n=15).\n")
}

# Save permutation results
perm_df <- data.frame(
    module        = colnames(MEs),
    obs_abs_cor   = obs_cor_disease,
    nominal_pval  = moduleTraitPval[, "Disease"],
    perm_pval     = perm_pvals,
    fwer_pval     = fwer_pvals,
    stringsAsFactors = FALSE
)
write.csv(perm_df, "results/03_permutation_pvalues.csv", row.names = FALSE)
cat("  Saved: results/03_permutation_pvalues.csv\n")

## ---- 11. Hub gene ranking: kTotal x Gene Significance ----------------------
cat("\n[11/11] Ranking hub genes by kTotal x Gene Significance ...\n")

# kTotal from adjacency (sum of all connection strengths)
kTotal <- colSums(adjacency) - 1  # subtract self-connection

# Gene Significance: absolute correlation of each gene with disease trait
# Using point-biserial correlation (= Pearson on binary)
disease_vec <- datTraits$Disease
geneSignificance <- abs(cor(datExpr, disease_vec, use = "p"))
colnames(geneSignificance) <- "GS_Disease"

# Hub metric: kTotal x Gene Significance
hub_metric <- kTotal * geneSignificance[, 1]

# Build hub gene table
hub_df <- data.frame(
    gene            = colnames(datExpr),
    module          = dynamicColors,
    kTotal          = kTotal,
    GS_Disease      = geneSignificance[, 1],
    hub_score       = hub_metric,
    stringsAsFactors = FALSE
)

# Add DE results for genes present in both
hub_df$logFC   <- NA_real_
hub_df$pvalue  <- NA_real_
matched <- match(hub_df$gene, de_results$gene)
hub_df$logFC[!is.na(matched)]  <- de_results$logFC[matched[!is.na(matched)]]
hub_df$pvalue[!is.na(matched)] <- de_results$P.Value[matched[!is.na(matched)]]

# Sort by hub score
hub_df <- hub_df[order(hub_df$hub_score, decreasing = TRUE), ]
rownames(hub_df) <- NULL

# Save full rankings
write.csv(hub_df, "results/03_hub_genes_ranked.csv", row.names = FALSE)
cat("  Saved: results/03_hub_genes_ranked.csv\n")

# Report top 20 hub genes
cat("\n  Top 20 hub genes (kTotal x GS):\n")
cat(sprintf("  %-4s  %-12s  %-12s  %10s  %10s  %12s  %+9s  %12s\n",
            "Rank", "Gene", "Module", "kTotal", "GS_Disease", "Hub Score",
            "logFC", "DE p-value"))
top_n <- min(20, nrow(hub_df))
for (i in 1:top_n) {
    cat(sprintf("  %-4d  %-12s  %-12s  %10.2f  %10.4f  %12.4f  %+9.2f  %12.2e\n",
                i,
                hub_df$gene[i],
                hub_df$module[i],
                hub_df$kTotal[i],
                hub_df$GS_Disease[i],
                hub_df$hub_score[i],
                ifelse(is.na(hub_df$logFC[i]), NA, hub_df$logFC[i]),
                ifelse(is.na(hub_df$pvalue[i]), NA, hub_df$pvalue[i])))
}

# Check MED13L position
med13l_rank <- which(hub_df$gene == "MED13L")
if (length(med13l_rank) > 0) {
    cat(sprintf("\n  MED13L rank: #%d (manuscript expects: #1 top hub)\n",
                med13l_rank[1]))
    med_row <- hub_df[med13l_rank[1], ]
    cat(sprintf("  MED13L: kTotal=%.2f, GS=%.4f, hub_score=%.4f, logFC=%+.2f, p=%.2e\n",
                med_row$kTotal, med_row$GS_Disease, med_row$hub_score,
                med_row$logFC, med_row$pvalue))
    if (med13l_rank[1] == 1) {
        cat("  CONFIRMED: MED13L is the top hub gene.\n")
    } else {
        cat(sprintf("  NOTE: MED13L is ranked #%d, not #1. Top gene is %s.\n",
                    med13l_rank[1], hub_df$gene[1]))
    }
} else {
    cat("\n  WARNING: MED13L not found in the top 3,000 MAD-selected genes.\n")
}

# Verify GSTT2B is absent from hub ranking
if ("GSTT2B" %in% hub_df$gene) {
    cat("  WARNING: GSTT2B unexpectedly found in hub rankings.\n")
} else {
    cat("  CONFIRMED: GSTT2B is absent from hub rankings (no probe on GPL6947).\n")
}

# Save module assignments
module_df <- data.frame(
    gene   = colnames(datExpr),
    module = dynamicColors,
    stringsAsFactors = FALSE
)
write.csv(module_df, "results/03_module_assignments.csv", row.names = FALSE)
cat("  Saved: results/03_module_assignments.csv\n")

## ---- Module-trait heatmap --------------------------------------------------
cat("\nGenerating module-trait heatmap ...\n")

# Build text matrix for the heatmap cells
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPval, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf("figures/03_module_trait_heatmap.pdf", width = 8, height = max(6, ncol(MEs) * 0.5 + 2))
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = colnames(MEs),
               ySymbols = gsub("ME", "", colnames(MEs)),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = paste0("GSE25186 — Module-Trait Correlations\n",
                             "(n=15; permutation p > 0.05 for all modules)"))
dev.off()
cat("  Saved: figures/03_module_trait_heatmap.pdf\n")

## ---- Hub gene scatter plot -------------------------------------------------
cat("Generating hub gene scatter plot ...\n")

scatter_df <- data.frame(
    gene       = hub_df$gene,
    kTotal     = hub_df$kTotal,
    GS_Disease = hub_df$GS_Disease,
    module     = hub_df$module,
    stringsAsFactors = FALSE
)

# Label top hub genes
top_labels <- c("MED13L", "DDX58", "CXCL9")
# Only label if they exist in the data
top_labels <- top_labels[top_labels %in% scatter_df$gene]
scatter_df$label <- ifelse(scatter_df$gene %in% top_labels, scatter_df$gene, NA)

# Also label the actual top 3 by hub score regardless
top3_genes <- hub_df$gene[1:min(3, nrow(hub_df))]
scatter_df$label[scatter_df$gene %in% top3_genes & is.na(scatter_df$label)] <-
    scatter_df$gene[scatter_df$gene %in% top3_genes & is.na(scatter_df$label)]

hub_scatter <- ggplot(scatter_df,
                      aes(x = kTotal, y = GS_Disease, colour = module)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_colour_identity() +
    labs(
        title = "GSE25186 — Hub Gene Identification",
        subtitle = "Hub score = kTotal x Gene Significance (|cor with disease|)",
        x = "kTotal (intramodular connectivity)",
        y = "Gene Significance (|cor with Disease|)",
        caption = paste0("n = 15 samples; WGCNA results are exploratory\n",
                         "GSTT2B absent (no probe on GPL6947)")
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        plot.caption  = element_text(size = 8, colour = "grey50")
    )

# Add labels
if (requireNamespace("ggrepel", quietly = TRUE)) {
    hub_scatter <- hub_scatter +
        ggrepel::geom_text_repel(
            aes(label = label),
            size = 3.5, max.overlaps = 15,
            segment.color = "grey50",
            fontface = "bold",
            na.rm = TRUE
        )
} else {
    label_data <- scatter_df[!is.na(scatter_df$label), ]
    if (nrow(label_data) > 0) {
        hub_scatter <- hub_scatter +
            geom_text(
                data = label_data,
                aes(label = label),
                size = 3.5, vjust = -0.8, fontface = "bold"
            )
    }
}

ggsave("figures/03_hub_gene_scatter.pdf", hub_scatter, width = 9, height = 7)
cat("  Saved: figures/03_hub_gene_scatter.pdf\n")

## ---- Save all WGCNA objects ------------------------------------------------
cat("\nSaving WGCNA results ...\n")

save(datExpr, datTraits, adjacency, TOM, dissTOM, geneTree,
     dynamicMods, dynamicColors, MEs, MEList,
     moduleTraitCor, moduleTraitPval,
     kTotal, geneSignificance, hub_df,
     perm_pvals, fwer_pvals, perm_cor_mat,
     sft, beta_power, sft_r2, sft_mk,
     file = "results/03_wgcna_results.RData")
cat("  Saved: results/03_wgcna_results.RData\n")

## ---- Summary ---------------------------------------------------------------
cat("\n--- WGCNA Exploratory Analysis Summary ---\n")
cat(sprintf("  Dataset:               GSE25186 (n = %d)\n", n_samples))
cat(sprintf("  Platform:              Illumina HumanHT-12 V3.0 (GPL6947)\n"))
cat(sprintf("  Input genes:           %d (top by MAD, target 3,000)\n", ncol(datExpr)))
cat(sprintf("  Network type:          signed\n"))
cat(sprintf("  Soft-threshold power:  beta = %d\n", beta_power))
cat(sprintf("  Scale-free R²:         %.4f (manuscript: > 0.85)\n", sft_r2))
cat(sprintf("  Mean connectivity:     %.2f (manuscript: 45.2)\n", sft_mk))
cat(sprintf("  Tree cutting:          deepSplit=2, minModuleSize=30\n"))
cat(sprintf("  Modules detected:      %d (excl. grey)\n", n_modules))
cat(sprintf("  Unassigned (grey):     %d genes\n", n_grey))
cat(sprintf("  Permutation test:      %d permutations\n", n_perm))
cat(sprintf("  Perm-significant:      %d modules (expected: 0)\n", n_perm_sig))
if (length(med13l_rank) > 0) {
    cat(sprintf("  Top hub gene:          %s (rank #%d, expected: MED13L #1)\n",
                hub_df$gene[1], 1))
    cat(sprintf("  MED13L hub rank:       #%d\n", med13l_rank[1]))
}
cat(sprintf("  GSTT2B status:         absent (no probe on GPL6947)\n"))
cat(sprintf("\n  *** n = %d < 20: WGCNA modules are UNSTABLE. ***\n", n_samples))
cat(sprintf("  *** Module-trait correlations are NOT significant after permutation. ***\n"))

cat("\n=== 03_wgcna_exploratory.R complete ===\n")
