##############################################################################
## Script:  07_nested_loocv.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Fully nested leave-one-out cross-validation (LOOCV) for
##          transcriptomic diagnostic signature assessment.
##
##          Implements BOTH a leaky (Pipeline A) and correct nested
##          (Pipeline B) approach, demonstrating the data-leakage gap.
##          This is the most methodologically critical script in the paper
##          and directly implements Section 2.4 of the manuscript.
##
## Dataset: GSE25186 — 6 PCD vs 9 controls (n = 15 total)
##          Platform: Illumina HumanHT-12 V3.0 (GPL6947)
##
## Methodology (manuscript Section 2.4):
##   "To assess the feasibility of a transcriptomic diagnostic signature, we
##    implemented fully nested leave-one-out cross-validation (LOOCV). In each
##    of the 15 outer folds, one sample was held out for testing. The entire
##    feature-selection pipeline was re-run on the remaining 14 training
##    samples only: (1) Random Forest importance ranking (500 bootstrap-
##    aggregated trees, Gini impurity) on training data; (2) selection of
##    top 50 candidates; (3) LASSO-regularised logistic regression with
##    5-fold inner CV to determine optimal lambda; (4) final feature set
##    locked; (5) Random Forest classification. Recursive feature elimination
##    was used as a secondary confirmation step. The held-out sample never
##    participates in feature selection, standardisation, or model fitting.
##    Statistical significance was assessed by 1,000 permutations of disease
##    labels and 1,000 bootstrap resamples for confidence intervals. Genes
##    selected in >= 8/15 outer folds constitute the stable candidate
##    signature."
##
## Two pipelines:
##   Pipeline A (LEAKY):  Feature selection on ALL 15 samples, then LOOCV
##                         with pre-selected features. Expected AUC ~ 0.991.
##   Pipeline B (NESTED): Feature selection INSIDE each outer fold on 14
##                         training samples only. Expected AUC ~ 0.750.
##   Leakage gap:         0.991 - 0.750 = 0.241 (reported as 0.24)
##
## Parameters:
##   - Random Forest: ntree = 500, importance = TRUE (Gini impurity)
##   - LASSO: glmnet, alpha = 1, family = "binomial"
##   - Inner CV: cv.glmnet, nfolds = 5, type.measure = "auc"
##   - Top candidates per fold: 50
##   - Stable signature threshold: >= 8/15 folds
##   - Permutation test: 1,000 permutations
##   - Bootstrap CI: 1,000 resamples
##   - Global seed: set.seed(42)
##   - Per-fold seed: set.seed(42 + fold_number)
##
## Inputs:
##   results/01_normalised_expression.RData — from 01_preprocessing_qc.R
##       Contains: expr_gene (gene x sample matrix), sample_meta (metadata)
##   results/02_DE_results.RData           — from 02_differential_expression.R
##       Contains: de_results, deg_strict, deg_broad (DEG lists)
##
## Outputs:
##   results/07_nested_loocv_results.RData  — all ML objects for downstream
##   results/07_fold_predictions.csv        — fold-by-fold predictions (both)
##   results/07_gene_selection_freq.csv     — gene selection frequency table
##   results/07_stable_signature.csv        — genes in >= 8/15 folds
##   results/07_permutation_results.csv     — permutation null distribution
##   results/07_bootstrap_results.csv       — bootstrap AUC distribution
##   figures/07_roc_comparison.pdf          — ROC curve: leaky vs nested (Fig 4a)
##   figures/07_gene_selection_barplot.pdf  — gene selection frequency (Fig 4b)
##   figures/07_permutation_histogram.pdf   — permutation null dist (Fig 4c)
##   figures/07_bootstrap_distribution.pdf  — bootstrap AUC dist (Fig 4d)
##
## Expected results (from manuscript / corrected analysis):
##   Pipeline A (leaky) AUC:   0.991 +/- 0.006
##   Pipeline B (nested) AUC:  0.750, 95% CI [0.43, 1.00]
##   Leakage gap:              0.24
##   Permutation p-value:      0.062
##   Stable signature (8 genes): PDLIM3, LOC650406, C7orf29, DTX3L,
##                                CHRNA2, SERPINB4, CXCL9, C1orf187
##
## R packages:
##   randomForest (v4.7-1.1): RF importance and classification
##   glmnet (v4.1-8):         LASSO logistic regression
##   pROC (v1.18.5):          ROC curves, AUC, confidence intervals
##
## Usage:   Rscript 07_nested_loocv.R
##############################################################################

cat("=== 07_nested_loocv.R ===\n")
cat("Starting nested LOOCV machine learning analysis for GSE25186 ...\n\n")

## ---- Global seed -----------------------------------------------------------
set.seed(42)

## ---- Load libraries --------------------------------------------------------
cat("[0/10] Loading libraries ...\n")
suppressPackageStartupMessages({
    library(randomForest)   # RF importance ranking and classification
    library(glmnet)         # LASSO-regularised logistic regression
    library(pROC)           # ROC curves, AUC, confidence intervals
    library(ggplot2)        # Publication-quality plots
})
cat("  randomForest: ", as.character(packageVersion("randomForest")), "\n")
cat("  glmnet:       ", as.character(packageVersion("glmnet")), "\n")
cat("  pROC:         ", as.character(packageVersion("pROC")), "\n")

## ---- Create output directories ---------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ============================================================================
## 1. LOAD DATA
## ============================================================================
cat("\n[1/10] Loading normalised expression data and DE results ...\n")

# Load normalised expression matrix from preprocessing
rdata_01 <- "results/01_normalised_expression.RData"
if (!file.exists(rdata_01)) {
    stop("Normalised data not found: ", rdata_01,
         "\n  Please run 01_preprocessing_qc.R first.")
}
load(rdata_01)  # loads: expr_gene, sample_meta, pca_result, var_explained, pval_pc1
cat(sprintf("  Expression matrix: %d genes x %d samples\n",
            nrow(expr_gene), ncol(expr_gene)))

# Load DE results for reference (not used in feature selection here —
# feature selection is performed de novo inside each fold)
rdata_02 <- "results/02_DE_results.RData"
if (!file.exists(rdata_02)) {
    stop("DE results not found: ", rdata_02,
         "\n  Please run 02_differential_expression.R first.")
}
load(rdata_02)  # loads: de_results, deg_strict, deg_broad, fit, fit_eb, design
cat(sprintf("  DE results: %d genes tested\n", nrow(de_results)))
cat(sprintf("  Strict DEGs (p<0.01, |lfc|>1): %d\n", nrow(deg_strict)))
cat(sprintf("  Broad DEGs  (p<0.05, |lfc|>1): %d\n", nrow(deg_broad)))

## ============================================================================
## 2. PREPARE DATA MATRIX
## ============================================================================
cat("\n[2/10] Preparing data matrix ...\n")

# Binary outcome: PCD = 1, Control = 0
# sample_meta$group is a factor with levels c("Control", "PCD") from script 02
labels <- as.character(sample_meta$group)
y <- ifelse(labels == "PCD", 1, 0)
names(y) <- rownames(sample_meta)
n_total <- length(y)
n_pcd   <- sum(y == 1)
n_ctrl  <- sum(y == 0)

cat(sprintf("  Total samples:  %d\n", n_total))
cat(sprintf("  PCD:            %d\n", n_pcd))
cat(sprintf("  Control:        %d\n", n_ctrl))
cat(sprintf("  Manuscript expects: n = 15 (6 PCD, 9 Control)\n"))

# Full expression matrix: genes x samples -> transpose to samples x genes
# Use ALL genes for the feature selection pipeline (RF will rank them)
X_full <- t(as.matrix(expr_gene))
cat(sprintf("  Feature matrix:  %d samples x %d genes\n",
            nrow(X_full), ncol(X_full)))

# Verify sample order matches labels
stopifnot(identical(rownames(X_full), names(y)))

## ============================================================================
## 3. PIPELINE A: LEAKY (feature selection on ALL samples, then LOOCV)
## ============================================================================
cat("\n[3/10] Pipeline A: LEAKY feature selection (data leakage) ...\n")
cat("  WARNING: This pipeline selects features on ALL 15 samples before\n")
cat("  cross-validation. The held-out sample has already influenced feature\n")
cat("  selection, inflating performance estimates.\n\n")

set.seed(42)

# ---- Step A1: RF importance ranking on ALL 15 samples (LEAKY) ---------------
# 500 bootstrap-aggregated trees with Gini impurity criterion.
# This is METHODOLOGICALLY WRONG — the test sample participates in feature
# selection, causing information leakage.
cat("  [A1] Random Forest importance on ALL 15 samples (500 trees, Gini) ...\n")
rf_all <- randomForest(
    x         = X_full,
    y         = factor(y, levels = c(0, 1), labels = c("Control", "PCD")),
    ntree     = 500,
    importance = TRUE
)
# Extract Mean Decrease Gini (Gini impurity-based importance)
rf_importance_all <- importance(rf_all, type = 2)  # type 2 = MeanDecreaseGini
rf_imp_sorted <- sort(rf_importance_all[, 1], decreasing = TRUE)
cat(sprintf("  Ranked %d genes by Gini importance\n", length(rf_imp_sorted)))
cat(sprintf("  Top 5: %s\n",
            paste(names(rf_imp_sorted)[1:5], collapse = ", ")))

# ---- Step A2: Select top 50 candidate genes (LEAKY) -------------------------
top50_leaky <- names(rf_imp_sorted)[1:50]
cat(sprintf("  [A2] Selected top 50 candidates (leaky, all-sample ranking)\n"))

# ---- Step A3: LASSO on ALL 15 samples to lock features (LEAKY) --------------
# alpha = 1 (pure LASSO / L1 penalty), family = "binomial" for logistic
# regression, 5-fold inner CV to select optimal lambda, AUC as the measure.
cat("  [A3] LASSO logistic regression on ALL 15 samples ...\n")
X_top50_leaky <- X_full[, top50_leaky, drop = FALSE]

set.seed(42)
lasso_all <- cv.glmnet(
    x            = X_top50_leaky,
    y            = y,
    family       = "binomial",
    alpha        = 1,              # LASSO (L1 penalty)
    nfolds       = 5,              # 5-fold inner CV
    type.measure = "auc"           # optimise for AUC
)

# Extract features with non-zero LASSO coefficients at lambda.min
lasso_coefs_all <- as.matrix(coef(lasso_all, s = "lambda.min"))
lasso_selected_leaky <- rownames(lasso_coefs_all)[
    lasso_coefs_all[, 1] != 0 & rownames(lasso_coefs_all) != "(Intercept)"
]
cat(sprintf("  LASSO selected %d features (lambda.min = %.4f)\n",
            length(lasso_selected_leaky), lasso_all$lambda.min))
cat(sprintf("  Selected genes: %s\n",
            paste(lasso_selected_leaky, collapse = ", ")))

# ---- Step A4: LOOCV with RF on pre-selected features (LEAKY) ----------------
# The held-out sample has ALREADY influenced the feature set via Steps A1-A3.
# This inflates the performance estimate.
cat("  [A4] LOOCV with pre-selected features (leaky) ...\n")
leaky_predictions <- data.frame(
    sample     = names(y),
    true_label = y,
    pred_prob  = NA_real_,
    stringsAsFactors = FALSE
)

# If LASSO selected zero features, fall back to top 50
if (length(lasso_selected_leaky) == 0) {
    cat("  WARNING: LASSO selected 0 features; using top 50 from RF\n")
    features_leaky <- top50_leaky
} else {
    features_leaky <- lasso_selected_leaky
}

X_leaky <- X_full[, features_leaky, drop = FALSE]

for (i in seq_len(n_total)) {
    set.seed(42 + i)

    # Hold out sample i — but features were already chosen using it (LEAKY)
    X_train <- X_leaky[-i, , drop = FALSE]
    y_train <- y[-i]
    X_test  <- X_leaky[i, , drop = FALSE]

    # Train RF classifier (500 trees) on pre-selected features
    rf_model <- randomForest(
        x     = X_train,
        y     = factor(y_train, levels = c(0, 1), labels = c("Control", "PCD")),
        ntree = 500
    )

    # Predict probability of PCD class for held-out sample
    pred <- predict(rf_model, X_test, type = "prob")
    leaky_predictions$pred_prob[i] <- pred[1, "PCD"]
}

# Compute ROC and AUC for leaky pipeline
roc_leaky <- roc(
    response  = leaky_predictions$true_label,
    predictor = leaky_predictions$pred_prob,
    levels    = c(0, 1),
    direction = "<",
    quiet     = TRUE
)
auc_leaky <- as.numeric(auc(roc_leaky))

cat(sprintf("\n  Pipeline A (LEAKY) results:\n"))
cat(sprintf("    AUC = %.3f\n", auc_leaky))
cat(sprintf("    Manuscript expects: AUC ~ 0.991\n"))
cat(sprintf("    Features used: %d (locked before CV — LEAKY)\n",
            length(features_leaky)))

## ============================================================================
## 4. PIPELINE B: FULLY NESTED LOOCV (correct methodology)
## ============================================================================
cat("\n[4/10] Pipeline B: FULLY NESTED LOOCV (correct methodology) ...\n")
cat("  Each outer fold performs the ENTIRE feature selection pipeline on\n")
cat("  training data only. The held-out sample NEVER participates in\n")
cat("  feature selection, standardisation, or model fitting.\n\n")
cat("  Per-fold pipeline:\n")
cat("    B1. Hold out 1 sample for testing\n")
cat("    B2. RF importance ranking on 14 training samples (500 trees, Gini)\n")
cat("    B3. Select top 50 candidate genes by importance\n")
cat("    B4. LASSO logistic regression (alpha=1, 5-fold inner CV, AUC)\n")
cat("    B5. Lock features with non-zero LASSO coefficients\n")
cat("    B6. Train RF classifier (500 trees) on locked features\n")
cat("    B7. Predict held-out sample probability\n")
cat("    B8. Record: prediction, true label, selected features\n\n")

# Storage for nested LOOCV results
nested_predictions <- data.frame(
    fold       = seq_len(n_total),
    sample     = names(y),
    true_label = y,
    pred_prob  = NA_real_,
    n_rf_top50 = NA_integer_,     # number of top-50 RF candidates
    n_lasso    = NA_integer_,     # number selected by LASSO
    lambda_min = NA_real_,        # optimal lambda from inner CV
    stringsAsFactors = FALSE
)

# Track gene selection across folds: list of character vectors
gene_selection_list <- vector("list", n_total)

for (fold in seq_len(n_total)) {
    set.seed(42 + fold)

    cat(sprintf("  Fold %2d/%d: held-out = %s (true = %s) ... ",
                fold, n_total, names(y)[fold],
                ifelse(y[fold] == 1, "PCD", "Control")))

    # ---- Step B1: Hold out one sample for testing ---------------------------
    X_train <- X_full[-fold, , drop = FALSE]
    y_train <- y[-fold]
    X_test  <- X_full[fold, , drop = FALSE]

    # ---- Step B2: RF importance ranking on 14 TRAINING samples --------------
    # 500 bootstrap-aggregated trees, Gini impurity criterion.
    # This is re-computed fresh for EACH fold on training data only —
    # the held-out sample plays no role.
    rf_train <- randomForest(
        x         = X_train,
        y         = factor(y_train, levels = c(0, 1), labels = c("Control", "PCD")),
        ntree     = 500,
        importance = TRUE
    )
    rf_imp_train <- importance(rf_train, type = 2)[, 1]  # MeanDecreaseGini
    rf_imp_train <- sort(rf_imp_train, decreasing = TRUE)

    # ---- Step B3: Select top 50 candidate genes by importance ---------------
    top50_genes <- names(rf_imp_train)[1:50]
    nested_predictions$n_rf_top50[fold] <- length(top50_genes)

    # ---- Step B4: LASSO logistic regression with 5-fold inner CV ------------
    # Uses only the top 50 candidates on TRAINING data to determine lambda.
    # alpha = 1 (pure LASSO / L1), family = "binomial", type.measure = "auc".
    # With 14 training samples (5 or 6 PCD depending on fold), 5-fold inner
    # CV may occasionally produce folds with only one class. We handle this
    # with tryCatch and fall back to deviance if AUC fails.
    X_train_top50 <- X_train[, top50_genes, drop = FALSE]
    X_test_top50  <- X_test[, top50_genes, drop = FALSE]

    lasso_fit <- tryCatch({
        cv.glmnet(
            x            = X_train_top50,
            y            = y_train,
            family       = "binomial",
            alpha        = 1,
            nfolds       = 5,
            type.measure = "auc"
        )
    }, error = function(e) {
        # If AUC measure fails (e.g., single-class fold in inner CV),
        # fall back to deviance which always works
        cv.glmnet(
            x            = X_train_top50,
            y            = y_train,
            family       = "binomial",
            alpha        = 1,
            nfolds       = 5,
            type.measure = "deviance"
        )
    })

    nested_predictions$lambda_min[fold] <- lasso_fit$lambda.min

    # Extract features with non-zero LASSO coefficients at lambda.min
    lasso_coefs <- as.matrix(coef(lasso_fit, s = "lambda.min"))
    selected_genes <- rownames(lasso_coefs)[
        lasso_coefs[, 1] != 0 & rownames(lasso_coefs) != "(Intercept)"
    ]
    nested_predictions$n_lasso[fold] <- length(selected_genes)
    gene_selection_list[[fold]] <- selected_genes

    # ---- Step B5: Lock final feature set ------------------------------------
    # If LASSO selected zero features (can happen with small n and high
    # regularisation), fall back to top 10 RF-ranked genes to avoid
    # an empty model.
    if (length(selected_genes) == 0) {
        cat("LASSO->0 features, fallback to top 10 RF ... ")
        selected_genes <- names(rf_imp_train)[1:10]
        gene_selection_list[[fold]] <- selected_genes
        nested_predictions$n_lasso[fold] <- length(selected_genes)
    }

    # ---- Step B6: Train Random Forest classifier on locked features ---------
    # 500 trees, using ONLY the LASSO-selected features and ONLY the 14
    # training samples. No data from the held-out sample has touched this.
    X_train_final <- X_train[, selected_genes, drop = FALSE]
    X_test_final  <- X_test[, selected_genes, drop = FALSE]

    rf_final <- randomForest(
        x     = X_train_final,
        y     = factor(y_train, levels = c(0, 1), labels = c("Control", "PCD")),
        ntree = 500
    )

    # ---- Step B7: Predict held-out sample probability -----------------------
    pred_prob <- predict(rf_final, X_test_final, type = "prob")
    nested_predictions$pred_prob[fold] <- pred_prob[1, "PCD"]

    # ---- Step B8: Record results --------------------------------------------
    cat(sprintf("pred=%.3f, n_feat=%d\n",
                nested_predictions$pred_prob[fold],
                length(selected_genes)))
}

## ---- Compute ROC and AUC for nested pipeline -------------------------------
cat("\n  Computing ROC curve for nested pipeline ...\n")
roc_nested <- roc(
    response  = nested_predictions$true_label,
    predictor = nested_predictions$pred_prob,
    levels    = c(0, 1),
    direction = "<",
    quiet     = TRUE
)
auc_nested <- as.numeric(auc(roc_nested))
ci_nested  <- ci.auc(roc_nested, conf.level = 0.95)

cat(sprintf("\n  Pipeline B (NESTED) results:\n"))
cat(sprintf("    AUC     = %.3f\n", auc_nested))
cat(sprintf("    95%% CI  = [%.2f, %.2f]  (DeLong method)\n",
            ci_nested[1], ci_nested[3]))
cat(sprintf("    Manuscript expects: AUC = 0.750, 95%% CI [0.43, 1.00]\n"))

## ---- Leakage gap -----------------------------------------------------------
leakage_gap <- auc_leaky - auc_nested
cat(sprintf("\n  LEAKAGE GAP:\n"))
cat(sprintf("    AUC_leaky  = %.3f\n", auc_leaky))
cat(sprintf("    AUC_nested = %.3f\n", auc_nested))
cat(sprintf("    Gap        = %.3f (manuscript reports 0.24)\n", leakage_gap))

## ============================================================================
## 5. GENE SELECTION FREQUENCY ANALYSIS
## ============================================================================
cat("\n[5/10] Gene selection frequency analysis ...\n")

# Build frequency table: how many of the 15 outer folds selected each gene?
all_selected <- unlist(gene_selection_list)
gene_freq_table <- sort(table(all_selected), decreasing = TRUE)
gene_freq_df <- data.frame(
    gene      = names(gene_freq_table),
    frequency = as.integer(gene_freq_table),
    frac      = as.numeric(gene_freq_table) / n_total,
    stringsAsFactors = FALSE
)

cat(sprintf("  Total unique genes selected across all folds: %d\n",
            nrow(gene_freq_df)))

# Stable signature: genes selected in >= 8/15 folds
stable_threshold <- 8
stable_genes <- gene_freq_df[gene_freq_df$frequency >= stable_threshold, ]
cat(sprintf("  Genes in >= %d/%d folds (stable signature): %d\n",
            stable_threshold, n_total, nrow(stable_genes)))
if (nrow(stable_genes) > 0) {
    cat(sprintf("  Stable genes: %s\n",
                paste(stable_genes$gene, collapse = ", ")))
}

# Expected stable signature from manuscript
expected_stable <- c("PDLIM3", "LOC650406", "C7orf29", "DTX3L",
                      "CHRNA2", "SERPINB4", "CXCL9", "C1orf187")
cat(sprintf("\n  Manuscript expects 8 stable genes:\n"))
cat(sprintf("    %s\n", paste(expected_stable, collapse = ", ")))
cat(sprintf("  Overlap with obtained: %d/%d\n",
            sum(stable_genes$gene %in% expected_stable),
            length(expected_stable)))

# Print fold-by-fold selection details
cat("\n  Per-fold selection details:\n")
for (fold in seq_len(n_total)) {
    cat(sprintf("    Fold %2d: %2d genes selected",
                fold, length(gene_selection_list[[fold]])))
    # Show which stable genes were in this fold
    if (nrow(stable_genes) > 0) {
        stable_in_fold <- intersect(gene_selection_list[[fold]],
                                     stable_genes$gene)
        if (length(stable_in_fold) > 0) {
            cat(sprintf(" [stable: %s]",
                        paste(stable_in_fold, collapse = ", ")))
        }
    }
    cat("\n")
}

## ============================================================================
## 6. PERMUTATION TEST (1,000 permutations)
## ============================================================================
cat("\n[6/10] Permutation test (1,000 permutations) ...\n")
cat("  Shuffling disease labels and re-running full nested LOOCV for each\n")
cat("  permutation. This tests the null hypothesis that the observed AUC\n")
cat("  could arise by chance.\n\n")

n_perm <- 1000
perm_aucs <- numeric(n_perm)

for (perm in seq_len(n_perm)) {
    set.seed(42 + n_total + perm)  # unique, reproducible seed per permutation

    # Shuffle labels (break the true disease-expression relationship)
    y_perm <- sample(y)
    names(y_perm) <- names(y)

    # Run full nested LOOCV with shuffled labels
    perm_preds <- numeric(n_total)
    for (fold in seq_len(n_total)) {
        X_train <- X_full[-fold, , drop = FALSE]
        y_train <- y_perm[-fold]
        X_test  <- X_full[fold, , drop = FALSE]

        # B2: RF importance on training data (500 trees, Gini)
        rf_p <- randomForest(
            x         = X_train,
            y         = factor(y_train, levels = c(0, 1),
                               labels = c("Control", "PCD")),
            ntree     = 500,
            importance = TRUE
        )
        rf_imp_p <- importance(rf_p, type = 2)[, 1]
        top50_p  <- names(sort(rf_imp_p, decreasing = TRUE))[1:50]

        # B4: LASSO on training data with top 50
        X_train_p <- X_train[, top50_p, drop = FALSE]
        X_test_p  <- X_test[, top50_p, drop = FALSE]

        lasso_p <- tryCatch({
            cv.glmnet(
                x            = X_train_p,
                y            = y_train,
                family       = "binomial",
                alpha        = 1,
                nfolds       = 5,
                type.measure = "auc"
            )
        }, error = function(e) {
            cv.glmnet(
                x            = X_train_p,
                y            = y_train,
                family       = "binomial",
                alpha        = 1,
                nfolds       = 5,
                type.measure = "deviance"
            )
        })

        # B5: Lock features
        coefs_p <- as.matrix(coef(lasso_p, s = "lambda.min"))
        sel_p <- rownames(coefs_p)[
            coefs_p[, 1] != 0 & rownames(coefs_p) != "(Intercept)"
        ]
        if (length(sel_p) == 0) {
            sel_p <- names(sort(rf_imp_p, decreasing = TRUE))[1:10]
        }

        # B6-B7: RF classification and prediction
        rf_cls_p <- randomForest(
            x     = X_train[, sel_p, drop = FALSE],
            y     = factor(y_train, levels = c(0, 1),
                           labels = c("Control", "PCD")),
            ntree = 500
        )
        perm_preds[fold] <- predict(rf_cls_p,
                                     X_test[, sel_p, drop = FALSE],
                                     type = "prob")[1, "PCD"]
    }

    # Compute AUC for this permutation
    perm_roc <- tryCatch({
        roc(response  = y_perm,
            predictor = perm_preds,
            levels    = c(0, 1),
            direction = "<",
            quiet     = TRUE)
    }, error = function(e) NULL)

    perm_aucs[perm] <- if (!is.null(perm_roc)) as.numeric(auc(perm_roc)) else 0.5

    # Progress report every 100 permutations
    if (perm %% 100 == 0) {
        cat(sprintf("  Permutation %4d/%d done (mean null AUC so far = %.3f)\n",
                    perm, n_perm, mean(perm_aucs[1:perm])))
    }
}

# Permutation p-value: fraction of permuted AUCs >= observed AUC
perm_p_value <- mean(perm_aucs >= auc_nested)
cat(sprintf("\n  Permutation test results:\n"))
cat(sprintf("    Observed AUC:       %.3f\n", auc_nested))
cat(sprintf("    Mean null AUC:      %.3f\n", mean(perm_aucs)))
cat(sprintf("    SD null AUC:        %.3f\n", sd(perm_aucs)))
cat(sprintf("    Permutation p:      %.3f\n", perm_p_value))
cat(sprintf("    Manuscript expects: p ~ 0.062\n"))
cat(sprintf("    >= observed: %d/%d permutations\n",
            sum(perm_aucs >= auc_nested), n_perm))

## ============================================================================
## 7. BOOTSTRAP CONFIDENCE INTERVALS (1,000 resamples)
## ============================================================================
cat("\n[7/10] Bootstrap confidence intervals (1,000 resamples) ...\n")
cat("  Resampling the 15 prediction-label pairs with replacement to\n")
cat("  construct a non-parametric 95%% CI for the nested AUC.\n\n")

n_boot <- 1000
boot_aucs <- numeric(n_boot)

for (b in seq_len(n_boot)) {
    set.seed(42 + n_total + n_perm + b)

    # Resample the 15 prediction-label pairs with replacement
    boot_idx <- sample(seq_len(n_total), replace = TRUE)
    boot_labels <- nested_predictions$true_label[boot_idx]
    boot_preds  <- nested_predictions$pred_prob[boot_idx]

    # Compute AUC only if both classes are represented in the bootstrap sample
    if (length(unique(boot_labels)) == 2) {
        boot_roc <- tryCatch({
            roc(response  = boot_labels,
                predictor = boot_preds,
                levels    = c(0, 1),
                direction = "<",
                quiet     = TRUE)
        }, error = function(e) NULL)
        boot_aucs[b] <- if (!is.null(boot_roc)) as.numeric(auc(boot_roc)) else NA
    } else {
        # Single-class bootstrap sample — AUC is undefined
        boot_aucs[b] <- NA
    }
}

# Remove NAs (single-class bootstrap samples where AUC is undefined)
boot_aucs_valid <- boot_aucs[!is.na(boot_aucs)]
boot_ci <- quantile(boot_aucs_valid, probs = c(0.025, 0.975))

cat(sprintf("  Bootstrap results:\n"))
cat(sprintf("    Valid resamples:   %d/%d\n", length(boot_aucs_valid), n_boot))
cat(sprintf("    Mean boot AUC:     %.3f\n", mean(boot_aucs_valid)))
cat(sprintf("    Median boot AUC:   %.3f\n", median(boot_aucs_valid)))
cat(sprintf("    95%% CI:            [%.2f, %.2f]\n", boot_ci[1], boot_ci[2]))
cat(sprintf("    Manuscript expects: [0.43, 1.00]\n"))

## ============================================================================
## 8. FIGURES
## ============================================================================
cat("\n[8/10] Generating publication-quality figures ...\n")

# ---- Fig 4a: ROC curve comparison (leaky vs nested) -------------------------
cat("  [Fig 4a] ROC curve comparison ...\n")

pdf("figures/07_roc_comparison.pdf", width = 7, height = 7)
par(mar = c(5, 5, 4, 2))

# Plot nested ROC first (correct methodology — red, prominent)
plot(roc_nested,
     col  = "#E41A1C",       # red for nested (correct)
     lwd  = 2.5,
     main = "Machine Learning Classification: Leaky vs. Nested LOOCV",
     cex.main = 1.1,
     legacy.axes = TRUE,     # 1-Specificity on x-axis
     print.auc = FALSE,
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)",
     cex.lab = 1.1)

# Overlay leaky ROC (data-leakage — blue, for comparison)
plot(roc_leaky,
     col  = "#377EB8",       # blue for leaky
     lwd  = 2.5,
     add  = TRUE,
     print.auc = FALSE)

# Diagonal reference line (AUC = 0.5, random classifier)
abline(a = 0, b = 1, lty = 2, col = "grey50")

# Legend with AUC values and CI
legend("bottomright",
       legend = c(
           sprintf("Pipeline B (Nested): AUC = %.3f [%.2f, %.2f]",
                   auc_nested, ci_nested[1], ci_nested[3]),
           sprintf("Pipeline A (Leaky):  AUC = %.3f", auc_leaky),
           sprintf("Leakage gap = %.3f", leakage_gap)
       ),
       col = c("#E41A1C", "#377EB8", NA),
       lwd = c(2.5, 2.5, NA),
       lty = c(1, 1, NA),
       cex = 0.85,
       bty = "n")

# Annotation highlighting the leakage gap
text(0.5, 0.3,
     sprintf("Leakage gap = %.2f", leakage_gap),
     col = "grey30", cex = 1.1, font = 2)

dev.off()
cat("  Saved: figures/07_roc_comparison.pdf\n")

# ---- Fig 4b: Gene selection frequency barplot --------------------------------
cat("  [Fig 4b] Gene selection frequency barplot ...\n")

# Show top 20 most frequently selected genes (or all if fewer)
n_show <- min(20, nrow(gene_freq_df))
plot_freq_df <- gene_freq_df[1:n_show, ]
plot_freq_df$gene <- factor(plot_freq_df$gene,
                             levels = rev(plot_freq_df$gene))
plot_freq_df$is_stable <- plot_freq_df$frequency >= stable_threshold

freq_plot <- ggplot(plot_freq_df,
                     aes(x = gene, y = frequency, fill = is_stable)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = stable_threshold, linetype = "dashed",
               colour = "grey40", linewidth = 0.6) +
    annotate("text",
             x = n_show * 0.85, y = stable_threshold + 0.3,
             label = sprintf(">= %d folds = stable", stable_threshold),
             size = 3.2, colour = "grey30", hjust = 0) +
    scale_fill_manual(
        values = c("TRUE" = "#E41A1C", "FALSE" = "#999999"),
        labels = c("TRUE" = "Stable signature", "FALSE" = "Below threshold"),
        name = NULL
    ) +
    coord_flip() +
    labs(
        title = "Gene Selection Frequency Across 15 Outer LOOCV Folds",
        x = NULL,
        y = sprintf("Number of folds selected (out of %d)", n_total),
        caption = sprintf("Stable signature: %d genes in >= %d/%d folds",
                           nrow(stable_genes), stable_threshold, n_total)
    ) +
    scale_y_continuous(breaks = seq(0, n_total, by = 2),
                        limits = c(0, n_total + 0.5)) +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        legend.position = c(0.8, 0.2),
        legend.background = element_rect(fill = "white", colour = "grey80")
    )

ggsave("figures/07_gene_selection_barplot.pdf", freq_plot,
       width = 8, height = 6)
cat("  Saved: figures/07_gene_selection_barplot.pdf\n")

# ---- Fig 4c: Permutation null distribution histogram -------------------------
cat("  [Fig 4c] Permutation null distribution histogram ...\n")

perm_df <- data.frame(auc = perm_aucs)

perm_plot <- ggplot(perm_df, aes(x = auc)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = "#CCCCCC", colour = "grey50",
                   linewidth = 0.3) +
    geom_density(colour = "grey30", linewidth = 0.7) +
    geom_vline(xintercept = auc_nested, colour = "#E41A1C",
               linewidth = 1.2, linetype = "solid") +
    annotate("text",
             x = auc_nested + 0.02, y = Inf,
             label = sprintf("Observed AUC = %.3f\np = %.3f",
                              auc_nested, perm_p_value),
             hjust = 0, vjust = 1.5,
             size = 3.5, colour = "#E41A1C", fontface = "bold") +
    labs(
        title = "Permutation Test: 1,000 Label Shuffles",
        x = "AUC (permuted labels)",
        y = "Density",
        caption = sprintf("Permutation p = %.3f (%d/%d >= observed)",
                           perm_p_value,
                           sum(perm_aucs >= auc_nested), n_perm)
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 12)
    )

ggsave("figures/07_permutation_histogram.pdf", perm_plot,
       width = 7, height = 5)
cat("  Saved: figures/07_permutation_histogram.pdf\n")

# ---- Fig 4d: Bootstrap AUC distribution --------------------------------------
cat("  [Fig 4d] Bootstrap AUC distribution ...\n")

boot_df <- data.frame(auc = boot_aucs_valid)

boot_plot <- ggplot(boot_df, aes(x = auc)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = "#B2D8E8", colour = "#377EB8",
                   linewidth = 0.3) +
    geom_density(colour = "#377EB8", linewidth = 0.7) +
    geom_vline(xintercept = auc_nested, colour = "#E41A1C",
               linewidth = 1.2, linetype = "solid") +
    geom_vline(xintercept = boot_ci[1], colour = "#E41A1C",
               linewidth = 0.8, linetype = "dashed") +
    geom_vline(xintercept = boot_ci[2], colour = "#E41A1C",
               linewidth = 0.8, linetype = "dashed") +
    annotate("text",
             x = auc_nested, y = Inf,
             label = sprintf("AUC = %.3f\n95%% CI [%.2f, %.2f]",
                              auc_nested, boot_ci[1], boot_ci[2]),
             hjust = -0.1, vjust = 1.5,
             size = 3.5, colour = "#E41A1C", fontface = "bold") +
    labs(
        title = "Bootstrap AUC Distribution (1,000 Resamples)",
        x = "AUC",
        y = "Density",
        caption = sprintf("%d valid resamples (excluded %d single-class samples)",
                           length(boot_aucs_valid),
                           n_boot - length(boot_aucs_valid))
    ) +
    scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.1)) +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 12)
    )

ggsave("figures/07_bootstrap_distribution.pdf", boot_plot,
       width = 7, height = 5)
cat("  Saved: figures/07_bootstrap_distribution.pdf\n")

## ============================================================================
## 9. SAVE RESULTS
## ============================================================================
cat("\n[9/10] Saving result tables and RData ...\n")

# ---- Fold-by-fold predictions ------------------------------------------------
# Combine leaky and nested predictions into one comparison table
fold_results <- data.frame(
    fold       = seq_len(n_total),
    sample     = names(y),
    true_label = y,
    true_group = ifelse(y == 1, "PCD", "Control"),
    pred_prob_leaky  = leaky_predictions$pred_prob,
    pred_prob_nested = nested_predictions$pred_prob,
    n_features_lasso = nested_predictions$n_lasso,
    lambda_min       = nested_predictions$lambda_min,
    stringsAsFactors = FALSE
)
write.csv(fold_results, "results/07_fold_predictions.csv", row.names = FALSE)
cat("  Saved: results/07_fold_predictions.csv\n")

# ---- Gene selection frequency ------------------------------------------------
write.csv(gene_freq_df, "results/07_gene_selection_freq.csv", row.names = FALSE)
cat(sprintf("  Saved: results/07_gene_selection_freq.csv  (%d genes)\n",
            nrow(gene_freq_df)))

# ---- Stable signature -------------------------------------------------------
write.csv(stable_genes, "results/07_stable_signature.csv", row.names = FALSE)
cat(sprintf("  Saved: results/07_stable_signature.csv  (%d genes)\n",
            nrow(stable_genes)))

# ---- Permutation results ----------------------------------------------------
perm_results_df <- data.frame(
    permutation = seq_len(n_perm),
    auc         = perm_aucs
)
write.csv(perm_results_df, "results/07_permutation_results.csv",
          row.names = FALSE)
cat(sprintf("  Saved: results/07_permutation_results.csv  (%d permutations)\n",
            n_perm))

# ---- Bootstrap results ------------------------------------------------------
boot_results_df <- data.frame(
    resample = seq_len(n_boot),
    auc      = boot_aucs
)
write.csv(boot_results_df, "results/07_bootstrap_results.csv",
          row.names = FALSE)
cat(sprintf("  Saved: results/07_bootstrap_results.csv  (%d resamples)\n",
            n_boot))

# ---- RData with all objects --------------------------------------------------
save(
    # Pipeline A (leaky)
    rf_all, rf_importance_all, top50_leaky, lasso_all,
    lasso_selected_leaky, features_leaky,
    leaky_predictions, roc_leaky, auc_leaky,
    # Pipeline B (nested)
    nested_predictions, gene_selection_list,
    roc_nested, auc_nested, ci_nested,
    # Gene selection analysis
    gene_freq_df, stable_genes, stable_threshold,
    # Permutation test
    perm_aucs, perm_p_value, n_perm,
    # Bootstrap
    boot_aucs, boot_aucs_valid, boot_ci, n_boot,
    # Summary statistics
    leakage_gap,
    file = "results/07_nested_loocv_results.RData"
)
cat("  Saved: results/07_nested_loocv_results.RData\n")

## ============================================================================
## 10. VERIFICATION AND SUMMARY
## ============================================================================
cat("\n[10/10] Verification against manuscript expectations ...\n")
cat("=========================================================\n\n")

# ---- Verification checks ----------------------------------------------------
# Each check compares obtained results to manuscript-stated values.
# Tolerance reflects the expected variability due to stochastic elements
# (RF bootstrap, LASSO path, permutation draws).

checks_passed <- 0
checks_total  <- 0

verify <- function(description, obtained, expected, tolerance = NULL) {
    checks_total <<- checks_total + 1
    if (is.null(tolerance)) {
        # Character / exact comparison
        match <- identical(obtained, expected)
        status <- ifelse(match, "PASS", "WARN")
        cat(sprintf("  [%s] %s\n", status, description))
        cat(sprintf("         Obtained: %s\n", toString(obtained)))
        cat(sprintf("         Expected: %s\n", toString(expected)))
    } else {
        # Numeric comparison with tolerance
        match <- abs(obtained - expected) <= tolerance
        status <- ifelse(match, "PASS", "WARN")
        cat(sprintf("  [%s] %s\n", status, description))
        cat(sprintf("         Obtained: %.4f\n", obtained))
        cat(sprintf("         Expected: %.4f (+/- %.4f)\n", expected, tolerance))
    }
    if (status == "PASS") checks_passed <<- checks_passed + 1
    cat("\n")
}

cat("--- PIPELINE A (LEAKY) ---\n")
verify("Leaky AUC ~ 0.991", auc_leaky, 0.991, tolerance = 0.05)

cat("--- PIPELINE B (NESTED) ---\n")
verify("Nested AUC ~ 0.750", auc_nested, 0.750, tolerance = 0.15)
verify("Nested CI lower ~ 0.43", as.numeric(ci_nested[1]), 0.43, tolerance = 0.15)
verify("Nested CI upper ~ 1.00", as.numeric(ci_nested[3]), 1.00, tolerance = 0.10)

cat("--- LEAKAGE GAP ---\n")
verify("Leakage gap ~ 0.24", leakage_gap, 0.24, tolerance = 0.10)

cat("--- PERMUTATION TEST ---\n")
verify("Permutation p ~ 0.062", perm_p_value, 0.062, tolerance = 0.05)

cat("--- BOOTSTRAP ---\n")
verify("Bootstrap CI lower ~ 0.43", as.numeric(boot_ci[1]), 0.43, tolerance = 0.15)
verify("Bootstrap CI upper ~ 1.00", as.numeric(boot_ci[2]), 1.00, tolerance = 0.10)

cat("--- STABLE SIGNATURE ---\n")
verify("Number of stable genes ~ 8", nrow(stable_genes), 8, tolerance = 3)
n_overlap <- sum(stable_genes$gene %in% expected_stable)
verify("Overlap with expected stable genes",
       n_overlap, length(expected_stable), tolerance = 4)

cat(sprintf("=========================================================\n"))
cat(sprintf("  Verification: %d/%d checks within tolerance\n",
            checks_passed, checks_total))
cat(sprintf("=========================================================\n"))

# ---- Final summary ----------------------------------------------------------
cat("\n--- Nested LOOCV Summary ---\n")
cat(sprintf("  Dataset:            GSE25186 (%d samples: %d PCD, %d Control)\n",
            n_total, n_pcd, n_ctrl))
cat(sprintf("  Pipeline A (leaky):\n"))
cat(sprintf("    AUC:              %.3f\n", auc_leaky))
cat(sprintf("    Features:         %d (selected on all %d samples)\n",
            length(features_leaky), n_total))
cat(sprintf("  Pipeline B (nested):\n"))
cat(sprintf("    AUC:              %.3f\n", auc_nested))
cat(sprintf("    95%% CI (DeLong):  [%.2f, %.2f]\n",
            as.numeric(ci_nested[1]), as.numeric(ci_nested[3])))
cat(sprintf("    95%% CI (boot):    [%.2f, %.2f]\n",
            as.numeric(boot_ci[1]), as.numeric(boot_ci[2])))
cat(sprintf("  Leakage gap:        %.3f\n", leakage_gap))
cat(sprintf("  Permutation p:      %.3f (%d/%d >= observed)\n",
            perm_p_value, sum(perm_aucs >= auc_nested), n_perm))
cat(sprintf("  Stable signature:   %d genes (in >= %d/%d folds)\n",
            nrow(stable_genes), stable_threshold, n_total))
if (nrow(stable_genes) > 0) {
    cat(sprintf("    Genes: %s\n", paste(stable_genes$gene, collapse = ", ")))
}
cat(sprintf("\n  Key methodological point:\n"))
cat(sprintf("    The %.2f AUC gap demonstrates that data leakage from\n",
            leakage_gap))
cat(sprintf("    pre-CV feature selection inflates classifier performance.\n"))
cat(sprintf("    The nested AUC of %.3f with p = %.3f represents an honest\n",
            auc_nested, perm_p_value))
cat(sprintf("    assessment of diagnostic feasibility in this cohort.\n"))

cat("\n=== 07_nested_loocv.R complete ===\n")

## ============================================================================
## SESSION INFO (for reproducibility)
## ============================================================================
cat("\n--- Session Info ---\n")
print(sessionInfo())
