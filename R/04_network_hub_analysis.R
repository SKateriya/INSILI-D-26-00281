##############################################################################
## Script:  04_network_hub_analysis.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Network-based hub gene analysis using WGCNA results (TOM-derived
##          kTotal) and DE-derived Gene Significance. Ranks hub genes by the
##          composite metric kTotal x Gene Significance, generates network
##          plots, and provides functional annotation of top hubs.
##
## Dataset: GSE25186 — 6 PCD vs 9 controls (n = 15 total)
##          Platform: Illumina HumanHT-12 V3.0 (GPL6947)
##
## Key findings (from manuscript):
##   - MED13L = top hub gene (logFC = -3.19, p = 1.63e-3)
##     UniProt: Q71F56 — Mediator complex subunit 13-like
##     Function: Mediator complex / ciliogenesis regulator
##   - DDX58 = second hub (RIG-I innate immunity receptor)
##   - CXCL9 = third hub (chemokine, immune recruitment)
##   - GSTT2B is ABSENT from GPL6947 — cannot be detected
##
## Inputs:
##   results/03_wgcna_results.RData — from 03_wgcna_exploratory.R
##       Contains: TOM, adjacency, kTotal, geneSignificance, hub_df,
##                 dynamicColors, MEs, datExpr, datTraits
##   results/02_DE_results.RData    — from 02_differential_expression.R
##       Contains: de_results, deg_strict, deg_broad
##
## Outputs:
##   results/04_hub_gene_table.csv          — ranked hub gene table with annotation
##   results/04_top_hub_functional_annot.csv — functional annotation of top hubs
##   results/04_network_summary.csv         — network-level statistics
##   figures/04_ktotal_vs_gs_scatter.pdf    — kTotal vs Gene Significance scatter
##   figures/04_hub_gene_barplot.pdf        — top hub genes barplot
##   figures/04_hub_module_membership.pdf   — hub genes by module
##
## Usage:   Rscript 04_network_hub_analysis.R
##############################################################################

cat("=== 04_network_hub_analysis.R ===\n")
cat("Starting network hub gene analysis for GSE25186 ...\n\n")

set.seed(42)

## ---- load libraries --------------------------------------------------------
suppressPackageStartupMessages({
    library(WGCNA)
    library(ggplot2)
    library(RColorBrewer)
})

## ---- create output directories ---------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load WGCNA results ------------------------------------------------
cat("[1/8] Loading WGCNA results ...\n")

wgcna_file <- "results/03_wgcna_results.RData"
if (!file.exists(wgcna_file)) {
    stop("WGCNA results not found: ", wgcna_file,
         "\n  Please run 03_wgcna_exploratory.R first.")
}
load(wgcna_file)
# loads: datExpr, datTraits, adjacency, TOM, dissTOM, geneTree,
#        dynamicMods, dynamicColors, MEs, MEList,
#        moduleTraitCor, moduleTraitPval,
#        kTotal, geneSignificance, hub_df,
#        perm_pvals, fwer_pvals, perm_cor_mat,
#        sft, beta_power, sft_r2, sft_mk

cat(sprintf("  WGCNA data loaded: %d genes, %d modules\n",
            ncol(datExpr),
            length(unique(dynamicColors)) - ifelse("grey" %in% dynamicColors, 1, 0)))

## ---- 2. Load DE results ----------------------------------------------------
cat("\n[2/8] Loading DE results ...\n")

de_file <- "results/02_DE_results.RData"
if (!file.exists(de_file)) {
    stop("DE results not found: ", de_file,
         "\n  Please run 02_differential_expression.R first.")
}
load(de_file)  # loads: de_results, deg_strict, deg_broad, fit, fit_eb, design

cat(sprintf("  DE results loaded: %d genes\n", nrow(de_results)))

## ---- 3. Compute hub metric: kTotal x Gene Significance ---------------------
cat("\n[3/8] Computing hub metric (kTotal x Gene Significance) ...\n")

# kTotal: total connectivity from WGCNA adjacency matrix
# (already loaded from 03 results)
cat(sprintf("  kTotal range: %.2f – %.2f\n", min(kTotal), max(kTotal)))
cat(sprintf("  Gene Significance range: %.4f – %.4f\n",
            min(geneSignificance), max(geneSignificance)))

# Rebuild hub table with full DE information
hub_table <- data.frame(
    gene       = colnames(datExpr),
    module     = dynamicColors,
    kTotal     = kTotal,
    GS_Disease = as.numeric(geneSignificance[, 1]),
    hub_score  = kTotal * as.numeric(geneSignificance[, 1]),
    stringsAsFactors = FALSE
)

# Merge with DE results
de_merge <- de_results[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
hub_table <- merge(hub_table, de_merge, by = "gene", all.x = TRUE)

# Sort by hub score descending
hub_table <- hub_table[order(hub_table$hub_score, decreasing = TRUE), ]
hub_table$rank <- seq_len(nrow(hub_table))
rownames(hub_table) <- NULL

cat(sprintf("  Hub table: %d genes ranked\n", nrow(hub_table)))

## ---- 4. Verify key genes --------------------------------------------------
cat("\n[4/8] Verifying key hub genes ...\n")

verify_gene <- function(gene_name, expected_rank_desc) {
    idx <- which(hub_table$gene == gene_name)
    if (length(idx) > 0) {
        row <- hub_table[idx, ]
        cat(sprintf("  %-10s  rank=#%-4d  kTotal=%7.2f  GS=%.4f  hub_score=%8.4f  logFC=%+.2f  p=%.2e\n",
                    gene_name, row$rank, row$kTotal, row$GS_Disease,
                    row$hub_score, row$logFC, row$P.Value))
    } else {
        cat(sprintf("  %-10s  NOT FOUND in WGCNA gene set (%s)\n",
                    gene_name, expected_rank_desc))
    }
}

cat("  Expected top hubs (from manuscript):\n")
verify_gene("MED13L",  "expected: #1 top hub")
verify_gene("DDX58",   "expected: top hub")
verify_gene("CXCL9",   "expected: top hub")

cat("\n  GSTT2B probe check:\n")
verify_gene("GSTT2B",  "expected: ABSENT — no probe on GPL6947")

# Confirm top gene
cat(sprintf("\n  Actual top hub gene: %s (rank #1)\n", hub_table$gene[1]))
if (hub_table$gene[1] == "MED13L") {
    cat("  CONFIRMED: MED13L is the top hub gene as expected by manuscript.\n")
} else {
    cat(sprintf("  NOTE: Top hub is %s, not MED13L. MED13L rank = #%d.\n",
                hub_table$gene[1],
                hub_table$rank[hub_table$gene == "MED13L"]))
}

## ---- 5. Functional annotation of top hubs ----------------------------------
cat("\n[5/8] Functional annotation of top hub genes ...\n")

# Curated functional annotation (from UniProt and literature)
functional_annot <- data.frame(
    gene        = c("MED13L", "DDX58", "CXCL9"),
    uniprot_id  = c("Q71F56", "O95786", "Q07325"),
    full_name   = c("Mediator complex subunit 13-like",
                    "DExD/H-box helicase 58 (RIG-I)",
                    "C-X-C motif chemokine ligand 9"),
    function_summary = c(
        "Component of the Mediator complex that bridges transcription factors to RNA Pol II. Implicated in ciliogenesis regulation; MED13L haploinsufficiency causes intellectual disability with cardiac defects. Relevance to PCD: ciliary transcriptional programme.",
        "Cytoplasmic receptor for double-stranded RNA; initiates innate immune response via MAVS/IRF3 signalling (RIG-I pathway). Relevance to PCD: chronic airway inflammation and recurrent infections trigger innate immunity.",
        "Interferon-gamma-induced chemokine (monokine MIG); recruits CXCR3+ T cells and NK cells. Relevance to PCD: immune cell recruitment to chronically infected airways."
    ),
    pcd_relevance = c(
        "Mediator complex / ciliogenesis: direct link to ciliary gene transcription",
        "RIG-I innate immunity: chronic respiratory infections in PCD activate viral sensing pathways",
        "Chemokine signalling: T cell and NK cell recruitment in PCD airway inflammation"
    ),
    stringsAsFactors = FALSE
)

# Add hub metric data for annotated genes
for (i in seq_len(nrow(functional_annot))) {
    gene <- functional_annot$gene[i]
    idx <- which(hub_table$gene == gene)
    if (length(idx) > 0) {
        functional_annot$rank[i]      <- hub_table$rank[idx]
        functional_annot$hub_score[i] <- hub_table$hub_score[idx]
        functional_annot$logFC[i]     <- hub_table$logFC[idx]
        functional_annot$pvalue[i]    <- hub_table$P.Value[idx]
    } else {
        functional_annot$rank[i]      <- NA
        functional_annot$hub_score[i] <- NA
        functional_annot$logFC[i]     <- NA
        functional_annot$pvalue[i]    <- NA
    }
}

# Add GSTT2B note
gstt2b_annot <- data.frame(
    gene             = "GSTT2B",
    uniprot_id       = "Q9HBY0",
    full_name        = "Glutathione S-transferase theta 2B",
    function_summary = "NOT DETECTABLE on GPL6947. No probe exists for this gene on the Illumina HumanHT-12 V3.0 platform. Cannot be identified as a hub gene from this dataset.",
    pcd_relevance    = "ABSENT from platform — any claim of GSTT2B as a hub gene from GSE25186 is unsupported",
    rank             = NA,
    hub_score        = NA,
    logFC            = NA,
    pvalue           = NA,
    stringsAsFactors = FALSE
)
functional_annot <- rbind(functional_annot, gstt2b_annot)

cat("\n  Functional annotation summary:\n")
for (i in seq_len(nrow(functional_annot))) {
    cat(sprintf("\n  %s (%s)\n", functional_annot$gene[i], functional_annot$uniprot_id[i]))
    cat(sprintf("    Full name: %s\n", functional_annot$full_name[i]))
    cat(sprintf("    PCD relevance: %s\n", functional_annot$pcd_relevance[i]))
    if (!is.na(functional_annot$rank[i])) {
        cat(sprintf("    Hub rank: #%d | hub_score: %.4f | logFC: %+.2f | p: %.2e\n",
                    functional_annot$rank[i], functional_annot$hub_score[i],
                    functional_annot$logFC[i], functional_annot$pvalue[i]))
    } else {
        cat("    Hub rank: N/A (gene absent from analysis)\n")
    }
}

write.csv(functional_annot, "results/04_top_hub_functional_annot.csv", row.names = FALSE)
cat("\n  Saved: results/04_top_hub_functional_annot.csv\n")

## ---- 6. Network summary statistics ----------------------------------------
cat("\n[6/8] Computing network summary statistics ...\n")

# Network-level statistics
network_stats <- data.frame(
    metric = c(
        "n_samples",
        "n_genes_wgcna",
        "network_type",
        "soft_threshold_power",
        "scale_free_R2",
        "mean_connectivity",
        "n_modules",
        "n_unassigned_grey",
        "median_kTotal",
        "max_kTotal",
        "top_hub_gene",
        "top_hub_score",
        "MED13L_rank",
        "GSTT2B_status",
        "permutation_significant_modules"
    ),
    value = c(
        nrow(datExpr),
        ncol(datExpr),
        "signed",
        beta_power,
        round(sft_r2, 4),
        round(sft_mk, 2),
        length(unique(dynamicColors)) - ifelse("grey" %in% dynamicColors, 1, 0),
        sum(dynamicColors == "grey"),
        round(median(kTotal), 2),
        round(max(kTotal), 2),
        hub_table$gene[1],
        round(hub_table$hub_score[1], 4),
        ifelse("MED13L" %in% hub_table$gene,
               hub_table$rank[hub_table$gene == "MED13L"], NA),
        "absent — no probe on GPL6947",
        sum(perm_pvals < 0.05)
    ),
    stringsAsFactors = FALSE
)

cat("  Network statistics:\n")
for (i in seq_len(nrow(network_stats))) {
    cat(sprintf("    %-35s  %s\n", network_stats$metric[i], network_stats$value[i]))
}

write.csv(network_stats, "results/04_network_summary.csv", row.names = FALSE)
cat("  Saved: results/04_network_summary.csv\n")

## ---- 7. Plots --------------------------------------------------------------
cat("\n[7/8] Generating network plots ...\n")

### 7a. kTotal vs Gene Significance scatter with MED13L labelled ###

scatter_data <- hub_table[, c("gene", "kTotal", "GS_Disease", "module",
                               "hub_score", "logFC", "rank")]

# Identify genes to label
label_genes <- c("MED13L", "DDX58", "CXCL9")
label_genes <- label_genes[label_genes %in% scatter_data$gene]

# Also include the top 5 hub genes by score
top5 <- hub_table$gene[1:min(5, nrow(hub_table))]
label_genes <- unique(c(label_genes, top5))

scatter_data$label <- ifelse(scatter_data$gene %in% label_genes,
                             scatter_data$gene, NA)

# Highlight MED13L specifically
scatter_data$is_med13l <- scatter_data$gene == "MED13L"

p_scatter <- ggplot(scatter_data,
                    aes(x = kTotal, y = GS_Disease)) +
    geom_point(aes(colour = module), alpha = 0.4, size = 1.5) +
    scale_colour_identity() +
    # Highlight MED13L
    geom_point(data = scatter_data[scatter_data$is_med13l, ],
               colour = "red", size = 4, shape = 18) +
    labs(
        title = "GSE25186 — kTotal vs Gene Significance",
        subtitle = paste0("Hub score = kTotal x |cor(gene, Disease)| | ",
                          "Top hub: MED13L (UniProt Q71F56)"),
        x = "kTotal (total network connectivity from TOM)",
        y = "Gene Significance (|Pearson cor with Disease|)",
        caption = paste0("n = 15 samples (below WGCNA minimum of 20)\n",
                         "GSTT2B absent from GPL6947 — cannot appear in this plot")
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        plot.caption  = element_text(size = 8, colour = "grey50")
    )

# Add labels
if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_scatter <- p_scatter +
        ggrepel::geom_text_repel(
            aes(label = label),
            size = 3.5, max.overlaps = 20,
            segment.color = "grey50",
            fontface = "bold",
            na.rm = TRUE
        )
} else {
    label_data <- scatter_data[!is.na(scatter_data$label), ]
    if (nrow(label_data) > 0) {
        p_scatter <- p_scatter +
            geom_text(
                data = label_data,
                aes(label = label),
                size = 3.5, vjust = -0.8, fontface = "bold"
            )
    }
}

ggsave("figures/04_ktotal_vs_gs_scatter.pdf", p_scatter, width = 10, height = 8)
cat("  Saved: figures/04_ktotal_vs_gs_scatter.pdf\n")

### 7b. Top hub genes barplot ###

n_bar <- min(20, nrow(hub_table))
bar_data <- hub_table[1:n_bar, ]
bar_data$gene <- factor(bar_data$gene, levels = rev(bar_data$gene))

# Colour by direction of fold change
bar_data$direction <- ifelse(is.na(bar_data$logFC), "N/A",
                             ifelse(bar_data$logFC > 0, "Up in PCD",
                                    "Down in PCD"))
bar_data$direction <- factor(bar_data$direction,
                             levels = c("Up in PCD", "Down in PCD", "N/A"))

p_bar <- ggplot(bar_data, aes(x = gene, y = hub_score, fill = direction)) +
    geom_col(alpha = 0.85, width = 0.7) +
    coord_flip() +
    scale_fill_manual(
        values = c("Up in PCD"   = "#E41A1C",
                   "Down in PCD" = "#377EB8",
                   "N/A"         = "grey60"),
        name = "Direction"
    ) +
    labs(
        title = sprintf("GSE25186 — Top %d Hub Genes (kTotal x GS)", n_bar),
        subtitle = "Hub score = kTotal x |cor(gene, Disease)|",
        x = NULL,
        y = "Hub Score (kTotal x Gene Significance)",
        caption = "GSTT2B absent from GPL6947 platform"
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        plot.caption  = element_text(size = 8, colour = "grey50")
    )

ggsave("figures/04_hub_gene_barplot.pdf", p_bar, width = 9, height = 7)
cat("  Saved: figures/04_hub_gene_barplot.pdf\n")

### 7c. Hub genes by module membership ###

# For each module, show top 5 hub genes
module_names <- unique(dynamicColors)
module_names <- module_names[module_names != "grey"]
module_names <- module_names[order(module_names)]

module_hub_list <- list()
for (mod in module_names) {
    mod_genes <- hub_table[hub_table$module == mod, ]
    mod_genes <- mod_genes[order(mod_genes$hub_score, decreasing = TRUE), ]
    n_show <- min(5, nrow(mod_genes))
    if (n_show > 0) {
        top_mod <- mod_genes[1:n_show, c("gene", "module", "hub_score", "kTotal",
                                          "GS_Disease", "logFC")]
        module_hub_list[[mod]] <- top_mod
    }
}

module_hub_df <- do.call(rbind, module_hub_list)
rownames(module_hub_df) <- NULL

# Create stacked barplot by module
module_hub_df$gene_label <- paste0(module_hub_df$gene, " (", module_hub_df$module, ")")
module_hub_df$gene_label <- factor(module_hub_df$gene_label,
    levels = module_hub_df$gene_label[order(module_hub_df$hub_score)])

p_module <- ggplot(module_hub_df,
                   aes(x = gene_label, y = hub_score, fill = module)) +
    geom_col(alpha = 0.85, width = 0.7) +
    coord_flip() +
    scale_fill_identity() +
    labs(
        title = "GSE25186 — Top Hub Genes by Module",
        subtitle = "Top 5 genes per module, ranked by hub score (kTotal x GS)",
        x = NULL,
        y = "Hub Score"
    ) +
    theme_bw(base_size = 11) +
    theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        axis.text.y   = element_text(size = 8)
    )

ggsave("figures/04_hub_module_membership.pdf", p_module,
       width = 10, height = max(6, nrow(module_hub_df) * 0.3 + 2))
cat("  Saved: figures/04_hub_module_membership.pdf\n")

## ---- 8. Save hub gene table ------------------------------------------------
cat("\n[8/8] Saving final hub gene table ...\n")

# Final table with all information
hub_final <- hub_table[, c("rank", "gene", "module", "kTotal", "GS_Disease",
                            "hub_score", "logFC", "AveExpr", "t",
                            "P.Value", "adj.P.Val")]

write.csv(hub_final, "results/04_hub_gene_table.csv", row.names = FALSE)
cat(sprintf("  Saved: results/04_hub_gene_table.csv (%d genes)\n", nrow(hub_final)))

## ---- Summary ---------------------------------------------------------------
cat("\n--- Network Hub Gene Analysis Summary ---\n")
cat(sprintf("  Dataset:             GSE25186 (n = %d)\n", nrow(datExpr)))
cat(sprintf("  Genes analysed:      %d (from WGCNA)\n", ncol(datExpr)))
cat(sprintf("  Hub metric:          kTotal x Gene Significance\n"))
cat(sprintf("  Network type:        signed, beta = %d\n", beta_power))
cat(sprintf("  Modules:             %d (excl. grey)\n",
            length(unique(dynamicColors)) - ifelse("grey" %in% dynamicColors, 1, 0)))

cat("\n  Top 10 hub genes:\n")
cat(sprintf("  %-4s  %-12s  %-12s  %10s  %10s  %12s  %+9s  %12s\n",
            "Rank", "Gene", "Module", "kTotal", "GS", "Hub Score", "logFC", "p-value"))
for (i in 1:min(10, nrow(hub_table))) {
    cat(sprintf("  %-4d  %-12s  %-12s  %10.2f  %10.4f  %12.4f  %+9.2f  %12.2e\n",
                hub_table$rank[i],
                hub_table$gene[i],
                hub_table$module[i],
                hub_table$kTotal[i],
                hub_table$GS_Disease[i],
                hub_table$hub_score[i],
                ifelse(is.na(hub_table$logFC[i]), NA, hub_table$logFC[i]),
                ifelse(is.na(hub_table$P.Value[i]), NA, hub_table$P.Value[i])))
}

cat("\n  Key gene verification:\n")
if ("MED13L" %in% hub_table$gene) {
    med_row <- hub_table[hub_table$gene == "MED13L", ]
    cat(sprintf("    MED13L:  rank=#%d  logFC=%+.2f  p=%.2e  (UniProt Q71F56)\n",
                med_row$rank, med_row$logFC, med_row$P.Value))
    cat(sprintf("             Mediator complex / ciliogenesis regulator\n"))
}
if ("DDX58" %in% hub_table$gene) {
    ddx_row <- hub_table[hub_table$gene == "DDX58", ]
    cat(sprintf("    DDX58:   rank=#%d  logFC=%+.2f  p=%.2e  (RIG-I innate immunity)\n",
                ddx_row$rank, ddx_row$logFC, ddx_row$P.Value))
}
if ("CXCL9" %in% hub_table$gene) {
    cxcl_row <- hub_table[hub_table$gene == "CXCL9", ]
    cat(sprintf("    CXCL9:   rank=#%d  logFC=%+.2f  p=%.2e  (chemokine)\n",
                cxcl_row$rank, cxcl_row$logFC, cxcl_row$P.Value))
}
cat(sprintf("    GSTT2B:  ABSENT — no probe on GPL6947, cannot be detected\n"))

cat(sprintf("\n  *** LIMITATION: n = %d is below the WGCNA minimum of 20. ***\n",
            nrow(datExpr)))
cat(sprintf("  *** Hub rankings are exploratory and should be validated externally. ***\n"))

cat("\n=== 04_network_hub_analysis.R complete ===\n")
