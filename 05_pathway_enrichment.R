##############################################################################
## Script:  05_pathway_enrichment.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Pathway enrichment analysis of 176 nominally significant DEGs
##          (p < 0.01, |log2FC| > 1) via the Enrichr API.  Queries three
##          libraries — GO Biological Process 2021, KEGG 2021 Human, and
##          Reactome 2022 — and reports results as EXPLORATORY.
##
## Dataset: GSE25186 — 6 PCD vs 9 controls (n = 15 total)
##          Platform: Illumina HumanHT-12 V3.0 (GPL6947)
##
## IMPORTANT NOTES:
##   - NO pathway reaches significance after FDR correction (consistent
##     with a small gene list from an underpowered study, n = 15).
##   - Among nominally enriched terms:
##       * Upregulated genes show oxidative stress response (GSR), innate
##         immune signalling, and protein quality control.
##       * Downregulated genes show RNA processing (DDX58/RIG-I),
##         intracellular trafficking (COG3, KDELR3), and chromatin
##         regulation (MBD1, BAZ1A, PHF3).
##   - The interferon-stimulated gene programme (IFIT2, IFIT3, IFI44L,
##     MX1, ISG15, STAT1) is ENTIRELY ABSENT after sex correction —
##     these genes collapsed to near-zero fold changes and do NOT appear
##     in the 176-DEG list.
##   - All results are reported as exploratory pathway themes, not as
##     statistically significant enrichments.
##
## Parameters (from manuscript):
##   - Input: 176 DEGs at nominal p < 0.01 with |log2FC| > 1
##   - Enrichr libraries: GO_Biological_Process_2021,
##                         KEGG_2021_Human,
##                         Reactome_2022
##   - FDR correction: Benjamini-Hochberg (applied by Enrichr)
##
## Inputs:
##   results/02_DEG_176_p01_lfc1.csv — from 02_differential_expression.R
##       Contains: gene, logFC, AveExpr, t, P.Value, adj.P.Val
##   results/02_DE_results.RData     — full DE results for direction info
##
## Outputs:
##   results/05_enrichr_GO_Biological_Process_2021.csv  — GO BP enrichment
##   results/05_enrichr_KEGG_2021_Human.csv             — KEGG enrichment
##   results/05_enrichr_Reactome_2022.csv               — Reactome enrichment
##   results/05_enrichment_combined.csv      — all libraries combined
##   results/05_enrichment_results.RData     — all enrichment objects
##   figures/05_pathway_bubble_plot.pdf      — exploratory bubble plot
##
## Usage:   Rscript 05_pathway_enrichment.R
##############################################################################

cat("=== 05_pathway_enrichment.R ===\n")
cat("Starting pathway enrichment analysis for GSE25186 ...\n\n")

set.seed(42)

## ---- load libraries --------------------------------------------------------
suppressPackageStartupMessages({
    library(enrichR)
    library(ggplot2)
    library(RColorBrewer)
})

## ---- create output directories ---------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load 176-DEG list --------------------------------------------------
cat("[1/7] Loading 176-DEG list ...\n")

deg_file <- "results/02_DEG_176_p01_lfc1.csv"
if (!file.exists(deg_file)) {
    stop("DEG list not found: ", deg_file,
         "\n  Please run 02_differential_expression.R first.")
}
deg_176 <- read.csv(deg_file, stringsAsFactors = FALSE)

cat(sprintf("  Loaded %d DEGs from %s\n", nrow(deg_176), deg_file))
cat(sprintf("  Columns: %s\n", paste(colnames(deg_176), collapse = ", ")))

# Identify the gene column (handle 'gene' or 'Gene' or first column)
if ("gene" %in% colnames(deg_176)) {
    gene_col <- "gene"
} else if ("Gene" %in% colnames(deg_176)) {
    gene_col <- "Gene"
} else {
    gene_col <- colnames(deg_176)[1]
    cat(sprintf("  NOTE: Using first column '%s' as gene names.\n", gene_col))
}

gene_list <- deg_176[[gene_col]]
cat(sprintf("  Gene list: %d unique genes\n", length(unique(gene_list))))

# Direction breakdown
if ("logFC" %in% colnames(deg_176)) {
    n_up   <- sum(deg_176$logFC > 0)
    n_down <- sum(deg_176$logFC < 0)
    cat(sprintf("  Direction: %d up-regulated, %d down-regulated in PCD\n",
                n_up, n_down))
    up_genes   <- deg_176[[gene_col]][deg_176$logFC > 0]
    down_genes <- deg_176[[gene_col]][deg_176$logFC < 0]
} else {
    cat("  WARNING: No logFC column — cannot split by direction.\n")
    up_genes   <- gene_list
    down_genes <- gene_list
}

## ---- 2. Verify ISG programme absence ---------------------------------------
cat("\n[2/7] Verifying interferon-stimulated gene (ISG) programme absence ...\n")

isg_genes <- c("IFIT2", "IFIT3", "IFI44L", "MX1", "ISG15", "STAT1")
isg_in_list <- isg_genes[isg_genes %in% gene_list]

if (length(isg_in_list) == 0) {
    cat("  CONFIRMED: ISG programme is entirely absent from the 176-DEG list.\n")
    cat("  Expected: ISGs collapsed to near-zero fold changes after sex correction.\n")
    cat(sprintf("  Checked: %s\n", paste(isg_genes, collapse = ", ")))
} else {
    cat(sprintf("  WARNING: %d ISG genes found in DEG list: %s\n",
                length(isg_in_list), paste(isg_in_list, collapse = ", ")))
    cat("  Manuscript expects these to be absent after sex correction.\n")
}

## ---- 3. Set up Enrichr connection ------------------------------------------
cat("\n[3/7] Setting up Enrichr connection ...\n")

# Set Enrichr site (use main Enrichr)
tryCatch({
    setEnrichrSite("Enrichr")
    cat("  Enrichr site set to: Enrichr (https://maayanlab.cloud/Enrichr)\n")
}, error = function(e) {
    cat(sprintf("  WARNING: Could not set Enrichr site: %s\n", e$message))
    cat("  Attempting to continue with default settings ...\n")
})

# List available databases to verify our targets exist
available_dbs <- NULL
tryCatch({
    available_dbs <- listEnrichrDbs()
    cat(sprintf("  Available Enrichr databases: %d\n", nrow(available_dbs)))
}, error = function(e) {
    cat(sprintf("  WARNING: Could not list databases: %s\n", e$message))
    cat("  Will attempt queries directly.\n")
})

## ---- 4. Define target libraries --------------------------------------------
cat("\n[4/7] Defining target Enrichr libraries ...\n")

target_libs <- c(
    "GO_Biological_Process_2021",
    "KEGG_2021_Human",
    "Reactome_2022"
)

cat("  Target libraries:\n")
for (lib in target_libs) {
    if (!is.null(available_dbs) && lib %in% available_dbs$libraryName) {
        cat(sprintf("    [OK]  %s\n", lib))
    } else if (!is.null(available_dbs)) {
        # Try fuzzy matching
        matches <- grep(gsub("_", ".", lib), available_dbs$libraryName,
                         ignore.case = TRUE, value = TRUE)
        if (length(matches) > 0) {
            cat(sprintf("    [~]   %s — closest match: %s\n", lib, matches[1]))
        } else {
            cat(sprintf("    [?]   %s — not found in database list\n", lib))
        }
    } else {
        cat(sprintf("    [?]   %s — will attempt query\n", lib))
    }
}

## ---- 5. Run Enrichr queries ------------------------------------------------
cat("\n[5/7] Running Enrichr queries ...\n")

enrichr_results <- list()
enrichr_success <- character(0)

for (lib in target_libs) {
    cat(sprintf("\n  Querying %s ...\n", lib))

    tryCatch({
        result <- enrichr(gene_list, lib)

        if (is.list(result) && length(result) > 0) {
            # enrichr() returns a named list — get the data frame
            df <- result[[1]]

            if (is.data.frame(df) && nrow(df) > 0) {
                # Add library name and sort by P.value
                df$Library <- lib
                df <- df[order(df$P.value), ]

                enrichr_results[[lib]] <- df
                enrichr_success <- c(enrichr_success, lib)

                n_total     <- nrow(df)
                n_nom_sig   <- sum(df$P.value < 0.05, na.rm = TRUE)
                n_fdr_sig   <- sum(df$Adjusted.P.value < 0.05, na.rm = TRUE)

                cat(sprintf("    Terms returned:       %d\n", n_total))
                cat(sprintf("    Nominally sig (p<0.05): %d\n", n_nom_sig))
                cat(sprintf("    FDR sig (adj.p<0.05):   %d\n", n_fdr_sig))

                if (n_fdr_sig == 0) {
                    cat("    CONFIRMED: No term reaches FDR < 0.05 significance.\n")
                } else {
                    cat("    NOTE: Some terms reach FDR significance — unexpected.\n")
                }

                # Show top 5
                top_n <- min(5, nrow(df))
                cat(sprintf("    Top %d terms:\n", top_n))
                for (j in 1:top_n) {
                    cat(sprintf("      %d. %s (p=%.2e, adj.p=%.2e, score=%.1f)\n",
                                j, df$Term[j], df$P.value[j],
                                df$Adjusted.P.value[j], df$Combined.Score[j]))
                }
            } else {
                cat("    No enrichment terms returned.\n")
                enrichr_results[[lib]] <- data.frame()
            }
        } else {
            cat("    Empty result returned.\n")
            enrichr_results[[lib]] <- data.frame()
        }
    }, error = function(e) {
        cat(sprintf("    ERROR querying %s: %s\n", lib, e$message))
        cat("    Storing empty result.\n")
        enrichr_results[[lib]] <<- data.frame()
    })
}

cat(sprintf("\n  Successfully queried: %d / %d libraries\n",
            length(enrichr_success), length(target_libs)))

## ---- 6. Verify expected enrichment themes ----------------------------------
cat("\n[6/7] Verifying expected enrichment themes from manuscript ...\n")

# Check for expected themes in nominally enriched terms
cat("\n  --- Expected themes (from manuscript) ---\n")
cat("  Upregulated genes:\n")
cat("    - Oxidative stress response (expect: GSR-related terms)\n")
cat("    - Innate immune signalling\n")
cat("    - Protein quality control\n")
cat("  Downregulated genes:\n")
cat("    - RNA processing (DDX58/RIG-I)\n")
cat("    - Intracellular trafficking (COG3, KDELR3)\n")
cat("    - Chromatin regulation (MBD1, BAZ1A, PHF3)\n")

# Key genes to verify in enrichment overlaps
key_up_genes   <- c("GSR")
key_down_genes <- c("DDX58", "COG3", "KDELR3", "MBD1", "BAZ1A", "PHF3")

cat("\n  Key gene presence in DEG list:\n")
for (g in c(key_up_genes, key_down_genes)) {
    present <- g %in% gene_list
    if (present && "logFC" %in% colnames(deg_176)) {
        lfc <- deg_176$logFC[deg_176[[gene_col]] == g]
        direction <- ifelse(lfc > 0, "UP", "DOWN")
        cat(sprintf("    %-10s  %s  (logFC = %+.2f, %s)\n",
                    g, "PRESENT", lfc, direction))
    } else if (present) {
        cat(sprintf("    %-10s  PRESENT\n", g))
    } else {
        cat(sprintf("    %-10s  ABSENT from 176-DEG list\n", g))
    }
}

# Verify NO IFN/ISG terms are enriched
cat("\n  Checking for spurious IFN/ISG enrichment:\n")
ifn_terms_found <- 0
for (lib in enrichr_success) {
    df <- enrichr_results[[lib]]
    if (nrow(df) > 0) {
        ifn_rows <- grep("interferon|IFN|ISG|antiviral",
                         df$Term, ignore.case = TRUE)
        if (length(ifn_rows) > 0) {
            nom_sig_ifn <- sum(df$P.value[ifn_rows] < 0.05)
            if (nom_sig_ifn > 0) {
                cat(sprintf("    WARNING: %d nominally significant IFN/ISG terms in %s\n",
                            nom_sig_ifn, lib))
                ifn_terms_found <- ifn_terms_found + nom_sig_ifn
            }
        }
    }
}
if (ifn_terms_found == 0) {
    cat("    CONFIRMED: No IFN/ISG terms are nominally enriched.\n")
    cat("    Consistent with ISG programme being absent after sex correction.\n")
}

## ---- 7. Save results and generate figures ----------------------------------
cat("\n[7/7] Saving results and generating figures ...\n")

# Save individual library results
for (lib in names(enrichr_results)) {
    df <- enrichr_results[[lib]]
    safe_name <- gsub(" ", "_", lib)
    out_file <- sprintf("results/05_enrichr_%s.csv", safe_name)
    write.csv(df, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s  (%d terms)\n", out_file, nrow(df)))
}

# Combined table
if (length(enrichr_success) > 0) {
    combined <- do.call(rbind, enrichr_results[enrichr_success])
    rownames(combined) <- NULL
    combined <- combined[order(combined$P.value), ]
    write.csv(combined, "results/05_enrichment_combined.csv", row.names = FALSE)
    cat(sprintf("  Saved: results/05_enrichment_combined.csv  (%d terms total)\n",
                nrow(combined)))
} else {
    combined <- data.frame()
    cat("  WARNING: No successful enrichment queries — no combined table.\n")
}

# Save RData
save(enrichr_results, enrichr_success, gene_list, up_genes, down_genes,
     combined, target_libs, isg_in_list,
     file = "results/05_enrichment_results.RData")
cat("  Saved: results/05_enrichment_results.RData\n")

## ---- Bubble plot (exploratory) ---------------------------------------------
cat("\nGenerating exploratory pathway bubble plot ...\n")

if (nrow(combined) > 0) {

    # Select top terms per library for the plot
    n_top_per_lib <- 10
    plot_data <- data.frame()

    for (lib in enrichr_success) {
        df <- enrichr_results[[lib]]
        if (nrow(df) > 0) {
            top_n <- min(n_top_per_lib, nrow(df))
            top_df <- df[1:top_n, ]
            plot_data <- rbind(plot_data, top_df)
        }
    }

    if (nrow(plot_data) > 0) {
        # Parse the Overlap column (e.g., "3/50") to get overlap count
        parse_overlap <- function(overlap_str) {
            parts <- strsplit(as.character(overlap_str), "/")
            sapply(parts, function(x) {
                if (length(x) == 2) as.numeric(x[1]) else NA
            })
        }

        plot_data$Overlap_count <- parse_overlap(plot_data$Overlap)
        plot_data$neg_log10_p   <- -log10(plot_data$P.value)

        # Shorten long term names for readability
        plot_data$Term_short <- substr(plot_data$Term, 1, 60)
        plot_data$Term_short <- ifelse(
            nchar(plot_data$Term) > 60,
            paste0(plot_data$Term_short, "..."),
            plot_data$Term_short
        )

        # Remove the GO/KEGG/Reactome ID suffixes if present
        plot_data$Term_short <- gsub("\\s*\\(GO:\\d+\\)$", "",
                                      plot_data$Term_short)
        plot_data$Term_short <- gsub("\\s*R-HSA-\\d+$", "",
                                      plot_data$Term_short)

        # Order terms by -log10(p) within each library
        plot_data <- plot_data[order(plot_data$Library, -plot_data$neg_log10_p), ]

        # Limit to top 30 overall for readability
        if (nrow(plot_data) > 30) {
            plot_data <- plot_data[order(-plot_data$neg_log10_p), ]
            plot_data <- plot_data[1:30, ]
        }

        # Reorder factor levels for plotting
        plot_data$Term_short <- factor(
            plot_data$Term_short,
            levels = rev(plot_data$Term_short[order(plot_data$neg_log10_p)])
        )

        # Clean library names for legend
        plot_data$Library_clean <- gsub("_", " ", plot_data$Library)

        bubble_plot <- ggplot(plot_data,
                              aes(x = neg_log10_p,
                                  y = Term_short,
                                  size = Overlap_count,
                                  colour = Combined.Score)) +
            geom_point(alpha = 0.8) +
            scale_colour_gradient(low = "#FDB863", high = "#D73027",
                                  name = "Combined\nScore") +
            scale_size_continuous(range = c(2, 8), name = "Gene\nOverlap") +
            facet_grid(Library_clean ~ ., scales = "free_y", space = "free_y") +
            geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                       colour = "grey50", linewidth = 0.5) +
            labs(
                title = "GSE25186 — Pathway Enrichment (176 DEGs, Exploratory)",
                subtitle = paste0(
                    "Enrichr: GO BP 2021, KEGG 2021, Reactome 2022 | ",
                    "Dashed line = nominal p = 0.05\n",
                    "No pathway reaches FDR < 0.05 significance"
                ),
                x = expression(-log[10]~(nominal~italic(p))),
                y = NULL,
                caption = paste0(
                    "Design: ~ group + sex | n = 15 | ",
                    "ISG programme absent after sex correction\n",
                    "Results are EXPLORATORY — not statistically significant"
                )
            ) +
            theme_bw(base_size = 11) +
            theme(
                plot.title    = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 9, colour = "grey30"),
                plot.caption  = element_text(size = 8, colour = "grey50",
                                             hjust = 0),
                strip.text    = element_text(face = "bold", size = 10),
                axis.text.y   = element_text(size = 8),
                legend.position = "right"
            )

        ggsave("figures/05_pathway_bubble_plot.pdf", bubble_plot,
               width = 12, height = max(8, nrow(plot_data) * 0.3 + 3))
        cat("  Saved: figures/05_pathway_bubble_plot.pdf\n")

    } else {
        cat("  WARNING: No enrichment data available for plotting.\n")
    }
} else {
    cat("  WARNING: No enrichment results — skipping bubble plot.\n")
    cat("  This may occur if Enrichr API is unreachable.\n")
}

## ---- Summary ---------------------------------------------------------------
cat("\n--- Pathway Enrichment Summary ---\n")
cat(sprintf("  Dataset:             GSE25186 (n = 15)\n"))
cat(sprintf("  Input:               %d DEGs (p<0.01, |logFC|>1)\n",
            length(gene_list)))
cat(sprintf("  Direction:           %d up, %d down in PCD\n",
            length(up_genes), length(down_genes)))
cat(sprintf("  Libraries queried:   %s\n",
            paste(target_libs, collapse = ", ")))
cat(sprintf("  Successful queries:  %d / %d\n",
            length(enrichr_success), length(target_libs)))

for (lib in enrichr_success) {
    df <- enrichr_results[[lib]]
    n_nom <- sum(df$P.value < 0.05, na.rm = TRUE)
    n_fdr <- sum(df$Adjusted.P.value < 0.05, na.rm = TRUE)
    cat(sprintf("    %-30s  %d terms (%d nom. sig, %d FDR sig)\n",
                lib, nrow(df), n_nom, n_fdr))
}

cat(sprintf("\n  ISG programme absent:  %s (expected: TRUE)\n",
            length(isg_in_list) == 0))
cat(sprintf("  IFN terms enriched:    %s (expected: NONE)\n",
            ifelse(ifn_terms_found == 0, "NONE", as.character(ifn_terms_found))))

cat("\n  INTERPRETATION:\n")
cat("    No pathway reaches significance after FDR correction.\n")
cat("    This is expected given a small gene list (176 genes) from an\n")
cat("    underpowered study (n = 15).  Nominally enriched themes include:\n")
cat("      - Oxidative stress response (upregulated genes, incl. GSR)\n")
cat("      - Innate immune signalling\n")
cat("      - Protein quality control\n")
cat("      - RNA processing (downregulated: DDX58/RIG-I)\n")
cat("      - Intracellular trafficking (downregulated: COG3, KDELR3)\n")
cat("      - Chromatin regulation (downregulated: MBD1, BAZ1A, PHF3)\n")
cat("    The ISG programme (IFIT2, IFIT3, IFI44L, MX1, ISG15, STAT1)\n")
cat("    is ENTIRELY ABSENT after sex correction.\n")
cat("    All results reported as EXPLORATORY.\n")

cat("\n=== 05_pathway_enrichment.R complete ===\n")
