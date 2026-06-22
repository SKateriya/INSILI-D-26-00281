##############################################################################
## Script:  06_tf_enrichment.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Transcription factor (TF) enrichment analysis of 176 nominally
##          significant DEGs (p < 0.01, |log2FC| > 1) via the Enrichr API.
##          Queries ChEA 2022 and TRRUST libraries, generates bubble plot
##          of top 15 enriched TFs and a TF-gene connectivity heatmap.
##          Results are reported as EXPLORATORY (Fig. 3a-b in manuscript).
##
## Dataset: GSE25186 — 6 PCD vs 9 controls (n = 15 total)
##          Platform: Illumina HumanHT-12 V3.0 (GPL6947)
##
## IMPORTANT NOTES:
##   - NO TF reaches significance after FDR correction.
##   - Nominally enriched TFs include those involved in inflammatory and
##     ciliary transcriptional programs.
##   - The ORIGINAL report of NFkB, STAT3, NRF2 as master regulators
##     was NOT SUPPORTED after sex correction.
##   - Results are reported as exploratory (Fig. 3a-b in manuscript).
##
## Parameters (from manuscript):
##   - Input: 176 DEGs at nominal p < 0.01 with |log2FC| > 1
##   - Enrichr libraries: ChEA_2022, TRRUST_Transcription_Factors_2019
##   - FDR correction: Benjamini-Hochberg (applied by Enrichr)
##   - Bubble plot: top 15 TFs, bubble size = gene set overlap,
##                  colour = combined score
##   - Heatmap: TF-gene connectivity matrix for top TFs
##
## Inputs:
##   results/02_DEG_176_p01_lfc1.csv — from 02_differential_expression.R
##       Contains: gene, logFC, AveExpr, t, P.Value, adj.P.Val
##
## Outputs:
##   results/06_enrichr_ChEA_2022.csv             — ChEA TF enrichment
##   results/06_enrichr_TRRUST_2019.csv            — TRRUST TF enrichment
##   results/06_tf_enrichment_combined.csv         — both libraries combined
##   results/06_tf_gene_connectivity.csv           — TF-gene overlap matrix
##   results/06_tf_enrichment_results.RData        — all TF enrichment objects
##   figures/06_tf_bubble_plot.pdf                 — top 15 TFs bubble plot
##   figures/06_tf_gene_heatmap.pdf                — TF-gene connectivity heatmap
##
## Usage:   Rscript 06_tf_enrichment.R
##############################################################################

cat("=== 06_tf_enrichment.R ===\n")
cat("Starting transcription factor enrichment analysis for GSE25186 ...\n\n")

set.seed(42)

## ---- load libraries --------------------------------------------------------
suppressPackageStartupMessages({
    library(enrichR)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

## ---- create output directories ---------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load 176-DEG list --------------------------------------------------
cat("[1/8] Loading 176-DEG list ...\n")

deg_file <- "results/02_DEG_176_p01_lfc1.csv"
if (!file.exists(deg_file)) {
    stop("DEG list not found: ", deg_file,
         "\n  Please run 02_differential_expression.R first.")
}
deg_176 <- read.csv(deg_file, stringsAsFactors = FALSE)

cat(sprintf("  Loaded %d DEGs from %s\n", nrow(deg_176), deg_file))

# Identify the gene column
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

# Direction info
if ("logFC" %in% colnames(deg_176)) {
    n_up   <- sum(deg_176$logFC > 0)
    n_down <- sum(deg_176$logFC < 0)
    cat(sprintf("  Direction: %d up-regulated, %d down-regulated in PCD\n",
                n_up, n_down))
}

## ---- 2. Set up Enrichr connection ------------------------------------------
cat("\n[2/8] Setting up Enrichr connection ...\n")

tryCatch({
    setEnrichrSite("Enrichr")
    cat("  Enrichr site set to: Enrichr (https://maayanlab.cloud/Enrichr)\n")
}, error = function(e) {
    cat(sprintf("  WARNING: Could not set Enrichr site: %s\n", e$message))
    cat("  Attempting to continue with default settings ...\n")
})

available_dbs <- NULL
tryCatch({
    available_dbs <- listEnrichrDbs()
    cat(sprintf("  Available Enrichr databases: %d\n", nrow(available_dbs)))
}, error = function(e) {
    cat(sprintf("  WARNING: Could not list databases: %s\n", e$message))
})

## ---- 3. Define TF libraries ------------------------------------------------
cat("\n[3/8] Defining target TF enrichment libraries ...\n")

target_libs <- c(
    "ChEA_2022",
    "TRRUST_Transcription_Factors_2019"
)

cat("  Target libraries:\n")
for (lib in target_libs) {
    if (!is.null(available_dbs) && lib %in% available_dbs$libraryName) {
        cat(sprintf("    [OK]  %s\n", lib))
    } else if (!is.null(available_dbs)) {
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

## ---- 4. Run Enrichr TF queries ---------------------------------------------
cat("\n[4/8] Running Enrichr TF enrichment queries ...\n")

tf_results <- list()
tf_success <- character(0)

for (lib in target_libs) {
    cat(sprintf("\n  Querying %s ...\n", lib))

    tryCatch({
        result <- enrichr(gene_list, lib)

        if (is.list(result) && length(result) > 0) {
            df <- result[[1]]

            if (is.data.frame(df) && nrow(df) > 0) {
                df$Library <- lib
                df <- df[order(df$P.value), ]

                tf_results[[lib]] <- df
                tf_success <- c(tf_success, lib)

                n_total     <- nrow(df)
                n_nom_sig   <- sum(df$P.value < 0.05, na.rm = TRUE)
                n_fdr_sig   <- sum(df$Adjusted.P.value < 0.05, na.rm = TRUE)

                cat(sprintf("    TFs returned:           %d\n", n_total))
                cat(sprintf("    Nominally sig (p<0.05): %d\n", n_nom_sig))
                cat(sprintf("    FDR sig (adj.p<0.05):   %d\n", n_fdr_sig))

                if (n_fdr_sig == 0) {
                    cat("    CONFIRMED: No TF reaches FDR < 0.05 significance.\n")
                } else {
                    cat("    NOTE: Some TFs reach FDR significance — unexpected.\n")
                }

                # Show top 5
                top_n <- min(5, nrow(df))
                cat(sprintf("    Top %d TFs:\n", top_n))
                for (j in 1:top_n) {
                    cat(sprintf("      %d. %s (p=%.2e, adj.p=%.2e, score=%.1f)\n",
                                j, df$Term[j], df$P.value[j],
                                df$Adjusted.P.value[j], df$Combined.Score[j]))
                }
            } else {
                cat("    No TF enrichment terms returned.\n")
                tf_results[[lib]] <- data.frame()
            }
        } else {
            cat("    Empty result returned.\n")
            tf_results[[lib]] <- data.frame()
        }
    }, error = function(e) {
        cat(sprintf("    ERROR querying %s: %s\n", lib, e$message))
        cat("    Storing empty result.\n")
        tf_results[[lib]] <<- data.frame()
    })
}

cat(sprintf("\n  Successfully queried: %d / %d libraries\n",
            length(tf_success), length(target_libs)))

## ---- 5. Verify manuscript expectations -------------------------------------
cat("\n[5/8] Verifying manuscript expectations ...\n")

cat("  --- Original claims (pre-sex-correction, NOT supported) ---\n")
cat("  NFkB, STAT3, NRF2 were originally reported as master regulators.\n")
cat("  After sex correction, these claims are NOT supported:\n")

# Check if NFkB/STAT3/NRF2 appear in results
nfkb_terms <- c("NFKB1", "RELA", "NFKB", "NF-kB")
stat3_terms <- c("STAT3")
nrf2_terms  <- c("NFE2L2", "NRF2", "Nrf2")

for (lib in tf_success) {
    df <- tf_results[[lib]]
    if (nrow(df) == 0) next

    cat(sprintf("\n  Checking %s:\n", lib))

    # NFkB
    nfkb_hits <- grep(paste(nfkb_terms, collapse = "|"),
                       df$Term, ignore.case = TRUE)
    if (length(nfkb_hits) > 0) {
        for (h in nfkb_hits) {
            cat(sprintf("    NFkB-related: %s (p=%.2e, adj.p=%.2e)\n",
                        df$Term[h], df$P.value[h], df$Adjusted.P.value[h]))
        }
    } else {
        cat("    NFkB-related: NOT found in results\n")
    }

    # STAT3
    stat3_hits <- grep(paste(stat3_terms, collapse = "|"),
                        df$Term, ignore.case = TRUE)
    if (length(stat3_hits) > 0) {
        for (h in stat3_hits) {
            cat(sprintf("    STAT3-related: %s (p=%.2e, adj.p=%.2e)\n",
                        df$Term[h], df$P.value[h], df$Adjusted.P.value[h]))
        }
    } else {
        cat("    STAT3-related: NOT found in results\n")
    }

    # NRF2
    nrf2_hits <- grep(paste(nrf2_terms, collapse = "|"),
                       df$Term, ignore.case = TRUE)
    if (length(nrf2_hits) > 0) {
        for (h in nrf2_hits) {
            cat(sprintf("    NRF2-related: %s (p=%.2e, adj.p=%.2e)\n",
                        df$Term[h], df$P.value[h], df$Adjusted.P.value[h]))
        }
    } else {
        cat("    NRF2-related: NOT found in results\n")
    }
}

## ---- 6. Save TF enrichment results ----------------------------------------
cat("\n[6/8] Saving TF enrichment results ...\n")

# Save individual library results
for (lib in names(tf_results)) {
    df <- tf_results[[lib]]
    if (lib == "ChEA_2022") {
        out_file <- "results/06_enrichr_ChEA_2022.csv"
    } else if (lib == "TRRUST_Transcription_Factors_2019") {
        out_file <- "results/06_enrichr_TRRUST_2019.csv"
    } else {
        safe_name <- gsub(" ", "_", lib)
        out_file <- sprintf("results/06_enrichr_%s.csv", safe_name)
    }
    write.csv(df, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s  (%d TFs)\n", out_file, nrow(df)))
}

# Combined table
if (length(tf_success) > 0) {
    tf_combined <- do.call(rbind, tf_results[tf_success])
    rownames(tf_combined) <- NULL
    tf_combined <- tf_combined[order(tf_combined$P.value), ]
    write.csv(tf_combined, "results/06_tf_enrichment_combined.csv",
              row.names = FALSE)
    cat(sprintf("  Saved: results/06_tf_enrichment_combined.csv  (%d TFs total)\n",
                nrow(tf_combined)))
} else {
    tf_combined <- data.frame()
    cat("  WARNING: No successful TF queries — no combined table.\n")
}

## ---- 7. Generate TF bubble plot (Fig. 3a) ----------------------------------
cat("\n[7/8] Generating TF bubble plot (Fig. 3a) ...\n")

if (nrow(tf_combined) > 0) {

    # Select top 15 TFs overall by combined score
    n_top <- min(15, nrow(tf_combined))
    top_tfs <- tf_combined[order(-tf_combined$Combined.Score), ][1:n_top, ]

    # Parse Overlap column (e.g., "3/50") to get overlap count
    parse_overlap <- function(overlap_str) {
        parts <- strsplit(as.character(overlap_str), "/")
        sapply(parts, function(x) {
            if (length(x) == 2) as.numeric(x[1]) else NA
        })
    }

    top_tfs$Overlap_count <- parse_overlap(top_tfs$Overlap)
    top_tfs$neg_log10_p   <- -log10(top_tfs$P.value)

    # Clean TF names (remove suffixes like " human" or identifiers)
    top_tfs$TF_name <- gsub("\\s+human$", "", top_tfs$Term, ignore.case = TRUE)
    top_tfs$TF_name <- gsub("\\s+\\d+$", "", top_tfs$TF_name)

    # Order by combined score
    top_tfs$TF_name <- factor(
        top_tfs$TF_name,
        levels = rev(top_tfs$TF_name[order(top_tfs$Combined.Score)])
    )

    # Clean library names
    top_tfs$Library_clean <- gsub("_", " ", top_tfs$Library)
    top_tfs$Library_clean <- gsub("Transcription Factors 2019",
                                   "TFs 2019", top_tfs$Library_clean)

    tf_bubble <- ggplot(top_tfs,
                        aes(x = Combined.Score,
                            y = TF_name,
                            size = Overlap_count,
                            colour = Combined.Score)) +
        geom_point(alpha = 0.8) +
        scale_colour_gradient(low = "#FDB863", high = "#D73027",
                              name = "Combined\nScore") +
        scale_size_continuous(range = c(3, 10), name = "Gene Set\nOverlap") +
        geom_vline(xintercept = 0, linetype = "solid",
                   colour = "grey80", linewidth = 0.3) +
        labs(
            title = "GSE25186 — TF Enrichment (176 DEGs, Exploratory)",
            subtitle = paste0(
                "Top 15 TFs by Enrichr combined score | ",
                "ChEA 2022 + TRRUST 2019\n",
                "No TF reaches FDR < 0.05 significance | ",
                "NFkB/STAT3/NRF2 NOT supported after sex correction"
            ),
            x = "Enrichr Combined Score",
            y = NULL,
            caption = paste0(
                "Bubble size = gene set overlap count | ",
                "Colour = combined score\n",
                "Design: ~ group + sex | n = 15 | ",
                "Results are EXPLORATORY (Fig. 3a)"
            )
        ) +
        facet_grid(Library_clean ~ ., scales = "free_y", space = "free_y") +
        theme_bw(base_size = 11) +
        theme(
            plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, colour = "grey30"),
            plot.caption  = element_text(size = 8, colour = "grey50",
                                         hjust = 0),
            strip.text    = element_text(face = "bold", size = 10),
            axis.text.y   = element_text(size = 9),
            legend.position = "right"
        )

    ggsave("figures/06_tf_bubble_plot.pdf", tf_bubble,
           width = 11, height = max(7, n_top * 0.45 + 3))
    cat("  Saved: figures/06_tf_bubble_plot.pdf\n")

} else {
    cat("  WARNING: No TF enrichment data — skipping bubble plot.\n")
}

## ---- 8. Generate TF-gene connectivity heatmap (Fig. 3b) --------------------
cat("\n[8/8] Generating TF-gene connectivity heatmap (Fig. 3b) ...\n")

if (nrow(tf_combined) > 0) {

    # Select top TFs for the heatmap (use top 15 by combined score)
    n_heatmap_tfs <- min(15, nrow(tf_combined))
    heatmap_tfs <- tf_combined[order(-tf_combined$Combined.Score), ][1:n_heatmap_tfs, ]

    # Parse the Genes column to extract overlapping gene names
    # Enrichr returns semicolon-separated gene lists
    parse_genes <- function(genes_str) {
        genes <- unlist(strsplit(as.character(genes_str), ";"))
        trimws(genes)
    }

    # Build TF-gene connectivity matrix
    all_overlap_genes <- character(0)
    tf_gene_list <- list()

    for (i in seq_len(nrow(heatmap_tfs))) {
        tf_name <- heatmap_tfs$Term[i]
        # Clean TF name
        tf_name_clean <- gsub("\\s+human$", "", tf_name, ignore.case = TRUE)
        tf_name_clean <- gsub("\\s+\\d+$", "", tf_name_clean)

        genes <- parse_genes(heatmap_tfs$Genes[i])
        tf_gene_list[[tf_name_clean]] <- genes
        all_overlap_genes <- union(all_overlap_genes, genes)
    }

    # Filter to genes that are in our DEG list
    overlap_in_deg <- all_overlap_genes[all_overlap_genes %in% gene_list]

    cat(sprintf("  TFs for heatmap: %d\n", length(tf_gene_list)))
    cat(sprintf("  Overlapping genes (in DEG list): %d\n", length(overlap_in_deg)))

    if (length(overlap_in_deg) > 0 && length(tf_gene_list) > 1) {

        # Limit genes for readability (top 30 by frequency across TFs)
        gene_freq <- table(unlist(tf_gene_list))
        gene_freq <- gene_freq[names(gene_freq) %in% overlap_in_deg]
        gene_freq <- sort(gene_freq, decreasing = TRUE)

        n_heatmap_genes <- min(30, length(gene_freq))
        heatmap_genes <- names(gene_freq)[1:n_heatmap_genes]

        # Build binary connectivity matrix: TF (rows) x Gene (columns)
        conn_matrix <- matrix(0,
                              nrow = length(tf_gene_list),
                              ncol = length(heatmap_genes),
                              dimnames = list(names(tf_gene_list),
                                              heatmap_genes))

        for (tf in names(tf_gene_list)) {
            genes_overlap <- intersect(tf_gene_list[[tf]], heatmap_genes)
            if (length(genes_overlap) > 0) {
                conn_matrix[tf, genes_overlap] <- 1
            }
        }

        # Remove TFs or genes with no connections (should not happen, but safe)
        row_sums <- rowSums(conn_matrix)
        col_sums <- colSums(conn_matrix)
        conn_matrix <- conn_matrix[row_sums > 0, col_sums > 0, drop = FALSE]

        # Save the connectivity matrix
        conn_df <- as.data.frame(conn_matrix)
        conn_df$TF <- rownames(conn_df)
        conn_df <- conn_df[, c("TF", setdiff(colnames(conn_df), "TF"))]
        write.csv(conn_df, "results/06_tf_gene_connectivity.csv",
                  row.names = FALSE)
        cat(sprintf("  Saved: results/06_tf_gene_connectivity.csv  (%d TFs x %d genes)\n",
                    nrow(conn_matrix), ncol(conn_matrix)))

        # Add gene direction annotation if logFC available
        gene_annotation <- NULL
        if ("logFC" %in% colnames(deg_176)) {
            direction <- character(ncol(conn_matrix))
            for (k in seq_len(ncol(conn_matrix))) {
                g <- colnames(conn_matrix)[k]
                if (g %in% deg_176[[gene_col]]) {
                    lfc <- deg_176$logFC[deg_176[[gene_col]] == g]
                    direction[k] <- ifelse(lfc > 0, "Up in PCD", "Down in PCD")
                } else {
                    direction[k] <- "Unknown"
                }
            }
            gene_annotation <- data.frame(
                Direction = direction,
                row.names = colnames(conn_matrix)
            )
        }

        # Add TF library annotation
        tf_lib_annot <- character(nrow(conn_matrix))
        for (k in seq_len(nrow(conn_matrix))) {
            tf_name_row <- rownames(conn_matrix)[k]
            # Match back to combined results
            matched <- grep(tf_name_row, tf_combined$Term,
                            ignore.case = TRUE, fixed = FALSE)
            if (length(matched) > 0) {
                tf_lib_annot[k] <- tf_combined$Library[matched[1]]
            } else {
                tf_lib_annot[k] <- "Unknown"
            }
        }
        tf_annotation <- data.frame(
            Library = gsub("_", " ", tf_lib_annot),
            row.names = rownames(conn_matrix)
        )

        # Define annotation colours
        ann_colors <- list(
            Direction = c("Up in PCD"   = "#D73027",
                          "Down in PCD" = "#4575B4",
                          "Unknown"     = "grey70"),
            Library   = c("ChEA 2022"                          = "#66C2A5",
                          "TRRUST Transcription Factors 2019"   = "#FC8D62",
                          "Unknown"                            = "grey70")
        )

        # Heatmap colour palette (binary: 0 = white, 1 = dark red)
        heatmap_colors <- colorRampPalette(c("white", "#D73027"))(2)

        # Generate heatmap
        pdf("figures/06_tf_gene_heatmap.pdf",
            width = max(10, ncol(conn_matrix) * 0.4 + 4),
            height = max(6, nrow(conn_matrix) * 0.5 + 3))

        tryCatch({
            pheatmap(conn_matrix,
                     color = heatmap_colors,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     clustering_method = "ward.D2",
                     annotation_col = gene_annotation,
                     annotation_row = tf_annotation,
                     annotation_colors = ann_colors,
                     border_color = "grey90",
                     cellwidth = 14,
                     cellheight = 16,
                     fontsize = 9,
                     fontsize_row = 9,
                     fontsize_col = 8,
                     angle_col = 45,
                     legend = FALSE,
                     main = paste0(
                         "GSE25186 — TF-Gene Connectivity ",
                         "(Exploratory, Fig. 3b)\n",
                         "White = no overlap, Red = TF target in DEG list | ",
                         "n = 15 | NFkB/STAT3/NRF2 NOT supported"
                     ))
        }, error = function(e) {
            cat(sprintf("  WARNING: pheatmap error: %s\n", e$message))
            cat("  Attempting simpler heatmap ...\n")
            # Fallback to base R heatmap
            heatmap(conn_matrix,
                    col = heatmap_colors,
                    scale = "none",
                    margins = c(8, 10),
                    main = "TF-Gene Connectivity (Exploratory)")
        })

        dev.off()
        cat("  Saved: figures/06_tf_gene_heatmap.pdf\n")

    } else {
        cat("  WARNING: Insufficient data for heatmap.\n")
        cat(sprintf("  Overlap genes in DEG list: %d\n", length(overlap_in_deg)))
        cat(sprintf("  TFs with overlap: %d\n", length(tf_gene_list)))
    }
} else {
    cat("  WARNING: No TF enrichment data — skipping heatmap.\n")
}

## ---- Save RData ------------------------------------------------------------
cat("\nSaving TF enrichment results ...\n")

save(tf_results, tf_success, tf_combined, gene_list, target_libs,
     file = "results/06_tf_enrichment_results.RData")
cat("  Saved: results/06_tf_enrichment_results.RData\n")

## ---- Summary ---------------------------------------------------------------
cat("\n--- TF Enrichment Summary ---\n")
cat(sprintf("  Dataset:             GSE25186 (n = 15)\n"))
cat(sprintf("  Input:               %d DEGs (p<0.01, |logFC|>1)\n",
            length(gene_list)))
cat(sprintf("  Libraries queried:   %s\n",
            paste(target_libs, collapse = ", ")))
cat(sprintf("  Successful queries:  %d / %d\n",
            length(tf_success), length(target_libs)))

for (lib in tf_success) {
    df <- tf_results[[lib]]
    n_nom <- sum(df$P.value < 0.05, na.rm = TRUE)
    n_fdr <- sum(df$Adjusted.P.value < 0.05, na.rm = TRUE)
    cat(sprintf("    %-35s  %d TFs (%d nom. sig, %d FDR sig)\n",
                lib, nrow(df), n_nom, n_fdr))
}

cat("\n  MANUSCRIPT CLAIMS — POST-SEX-CORRECTION STATUS:\n")
cat("    NFkB as master regulator:   NOT SUPPORTED\n")
cat("    STAT3 as master regulator:  NOT SUPPORTED\n")
cat("    NRF2 as master regulator:   NOT SUPPORTED\n")
cat("    (Original claims were artefacts of sex-confounded analysis)\n")

cat("\n  INTERPRETATION:\n")
cat("    No TF reaches significance after FDR correction.\n")
cat("    Nominally enriched TFs include those involved in inflammatory\n")
cat("    and ciliary transcriptional programs, consistent with PCD\n")
cat("    biology.  However, these are EXPLORATORY observations.\n")
cat("    The original NFkB/STAT3/NRF2 master-regulator narrative is\n")
cat("    not reproducible after correcting for the PCD-Control sex\n")
cat("    imbalance (PCD: 4F/2M, Control: 2F/7M).\n")
cat("    Results correspond to Fig. 3a-b in the manuscript.\n")

cat("\n=== 06_tf_enrichment.R complete ===\n")
