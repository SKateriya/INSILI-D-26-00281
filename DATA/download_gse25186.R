#!/usr/bin/env Rscript
# =============================================================================
# Download and preprocess GSE25186 expression data from NCBI GEO
# Dataset: Geremek et al. 2014, PLoS ONE 9:e88216
# Platform: Illumina HumanHT-12 V3.0 BeadChip (GPL6947)
# Samples: 6 PCD patients + 9 healthy controls (nasal epithelial brushings)
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
})

cat("=== Downloading GSE25186 from NCBI GEO ===\n")

# Create data directory if needed
data_dir <- file.path(dirname(sys.frame(1)$ofile), ".")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Download the GEO dataset
gse <- getGEO("GSE25186", destdir = data_dir, GSEMatrix = TRUE, AnnotGPL = TRUE)
eset <- gse[[1]]

# Extract expression matrix
expr_matrix <- exprs(eset)
cat(sprintf("Expression matrix: %d probes x %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))

# Extract feature (probe) annotation
feature_data <- fData(eset)
probe_to_gene <- data.frame(
  probe_id = rownames(feature_data),
  gene_symbol = feature_data$"Gene symbol",
  gene_name = feature_data$"Gene title",
  stringsAsFactors = FALSE
)

# Save expression matrix
write.csv(expr_matrix, file.path(data_dir, "GSE25186_expression_matrix.csv"),
          row.names = TRUE)
cat("Saved: GSE25186_expression_matrix.csv\n")

# Save probe-to-gene mapping
write.csv(probe_to_gene, file.path(data_dir, "GPL6947_probe_annotation.csv"),
          row.names = FALSE)
cat("Saved: GPL6947_probe_annotation.csv\n")

# Load and display sample metadata
meta <- read.csv(file.path(data_dir, "sample_metadata.csv"))
cat("\nSample metadata:\n")
cat(sprintf("  PCD: n=%d (%d Female, %d Male)\n",
    sum(meta$group == "PCD"),
    sum(meta$group == "PCD" & meta$sex == "Female"),
    sum(meta$group == "PCD" & meta$sex == "Male")))
cat(sprintf("  Control: n=%d (%d Female, %d Male)\n",
    sum(meta$group == "Control"),
    sum(meta$group == "Control" & meta$sex == "Female"),
    sum(meta$group == "Control" & meta$sex == "Male")))

cat("\n=== Data download complete ===\n")
cat("NOTE: Sex assignments are determined from XIST/RPS4Y1 expression in\n")
cat("01_preprocessing_qc.R. The metadata file provides the expected assignments.\n")
