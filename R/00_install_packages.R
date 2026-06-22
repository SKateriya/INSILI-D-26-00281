##############################################################################
## Script:  00_install_packages.R
## Project: INSILI-D-26-00281 — PCD transcriptomics (GSE25186)
## Purpose: Install all R packages required by the analysis pipeline.
##          Bioconductor packages are installed via BiocManager; CRAN
##          packages are installed from a fixed mirror.  Each package is
##          checked before installation so the script is safe to re-run.
##
## Inputs:  None (internet access required)
## Outputs: Side-effect only — packages installed into the user library
##
## Usage:   Rscript 00_install_packages.R
##############################################################################

cat("=== 00_install_packages.R ===\n")
cat("Checking and installing required packages ...\n\n")

## ---- helper ----------------------------------------------------------------
install_if_missing <- function(pkg, version = NULL, bioc = FALSE) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        installed_ver <- as.character(packageVersion(pkg))
        cat(sprintf("  [OK]  %-20s  (v%s already installed)\n", pkg, installed_ver))
    } else {
        cat(sprintf("  [INSTALL]  %-20s ...\n", pkg))
        if (bioc) {
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
        } else {
            install.packages(pkg)
        }
        if (!requireNamespace(pkg, quietly = TRUE)) {
            warning(sprintf("Installation of '%s' may have failed — please check manually.", pkg))
        } else {
            cat(sprintf("  [DONE]  %-20s  (v%s)\n", pkg, as.character(packageVersion(pkg))))
        }
    }
}

## ---- set CRAN mirror -------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

## ---- BiocManager first (needed for Bioconductor packages) ------------------
cat("--- Checking BiocManager ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("  [INSTALL]  BiocManager ...\n")
    install.packages("BiocManager")
} else {
    cat(sprintf("  [OK]  BiocManager  (v%s already installed)\n",
                as.character(packageVersion("BiocManager"))))
}

## ---- Bioconductor packages -------------------------------------------------
cat("\n--- Bioconductor packages ---\n")
bioc_pkgs <- c(
    "limma",
    "GEOquery",
    "org.Hs.eg.db"
)
for (pkg in bioc_pkgs) {
    install_if_missing(pkg, bioc = TRUE)
}

## ---- CRAN packages ---------------------------------------------------------
cat("\n--- CRAN packages ---\n")

# Packages with specific version expectations from the manuscript
# (versions listed for reproducibility reference; install.packages pulls latest
#  compatible release — pin manually if exact version matching is critical)
#   glmnet      4.1-8
#   randomForest 4.7-1.1
#   pROC        1.18.5

cran_pkgs <- c(
    "WGCNA",
    "glmnet",
    "randomForest",
    "pROC",
    "enrichR",
    "ggplot2",
    "pheatmap",
    "RColorBrewer"
)
for (pkg in cran_pkgs) {
    install_if_missing(pkg, bioc = FALSE)
}

## ---- version report --------------------------------------------------------
cat("\n--- Installed version summary ---\n")
all_pkgs <- c("BiocManager", bioc_pkgs, cran_pkgs)
for (pkg in all_pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  %-20s  v%s\n", pkg, as.character(packageVersion(pkg))))
    } else {
        cat(sprintf("  %-20s  ** NOT INSTALLED **\n", pkg))
    }
}

## ---- version warnings for manuscript-pinned packages -----------------------
cat("\n--- Manuscript version checks ---\n")
version_check <- function(pkg, expected) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        got <- as.character(packageVersion(pkg))
        if (got == expected) {
            cat(sprintf("  %-20s  v%s  [matches manuscript]\n", pkg, got))
        } else {
            cat(sprintf("  %-20s  v%s  [manuscript specifies v%s — results may differ]\n",
                        pkg, got, expected))
        }
    }
}
version_check("glmnet",       "4.1-8")
version_check("randomForest", "4.7-1.1")
version_check("pROC",         "1.18.5")

cat("\n=== 00_install_packages.R complete ===\n")
