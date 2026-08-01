#!/usr/bin/env Rscript
# =============================================================================
# Install all R/Bioconductor packages needed for the pipeline
# Run inside the 'bioinformatics' conda environment:
#   conda activate bioinformatics
#   Rscript scripts/00_install_packages.R
# =============================================================================

packages_to_check <- list(
  CRAN = c(
    # Data manipulation
    "tidyverse", "data.table",
    # Survival
    "survival", "survminer",
    # Machine learning
    "glmnet", "caret",
    # Gene sets
    "msigdbr",
    # Visualization
    "ggplot2", "pheatmap", "viridis", "ggrepel", "gridExtra",
    # Utilities
    "logger"
  ),
  Bioc = c(
    # Differential expression
    "DESeq2", "limma", "edgeR",
    # Data retrieval
    "TCGAbiolinks", "GEOquery",
    # Co-expression network
    "WGCNA",
    # Functional enrichment
    "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
    "pathview", "DOSE", "enrichplot",
    # Visualization
    "ComplexHeatmap"
  )
)

cat("=== Installing CRAN packages ===\n")
for (pkg in packages_to_check$CRAN) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s ...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

cat("\n=== Installing/updating BiocManager ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}

cat("\n=== Installing Bioconductor packages ===\n")
for (pkg in packages_to_check$Bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s ...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

cat("\n=== Verifying all packages load ===\n")
all_ok <- TRUE
for (pkg in packages_to_check$CRAN) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  MISSING: %s\n", pkg))
    all_ok <- FALSE
  }
}
for (pkg in packages_to_check$Bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  MISSING: %s\n", pkg))
    all_ok <- FALSE
  }
}

if (all_ok) {
  cat("\n✅ All packages installed successfully!\n")
} else {
  cat("\n❌ Some packages are missing. Check messages above.\n")
  quit(status = 1)
}
