#!/usr/bin/env Rscript
# =============================================================================
# 02b_parse_geo_supplementary.R — Parse GEO supplementary files into matrices
#
# For GSE107422: extract individual sample .txt.gz files from RAW.tar,
#   combine into expression matrix
# (GSE220148 parsing removed — expression-only reference cohort excluded from study)
# 
# Outputs:
#   data/raw/geo/GSE107422/GSE107422_expression_matrix.rds
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(readxl)
  library(SummarizedExperiment)
})

GEO_DIR <- file.path("data", "raw", "geo")

# =============================================================================
# 1. GSE107422 — Parse individual sample files from RAW.tar
# =============================================================================
message("==============================================")
message("  Processing GSE107422 expression data")
message("==============================================")

gse1_dir <- file.path(GEO_DIR, "GSE107422")
supp_dir <- file.path(gse1_dir, "supplementary")
tar_file <- file.path(supp_dir, "GSE107422_RAW.tar")
extract_dir <- file.path(supp_dir, "extracted")
dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(tar_file) && !file.exists(file.path(gse1_dir, "GSE107422_expression_matrix.rds"))) {
  message(">>> Extracting GSE107422_RAW.tar ...")
  untar(tar_file, exdir = extract_dir)

  gz_files <- list.files(extract_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  message("   Found ", length(gz_files), " sample files")

  first <- read.table(gzfile(gz_files[1]), header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE)
  genes <- first[, 1]
  n_genes <- length(genes)

  sample_names <- gsub("_log2_TPM\\.txt\\.gz$", "", basename(gz_files))
  expr_mat <- matrix(NA, nrow = n_genes, ncol = length(gz_files))
  colnames(expr_mat) <- sample_names
  rownames(expr_mat) <- genes

  for (i in seq_along(gz_files)) {
    dat <- read.table(gzfile(gz_files[i]), header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE)
    expr_mat[, i] <- as.numeric(dat[, 2])
  }

  message("   Expression matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
  message("   Value range: ", round(min(expr_mat, na.rm = TRUE), 2), " — ",
          round(max(expr_mat, na.rm = TRUE), 2))

  saveRDS(expr_mat, file.path(gse1_dir, "GSE107422_expression_matrix.rds"))
  message("   → Saved GSE107422_expression_matrix.rds")
} else {
  message("   Already processed or tar file missing — skipping")
}

# =============================================================================
# 2. GSE220148 — REMOVED
# =============================================================================
message("\n==============================================")
message("  GSE220148 (REMOVED — expression-only reference)")
message("==============================================")
message("  This cohort was removed from the study.")
message("  See 02d_download_geo_microarray.R for additional validation cohorts.")

message("\n>>> Supplementary data parsing complete.")
