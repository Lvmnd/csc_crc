#!/usr/bin/env Rscript
# =============================================================================
# 02c_align_sample_names.R — Align expression matrix colnames with pheno rownames
#
# GEO metadata uses GSM IDs as sample identifiers while supplementary
# files use various naming schemes. This script creates a mapping.
#
# Outputs: updated GSE*_expression_matrix.rds with GSM ID colnames
# =============================================================================

source("scripts/00_set_seed.R")

GEO_DIR <- file.path("data", "raw", "geo")

# =============================================================================
# 1. GSE107422 — Strip suffix from matrix colnames
#    Matrix: GSM2866686_AMC.R1 → Pheno: GSM2866686
#    Also create a title → GSM mapping for metadata merging
# =============================================================================
message("=== GSE107422: Aligning sample names ===")

mat1 <- readRDS(file.path(GEO_DIR, "GSE107422", "GSE107422_expression_matrix.rds"))
pheno1 <- read.csv(file.path(GEO_DIR, "GSE107422", "GSE107422_metadata.csv"),
                   row.names = 1, check.names = FALSE)

# Current colnames: GSM2866686_AMC.R1, ...
# Strip to just the GSM prefix
new_names <- gsub("^(GSM\\d+)_.*$", "\\1", colnames(mat1))
message("  Before: ", colnames(mat1)[1])
message("  After : ", new_names[1])

# Check uniqueness
if (any(duplicated(new_names))) {
  warning("  Duplicate names after stripping!")
}
if (sum(new_names %in% rownames(pheno1)) == 0) {
  # Try alternative: match via title
  message("  No direct match, trying title-based mapping...")
  # Matrix names have format GSM2866686_AMC.R1
  # Title has format "primary colorectal cancer AMC-R1"
  # Extract the sample code after the last underscore from matrix names
  mat_codes <- gsub("^GSM\\d+_(.*)$", "\\1", colnames(mat1))
  # Replace . with - in matrix codes (AMC.R1 → AMC-R1)
  mat_codes <- gsub("\\.", "-", mat_codes)
  
  # Extract codes from pheno titles
  pheno_titles <- pheno1$title
  pheno_codes <- gsub("^primary colorectal cancer ", "", pheno_titles)
  
  # Build mapping
  mapping <- data.frame(
    mat_name = colnames(mat1),
    mat_code = mat_codes,
    stringsAsFactors = FALSE
  )
  mapping$pheno_idx <- match(mapping$mat_code, pheno_codes)
  mapping$gsm_id <- rownames(pheno1)[mapping$pheno_idx]
  
  n_match <- sum(!is.na(mapping$gsm_id))
  message("  Matched ", n_match, " / ", nrow(mapping), " samples via title")
  
  if (n_match > 0) {
    new_names <- mapping$gsm_id
  }
}

colnames(mat1) <- new_names
saveRDS(mat1, file.path(GEO_DIR, "GSE107422", "GSE107422_expression_matrix.rds"))
message("  ✅ Updated GSE107422 expression matrix colnames")
message("  Match with pheno: ", sum(colnames(mat1) %in% rownames(pheno1)), " / ", ncol(mat1))

# (GSE220148 removed — expression-only reference cohort excluded from study)

message("\nDone. Proceed to scripts/03_qc_validation.R")
