#!/usr/bin/env Rscript
# =============================================================================
# 02_download_geo_validation.R — Download GEO validation cohorts
#
# Purpose:
#   Download RNA-seq expression data and metadata from two GEO series as
#   independent Asian-population validation cohorts for CRC.
#
# Cohorts:
#   Primary:   GSE107422 — Korean CRC, HiSeq 2000, 110 samples, recurrence data
#   (See 02d_download_geo_microarray.R for additional cohorts)
#              (expression-only replication, no outcome data available)
#
# Outputs per dataset:
#   data/raw/geo/<GSE_ID>/
#   ├── <GSE_ID>_series_matrix.txt.gz
#   ├── <GSE_ID>_metadata.csv
#   ├── <GSE_ID>_expression_matrix.rds    # parsed expression matrix
#   └── <GSE_ID>_supplementary/           # raw supplementary files
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(SummarizedExperiment)
})

OUT_DIR <- file.path("data", "raw", "geo")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Helper function to download and process one GEO series
# =============================================================================
download_geo_series <- function(gse_id, out_dir) {
  message("\n==============================================")
  message("   Processing ", gse_id)
  message("==============================================")

  series_dir <- file.path(out_dir, gse_id)
  dir.create(series_dir, showWarnings = FALSE, recursive = TRUE)

  # -- 1. Download series matrix --------------------------------------------
  message(">>> Downloading series matrix ...")
  geo_series <- getGEO(GEO = gse_id, destdir = series_dir, AnnotGPL = TRUE)

  if (length(geo_series) == 0) {
    warning("No ExpressionSet returned for ", gse_id)
    return(invisible(NULL))
  }

  eset <- geo_series[[1]]
  pheno <- pData(eset)
  write.csv(pheno, file.path(series_dir, paste0(gse_id, "_metadata.csv")))

  # -- 2. Expression data (if available in ExpressionSet) -------------------
  expr_mat <- tryCatch(exprs(eset), error = function(e) NULL)

  if (!is.null(expr_mat)) {
    message("   Expression matrix in ExpressionSet: ", nrow(expr_mat), " x ", ncol(expr_mat))
  } else {
    message("   No expression matrix in ExpressionSet (typical for RNA-seq on GEO).")
  }

  # -- 3. Supplementary files -----------------------------------------------
  message(">>> Downloading supplementary files ...")
  supp_dir <- file.path(series_dir, "supplementary")
  dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

  tryCatch({
    supp_files <- getGEOSuppFiles(gse_id, baseDir = series_dir, fetch_files = TRUE)
    # Move to supplementary subdir
    for (f in list.files(file.path(series_dir), pattern = "^", full.names = TRUE)) {
      if (!dir.exists(f) && !grepl("\\.(csv|rds|txt)", f)) {
        file.rename(f, file.path(supp_dir, basename(f)))
      }
    }
  }, error = function(e) {
    message("   Warning: could not download supplementary files: ", e$message)
  })

  # -- 4. Save ExpressionSet as RDS -----------------------------------------
  saveRDS(eset, file.path(series_dir, paste0(gse_id, "_eset.rds")))

  # -- 5. Print summary -----------------------------------------------------
  message("")
  message("   Title      : ", eset@experimentData@title)
  message("   Platform   : ", annotation(eset))
  message("   Samples    : ", ncol(eset))

  # Print non-empty phenotype columns
  for (cn in colnames(pheno)) {
    vals <- unique(pheno[[cn]])
    n_na <- sum(is.na(pheno[[cn]]) | pheno[[cn]] == "")
    if (length(vals) > 1 && n_na < ncol(eset)) {
      message(sprintf("     %-40s : %d unique", cn, length(vals)))
    }
  }

  # Save the ExpressionSet with a clean name
  assign(gsub("-", "_", gse_id), eset, envir = .GlobalEnv)

  message(">>> ", gse_id, " complete.")
  invisible(eset)
}

# =============================================================================
# Download validation cohorts
# =============================================================================

# -- Primary validation cohort ------------------------------------------------
# GSE107422 — Korean CRC, 110 samples, recurrence outcomes
eset1 <- download_geo_series("GSE107422", OUT_DIR)

# -- Secondary validation cohort (expression-only replication) ----------------
# GSE220148 removed — expression-only reference cohort excluded from study

# -- Combined summary ---------------------------------------------------------
message("")
message("==============================================")
message("   VALIDATION COHORTS — SUMMARY               ")
message("==============================================")
message("")
message(sprintf("  %-12s  %-20s  %3d  %s",
                "GSE107422", "Korean (Asan MC)", 110, "Recurrence outcome"))

message("")
message(">>> GSE107422 downloaded. See 02d_download_geo_microarray.R for remaining cohorts. Proceed to scripts/03_qc_validation.R")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))
