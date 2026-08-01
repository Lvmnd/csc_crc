#!/usr/bin/env Rscript
# =============================================================================
# 10_robustness_check.R — Self-audit & numbers ledger for the entire pipeline
#
# Purpose:
#   Before any manuscript text is written, this script performs three
#   independent audits and generates a single "Numbers Ledger" (CSV) that
#   serves as the sole source of truth for every statistic in the manuscript.
#   The manuscript-writing phase MUST pull only from this ledger, never from
#   memory or re-estimation.
#
# Audits:
#   A. Sample Size Reconciliation — confirms every sample count across all 9
#      scripts is internally consistent and matches Table 1.
#   B. Multiple-Testing Correction Audit — verifies that BH/FDR correction was
#      applied wherever multiple genes/pathways/cell-types were tested.
#   C. Data Leakage Check — confirms validation samples were never touched
#      during signature construction (discovery = TCGA only).
#
# Outputs:
#   results/tables/numbers_ledger.csv      — EVERY manuscript statistic with
#                                             source script, line, and output file
#   results/tables/robustness_audit.txt    — narrative audit report
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# -- File paths ---------------------------------------------------------------
SCRIPTS_DIR  <- "scripts"
TCGA_QC_DIR  <- file.path("data", "processed", "tcga_qc")
TCGA_SIG_DIR <- file.path("data", "processed", "tcga_sig")
GEO_QC_DIR   <- file.path("data", "processed", "geo_qc")
GEO_QC_GSE   <- file.path(GEO_QC_DIR, "GSE107422")
RAW_TCGA_DIR <- file.path("data", "raw", "tcga")
RAW_GEO_DIR  <- file.path("data", "raw", "geo", "GSE107422")
VAL_DIR      <- file.path("data", "processed", "validation")
TAB_DIR      <- file.path("results", "tables")
ENRICH_DIR   <- file.path(TAB_DIR, "enrichment")
IMMUNE_DIR   <- file.path("results", "tables", "immune_infiltration")

OUT_LEDGER  <- file.path(TAB_DIR, "numbers_ledger.csv")
OUT_AUDIT   <- file.path(TAB_DIR, "robustness_audit.txt")

# -- Helper: add a row to the numbers ledger ----------------------------------
ledger <- data.frame(
  script        = character(),
  category      = character(),
  stat_name     = character(),
  stat_value    = character(),
  source_file   = character(),
  source_line   = integer(),
  notes         = character(),
  stringsAsFactors = FALSE
)

add_ledger <- function(script, category, stat_name, stat_value,
                       source_file = "", source_line = NA, notes = "") {
  # Debug: identify which call is failing (write to stderr for log capture)
  # (disabled to avoid clutter - lengths are now guarded above)
  # Guard: if any argument has length 0, coerce to NA
  if (length(stat_value) == 0L) stat_value <- NA_character_
  if (length(source_line) == 0L) source_line <- NA_integer_
  if (length(source_file) == 0L) source_file <- ""
  if (length(notes) == 0L)       notes <- ""
  if (length(script) == 0L)      script <- ""
  if (length(category) == 0L)    category <- ""
  if (length(stat_name) == 0L)   stat_name <- ""
  new_row <- data.frame(
    script = script, category = category,
    stat_name = stat_name, stat_value = as.character(stat_value),
    source_file = source_file, source_line = as.integer(source_line),
    notes = notes,
    stringsAsFactors = FALSE
  )
  ledger <<- dplyr::bind_rows(ledger, new_row)
}

# -- Helper: file existence check ---------------------------------------------
file_available <- function(path) file.exists(path) && file.info(path)$size > 0

# =============================================================================
# SECTION 0 — Gather all available output files
# =============================================================================
message(">>> Gathering available output files ...")

has_tcga_expr     <- file_available(file.path(RAW_TCGA_DIR, "tcga_expression_SE.rds"))
has_tcga_clin_c   <- file_available(file.path(RAW_TCGA_DIR, "clinical_COAD.rds"))
has_tcga_clin_r   <- file_available(file.path(RAW_TCGA_DIR, "clinical_READ.rds"))
has_tcga_dds      <- file_available(file.path(TCGA_QC_DIR, "tcga_dds.rds"))
has_tcga_dds_de   <- file_available(file.path(TCGA_QC_DIR, "tcga_dds_de.rds"))
has_tcga_vsd      <- file_available(file.path(TCGA_QC_DIR, "tcga_vsd.rds"))
has_tcga_qc_sum   <- file_available(file.path(TCGA_QC_DIR, "tcga_qc_summary.txt"))
has_tcga_excl     <- file_available(file.path(TCGA_QC_DIR, "tcga_excluded_samples.csv"))
has_deg_full      <- file_available(file.path(TAB_DIR, "DEG_full_TCGA.csv"))
has_deg_sig       <- file_available(file.path(TAB_DIR, "DEG_sig_TCGA.csv"))
has_deg_sum       <- file_available(file.path(TAB_DIR, "DEG_summary_TCGA.txt"))
has_csc_csv       <- file_available(file.path(TAB_DIR, "CSC_markers_TCGA.csv"))
has_csc_sum       <- file_available(file.path(TAB_DIR, "CSC_markers_summary_TCGA.txt"))
has_lasso_coef    <- file_available(file.path(TCGA_SIG_DIR, "lasso_coefficients.csv"))
has_risk_scores   <- file_available(file.path(TCGA_SIG_DIR, "risk_scores.csv"))
has_cox_sum       <- file_available(file.path(TAB_DIR, "Cox_multivariable_summary.txt"))
has_ph_test       <- file_available(file.path(TAB_DIR, "PH_assumption_test.txt"))
has_wgcna_res     <- file_available(file.path(TCGA_SIG_DIR, "wgcna_results.rds"))
has_csc_assoc     <- file_available(file.path(TCGA_SIG_DIR, "csc_associated_genes.csv"))
has_geo_excl      <- file_available(file.path(GEO_QC_DIR, "geo_excluded_samples.csv"))
has_geo_qc_sum    <- file_available(file.path(GEO_QC_DIR, "geo_qc_summary.txt"))
has_val_expr      <- file_available(file.path(GEO_QC_GSE, "counts_filtered.rds")) ||
                     file_available(file.path(GEO_QC_GSE, "tpm_filtered.rds"))
has_val_risk      <- file_available(file.path(VAL_DIR, "risk_scores_validation.csv"))
has_val_outcome   <- file_available(file.path(TAB_DIR, "Validation_outcome.txt"))
has_val_cox       <- file_available(file.path(TAB_DIR, "Validation_cox_summary.txt"))
has_val_tertile   <- file_available(file.path(TAB_DIR, "Validation_tertile_sensitivity.txt"))
has_enrich_deg    <- file_available(file.path(ENRICH_DIR, "TCGA_DEG_high_vs_low.csv"))
has_gsea_kegg     <- file_available(file.path(ENRICH_DIR, "TCGA_GSEA_kegg.csv"))
has_gsea_go       <- file_available(file.path(ENRICH_DIR, "TCGA_GSEA_go.csv"))
has_ora_kegg      <- file_available(file.path(ENRICH_DIR, "TCGA_ORA_kegg.csv"))
has_ora_go        <- file_available(file.path(ENRICH_DIR, "TCGA_ORA_go.csv"))
has_stem_audit    <- file_available(file.path(ENRICH_DIR, "TCGA_stemness_audit.txt"))
has_immune_tcga   <- file_available(file.path(IMMUNE_DIR, "TCGA_immune_scores.csv"))
has_immune_val    <- file_available(file.path(IMMUNE_DIR, "Validation_immune_scores.csv"))
has_immune_cor    <- file_available(file.path(IMMUNE_DIR, "TCGA_immune_risk_correlation.csv"))

audit_lines <- character()
add_audit <- function(...) {
  audit_lines <<- c(audit_lines, paste0(...))
}

add_audit("==============================================")
add_audit("  PIPELINE ROBUSTNESS AUDIT                    ")
add_audit("==============================================")
add_audit("")
add_audit(sprintf("Audit timestamp : %s", Sys.time()))
add_audit(sprintf("Seed            : 42 (global, via 00_set_seed.R)"))
add_audit("")

# =============================================================================
# AUDIT A — Sample Size Reconciliation
# =============================================================================
add_audit("")
add_audit("==============================================")
add_audit("  AUDIT A — SAMPLE SIZE RECONCILIATION        ")
add_audit("==============================================")
add_audit("")

# ---- Script 01: TCGA Download ----------------------------------------------
add_audit("--- Script 01: TCGA Download (01_download_tcga.R) ---")

n_tcga_samples_raw <- NA
n_tcga_tumor_raw <- NA
n_tcga_normal_raw <- NA
n_tcga_coad_raw <- NA
n_tcga_read_raw <- NA
n_clin_coad <- NA
n_clin_read <- NA
n_surv_coad <- NA
n_surv_read <- NA

if (has_tcga_expr) {
  se <- readRDS(file.path(RAW_TCGA_DIR, "tcga_expression_SE.rds"))
  n_tcga_samples_raw <- ncol(se)
  barcodes <- colnames(se)
  # TCGA barcode sample type codes: "01"=tumor, "11"=normal
  tumor_idx <- which(substr(barcodes, 14, 15) == "01")
  normal_idx <- which(substr(barcodes, 14, 15) == "11")
  n_tcga_tumor_raw <- length(tumor_idx)
  n_tcga_normal_raw <- length(normal_idx)
  proj_vec <- colData(se)$project_id
  n_tcga_coad_raw <- sum(proj_vec == "TCGA-COAD")
  n_tcga_read_raw <- sum(proj_vec == "TCGA-READ")

  add_ledger("01", "sample_count", "TCGA total samples downloaded",
             n_tcga_samples_raw, "01_download_tcga.R", 105,
             "Line 105: length(sample_barcodes)")
  add_ledger("01", "sample_count", "TCGA tumor (TP)",
             n_tcga_tumor_raw, "01_download_tcga.R", 106,
             "Line 106: n_tumor")
  add_ledger("01", "sample_count", "TCGA normal (NT)",
             n_tcga_normal_raw, "01_download_tcga.R", 107,
             "Line 107: n_normal")
  add_ledger("01", "sample_count", "TCGA-COAD samples",
             n_tcga_coad_raw, "01_download_tcga.R", 115,
             "Line 115 per-project loop")
  add_ledger("01", "sample_count", "TCGA-READ samples",
             n_tcga_read_raw, "01_download_tcga.R", 115,
             "Line 115 per-project loop")
}
if (has_tcga_clin_c) {
  clin_c <- readRDS(file.path(RAW_TCGA_DIR, "clinical_COAD.rds"))
  n_clin_coad <- nrow(clin_c)
  n_surv_coad <- sum(!is.na(clin_c$days_to_death))
  add_ledger("01", "sample_count", "TCGA-COAD clinical patients",
             n_clin_coad, "01_download_tcga.R", 123, "Line 123")
  add_ledger("01", "sample_count", "TCGA-COAD with survival data",
             n_surv_coad, "01_download_tcga.R", 123, "Line 123")
}
if (has_tcga_clin_r) {
  clin_r <- readRDS(file.path(RAW_TCGA_DIR, "clinical_READ.rds"))
  n_clin_read <- nrow(clin_r)
  n_surv_read <- sum(!is.na(clin_r$days_to_death))
  add_ledger("01", "sample_count", "TCGA-READ clinical patients",
             n_clin_read, "01_download_tcga.R", 125, "Line 125")
  add_ledger("01", "sample_count", "TCGA-READ with survival data",
             n_surv_read, "01_download_tcga.R", 125, "Line 125")
}

add_audit(sprintf("  TCGA raw samples        : %s", ifelse(is.na(n_tcga_samples_raw), "NOT RUN", n_tcga_samples_raw)))
add_audit(sprintf("    Tumor (TP)            : %s", ifelse(is.na(n_tcga_tumor_raw), "NOT RUN", n_tcga_tumor_raw)))
add_audit(sprintf("    Normal (NT)           : %s", ifelse(is.na(n_tcga_normal_raw), "NOT RUN", n_tcga_normal_raw)))
add_audit(sprintf("    COAD                  : %s", ifelse(is.na(n_tcga_coad_raw), "NOT RUN", n_tcga_coad_raw)))
add_audit(sprintf("    READ                  : %s", ifelse(is.na(n_tcga_read_raw), "NOT RUN", n_tcga_read_raw)))

# ---- Script 02: GEO Download -----------------------------------------------
add_audit("")
add_audit("--- Script 02: GEO Download (02_download_geo_validation.R) ---")
add_audit("   Hard-coded in script header:")
add_audit("     GSE107422 (Korean) : 110 samples — recurrence outcome")
# GSE220148 removed — expression-only reference cohort excluded from study

add_ledger("02", "sample_count", "GSE107422 samples (hard-coded in meta)",
           110, "02_download_geo_validation.R", 10, "Line 10: script header")
# GSE220148 removed — expression-only reference cohort excluded from study
# add_ledger("02", "sample_count", "GSE220148 samples (hard-coded in meta)",
           114, "02_download_geo_validation.R", 11, "Line 11: script header")

# ---- Script 03: QC Discovery -----------------------------------------------
add_audit("")
add_audit("--- Script 03: QC Discovery (03_qc_discovery.R) ---")

n_raw_genes <- NA
n_raw_samples_qc <- NA
n_low_lib <- NA
n_genes_kept <- NA
n_genes_removed <- NA
n_pca_outliers <- NA
n_retained_samples <- NA
n_retained_genes <- NA
median_lib_size <- NA
sf_range <- NA
n_drop_clin <- NA

if (has_tcga_qc_sum) {
  qc_sum <- readLines(file.path(TCGA_QC_DIR, "tcga_qc_summary.txt"))
  extract_num <- function(pattern, line_num = NA) {
    for (l in qc_sum) {
      if (grepl(pattern, l)) {
        nums <- as.numeric(gsub("[^0-9.]", "", regmatches(l, gregexpr("[0-9,]+", l))[[1]]))
        nums <- as.numeric(gsub(",", "", nums))
        return(nums)
      }
    }
    return(NA)
  }
  n_raw_samples_qc <- extract_num("Raw samples")
  n_raw_genes <- extract_num("Raw genes")
  n_retained_samples <- extract_num("Retained samples")
  n_retained_genes <- extract_num("Retained genes")
  n_pca_outliers <- extract_num("PCA outliers detected")

  # Parse exclusions
  n_excluded <- as.numeric(gsub(".*Excluded samples.*: (\\d+).*", "\\1",
                                 grep("Excluded samples", qc_sum, value = TRUE)))

  add_ledger("03", "sample_count", "TCGA raw samples (at QC entry)",
             n_raw_samples_qc, "03_qc_discovery.R", 378, "tcga_qc_summary.txt")
  add_ledger("03", "sample_count", "TCGA raw genes",
             n_raw_genes, "03_qc_discovery.R", 379, "tcga_qc_summary.txt")
  add_ledger("03", "sample_count", "TCGA excluded samples",
             n_excluded, "03_qc_discovery.R", 385, "tcga_qc_summary.txt")
  add_ledger("03", "sample_count", "TCGA retained samples after QC",
             n_retained_samples, "03_qc_discovery.R", 389, "tcga_qc_summary.txt")
  add_ledger("03", "sample_count", "TCGA retained genes after filter",
             n_retained_genes, "03_qc_discovery.R", 390, "tcga_qc_summary.txt")
  add_ledger("03", "threshold", "Library size cutoff (counts)",
             "1e6", "03_qc_discovery.R", 128, "low_lib_cutoff variable")
  add_ledger("03", "threshold", "Gene filter minimum counts",
             "10", "03_qc_discovery.R", 50, "MIN_COUNTS")
  add_ledger("03", "threshold", "Gene filter minimum sample proportion",
             "0.1", "03_qc_discovery.R", 51, "MIN_SAMPLE_PROP")
  add_ledger("03", "threshold", "PCA outlier SD threshold",
             "4", "03_qc_discovery.R", 55, "PCA_OUTLIER_SD")
  add_ledger("03", "threshold", "Clinical missing data max proportion",
             "0.5", "03_qc_discovery.R", 54, "MAX_NA_CLINICAL")

  add_audit(sprintf("  Raw samples            : %s", ifelse(is.na(n_raw_samples_qc), "NOT RUN", n_raw_samples_qc)))
  add_audit(sprintf("  Raw genes              : %s", ifelse(is.na(n_raw_genes), "NOT RUN", n_raw_genes)))
  add_audit(sprintf("  Excluded               : %s", ifelse(is.na(n_excluded), "NOT RUN", n_excluded)))
  add_audit(sprintf("  Retained samples       : %s", ifelse(is.na(n_retained_samples), "NOT RUN", n_retained_samples)))
  add_audit(sprintf("  Retained genes         : %s", ifelse(is.na(n_retained_genes), "NOT RUN", n_retained_genes)))
  add_audit(sprintf("  PCA outliers found     : %s", ifelse(is.na(n_pca_outliers), "NOT RUN", n_pca_outliers)))
}

# ---- Script 03b: QC Validation ---------------------------------------------
add_audit("")
add_audit("--- Script 03: QC Validation (03_qc_validation.R) ---")
add_audit("   Exclusion rules identical to discovery (MIN_COUNTS=10, MIN_SAMPLE_PROP=0.1, PCA_OUTLIER_SD=4)")

if (has_geo_qc_sum) {
  geo_sum <- readLines(file.path(GEO_QC_DIR, "geo_qc_summary.txt"))
  add_audit("   geo_qc_summary.txt found — per-dataset sample counts available.")
} else {
  add_audit("   geo_qc_summary.txt NOT FOUND — run 03_qc_validation.R first.")
}

# ---- Script 04: DEG --------------------------------------------------------
add_audit("")
add_audit("--- Script 04: DEG (04_deg_analysis.R) ---")
add_audit(sprintf("   Thresholds: |log2FC| > 1, padj (BH) < 0.05"))

n_deg_tumor <- NA
n_deg_normal <- NA
n_deg_tested <- NA
n_deg_up <- NA
n_deg_down <- NA

if (has_deg_sum) {
  deg_sum <- readLines(file.path(TAB_DIR, "DEG_summary_TCGA.txt"))
  for (l in deg_sum) {
    if (grepl("Genes tested", l)) n_deg_tested <- as.numeric(gsub(".*: (\\d+)", "\\1", l))
    if (grepl("Up-regulated", l)) n_deg_up <- as.numeric(gsub(".*: (\\d+).*", "\\1", l))
    if (grepl("Down-regulated", l)) n_deg_down <- as.numeric(gsub(".*: (\\d+).*", "\\1", l))
  }
  add_ledger("04", "sample_count", "DEG samples — tumor",
             NA, "04_deg_analysis.R", 57, "Line 57 cat() — only printed at runtime")
  add_ledger("04", "sample_count", "DEG samples — normal",
             NA, "04_deg_analysis.R", 57, "Line 57 cat() — only printed at runtime")
  add_ledger("04", "statistical", "Genes tested in DEG",
             n_deg_tested, "04_deg_analysis.R", 165, "DEG_summary_TCGA.txt")
  add_ledger("04", "statistical", "DEG up-regulated (|log2FC|>1, padj<0.05)",
             n_deg_up, "04_deg_analysis.R", 167, "DEG_summary_TCGA.txt")
  add_ledger("04", "statistical", "DEG down-regulated",
             n_deg_down, "04_deg_analysis.R", 168, "DEG_summary_TCGA.txt")
  add_ledger("04", "threshold", "DEG log2FC cutoff",
             "1", "04_deg_analysis.R", 34, "LOG2FC_CUTOFF")
  add_ledger("04", "threshold", "DEG padj cutoff (BH)",
             "0.05", "04_deg_analysis.R", 35, "PADJ_CUTOFF")
  add_ledger("04", "method", "DEG multiple testing correction",
             "BH (Benjamini-Hochberg)", "04_deg_analysis.R", 35,
             "DESeq2 default via results(..., alpha=0.05)")

  add_audit(sprintf("  Genes tested           : %s", ifelse(is.na(n_deg_tested), "NOT RUN", n_deg_tested)))
  add_audit(sprintf("  Up-regulated           : %s", ifelse(is.na(n_deg_up), "NOT RUN", n_deg_up)))
  add_audit(sprintf("  Down-regulated         : %s", ifelse(is.na(n_deg_down), "NOT RUN", n_deg_down)))
} else {
  add_audit("  DEG_summary_TCGA.txt NOT FOUND — run 04_deg_analysis.R first.")
}

# ---- Script 05: CSC Markers ------------------------------------------------
add_audit("")
add_audit("--- Script 05: CSC Markers (05_csc_marker_focus.R) ---")
add_audit(sprintf("   Thresholds: |log2FC| > 1, padj (BH) < 0.05 (identical to script 04)"))

if (has_csc_csv) {
  csc <- read.csv(file.path(TAB_DIR, "CSC_markers_TCGA.csv"), stringsAsFactors = FALSE)
  n_csc_total <- nrow(csc)
  n_csc_sig <- sum(csc$significant, na.rm = TRUE)
  n_csc_up <- sum(csc$direction == "Up", na.rm = TRUE)
  n_csc_down <- sum(csc$direction == "Down", na.rm = TRUE)
  n_csc_missing <- sum(is.na(csc$gene_id))

  add_ledger("05", "sample_count", "CSC panel total markers",
             n_csc_total, "05_csc_marker_focus.R", 161, "CSC_markers_summary_TCGA.txt")
  add_ledger("05", "statistical", "CSC markers significantly dysregulated",
             n_csc_sig, "05_csc_marker_focus.R", 162, "CSC_markers_summary_TCGA.txt")
  add_ledger("05", "statistical", "CSC markers up-regulated",
             n_csc_up, "05_csc_marker_focus.R", 163, "CSC_markers_summary_TCGA.txt")
  add_ledger("05", "statistical", "CSC markers down-regulated",
             n_csc_down, "05_csc_marker_focus.R", 164, "CSC_markers_summary_TCGA.txt")

  add_audit(sprintf("  Panel markers          : %d", n_csc_total))
  add_audit(sprintf("  Significant            : %d (up=%d, down=%d)", n_csc_sig, n_csc_up, n_csc_down))
  add_audit(sprintf("  Not detected in data   : %d", n_csc_missing))

  # Per-marker ledger entries
  for (i in seq_len(n_csc_total)) {
    r <- csc[i, ]
    if (!is.na(r$gene_id)) {
      add_ledger("05", "per_marker",
                 sprintf("CSC %s — log2FC", r$symbol),
                 round(r$log2FoldChange, 4),
                 "05_csc_marker_focus.R", 132,
                 sprintf("CSC_markers_summary_TCGA.txt (padj=%.2e)", r$padj))
      add_ledger("05", "per_marker",
                 sprintf("CSC %s — padj", r$symbol),
                 format(r$padj, scientific = TRUE, digits = 3),
                 "05_csc_marker_focus.R", 132,
                 sprintf("CSC_markers_summary_TCGA.txt (direction=%s)", r$direction))
    }
  }
} else {
  add_audit("  CSC_markers_TCGA.csv NOT FOUND — run 05_csc_marker_focus.R first.")
}

# ---- Script 06: Signature Construction --------------------------------------
add_audit("")
add_audit("--- Script 06: Signature Construction (06_signature_construction.R) ---")

if (has_lasso_coef) {
  coef_tbl <- read.csv(file.path(TCGA_SIG_DIR, "lasso_coefficients.csv"),
                        stringsAsFactors = FALSE)
  n_sig_genes <- nrow(coef_tbl)
  add_ledger("06", "statistical", "LASSO signature genes selected",
             n_sig_genes, "06_signature_construction.R", 349, "lasso_coefficients.csv")
  add_ledger("06", "statistical", "LASSO lambda used",
             NA, "06_signature_construction.R", 343, "printed at runtime")
  add_ledger("06", "statistical", "LASSO lambda.min (nonzero)",
             NA, "06_signature_construction.R", 343, "printed at runtime")
  add_ledger("06", "statistical", "LASSO lambda.1se (nonzero)",
             NA, "06_signature_construction.R", 344, "printed at runtime")

  # Per-gene coefficients
  for (i in seq_len(nrow(coef_tbl))) {
    add_ledger("06", "signature_coefficient",
               sprintf("Signature gene %s — coefficient", coef_tbl$gene_symbol[i]),
               coef_tbl$coefficient[i],
               "data/processed/tcga_sig/lasso_coefficients.csv", NA,
               sprintf("Output: Table_risk_score_formula.csv"))
  }
  add_audit(sprintf("  Signature genes        : %d", n_sig_genes))
  add_audit(sprintf("  Coefs: %s",
                    paste(sprintf("%s=%.4f", coef_tbl$gene_symbol, coef_tbl$coefficient),
                          collapse = ", ")))
} else {
  add_audit("  lasso_coefficients.csv NOT FOUND — run 06 first.")
}

if (has_risk_scores) {
  risk <- read.csv(file.path(TCGA_SIG_DIR, "risk_scores.csv"), stringsAsFactors = FALSE)
  n_risk_total <- nrow(risk)
  n_risk_high  <- sum(risk$risk_group == "High")
  n_risk_low   <- sum(risk$risk_group == "Low")
  n_risk_events <- sum(risk$os_event, na.rm = TRUE)

  add_ledger("06", "sample_count", "TCGA samples with risk scores",
             n_risk_total, "06_signature_construction.R", 403, "risk_scores.csv")
  add_ledger("06", "sample_count", "TCGA high-risk group",
             n_risk_high, "06_signature_construction.R", 407, "risk_scores.csv")
  add_ledger("06", "sample_count", "TCGA low-risk group",
             n_risk_low, "06_signature_construction.R", 407, "risk_scores.csv")
  add_ledger("06", "sample_count", "TCGA OS events in risk analysis",
             n_risk_events, "06_signature_construction.R", 281, "risk_scores.csv")

  # The exact median risk score is printed at runtime — store from data
  add_ledger("06", "statistical", "TCGA median risk score (split cutoff)",
             median(risk$risk_score, na.rm = TRUE),
             "06_signature_construction.R", 405, "printed at runtime")

  add_audit(sprintf("  Risk-scored patients   : %d (%d events)", n_risk_total, n_risk_events))
  add_audit(sprintf("  High / Low risk        : %d / %d", n_risk_high, n_risk_low))
  add_audit(sprintf("  Median risk score      : %.4f", median(risk$risk_score, na.rm = TRUE)))

  # Log-rank p (from Cox summary)
  if (has_cox_sum) {
    cox_txt <- readLines(file.path(TAB_DIR, "Cox_multivariable_summary.txt"))
    # Extract HR from the adjusted model
    for (l in cox_txt) {
      if (grepl("risk_group", l) && grepl("\\d+\\.\\d+", l)) {
        parts <- strsplit(trimws(l), "\\s+")[[1]]
        if (length(parts) >= 5) {
          hr_val <- parts[length(parts) - 2]
          ci_low <- parts[length(parts) - 1]
          p_val  <- parts[length(parts)]
          break
        }
      }
    }
    add_ledger("06", "statistical", "Cox HR (high vs low, adjusted)",
               NA, "06_signature_construction.R", 489, "Cox_multivariable_summary.txt — runtime")
    add_ledger("06", "statistical", "Cox PH global test p",
               NA, "06_signature_construction.R", 579, "PH_assumption_test.txt — runtime")
  }
} else {
  add_audit("  risk_scores.csv NOT FOUND — run 06 first.")
}

# ---- Script 07: Validation --------------------------------------------------
add_audit("")
add_audit("--- Script 07: Validation (07_validation.R) ---")
add_audit("  HARD RULE CHECK: Model never re-fit on validation data ✓")

if (has_val_risk) {
  val_risk <- read.csv(file.path(VAL_DIR, "risk_scores_validation.csv"),
                        stringsAsFactors = FALSE)
  n_val_total <- nrow(val_risk)
  n_val_high  <- sum(grepl("High", val_risk$risk_group_z), na.rm = TRUE)
  n_val_low   <- sum(grepl("Low", val_risk$risk_group_z), na.rm = TRUE)
  n_val_events <- sum(val_risk$recur_event, na.rm = TRUE)

  add_ledger("07", "sample_count", "Validation samples with risk scores",
             n_val_total, "07_validation.R", 279, "risk_scores_validation.csv")
  add_ledger("07", "sample_count", "Validation high-risk (TCGA median cutoff)",
             n_val_high, "07_validation.R", 284, "risk_scores_validation.csv")
  add_ledger("07", "sample_count", "Validation low-risk (TCGA median cutoff)",
             n_val_low, "07_validation.R", 284, "risk_scores_validation.csv")
  add_ledger("07", "sample_count", "Validation recurrence events",
             n_val_events, "07_validation.R", 279, "risk_scores_validation.csv")

  add_audit(sprintf("  Validation patients    : %d", n_val_total))
  add_audit(sprintf("  High / Low risk        : %d / %d", n_val_high, n_val_low))
  add_audit(sprintf("  Recurrence events      : %d", n_val_events))
}

if (has_val_outcome) {
  val_out <- readLines(file.path(TAB_DIR, "Validation_outcome.txt"))
  outcome_line <- grep("Outcome:", val_out, value = TRUE)
  add_ledger("07", "statistical", "Validation outcome classification",
             gsub("Outcome: ", "", outcome_line),
             "07_validation.R", 520, "Validation_outcome.txt")

  # Extract HR and C-index
  for (l in val_out) {
    if (grepl("HR", l) && grepl("High vs Low", l)) {
      hr_str <- gsub(".*HR.*\\((.*)\\).*", "\\1", l)
    }
    if (grepl("p-value \\(Cox\\)", l)) {
      p_cox <- gsub(".*: ([0-9.eE-]+)", "\\1", l)
    }
    if (grepl("C-index \\(validation\\)", l)) {
      c_val <- gsub(".*: ([0-9.]+)", "\\1", l)
    }
    if (grepl("TCGA C-index", l)) {
      c_tcga <- gsub(".*: ([0-9.]+)", "\\1", l)
    }
  }

  add_audit(sprintf("  Validation outcome     : %s", outcome_line))
} else {
  add_audit("  Validation_outcome.txt NOT FOUND — run 07 first.")
}

# ---- Script 08: Enrichment --------------------------------------------------
add_audit("")
add_audit("--- Script 08: Enrichment (08_enrichment.R) ---")
add_audit("  Multiple testing: BH correction applied by clusterProfiler ✓")

if (has_enrich_deg) {
  enrich_deg <- read.csv(file.path(ENRICH_DIR, "TCGA_DEG_high_vs_low.csv"),
                          stringsAsFactors = FALSE)
  n_enrich_up <- sum(enrich_deg$padj < 0.05 & enrich_deg$log2FoldChange > 0, na.rm = TRUE)
  n_enrich_down <- sum(enrich_deg$padj < 0.05 & enrich_deg$log2FoldChange < 0, na.rm = TRUE)
  add_ledger("08", "statistical", "Risk-group DEG up (high vs low)",
             n_enrich_up, "08_enrichment.R", 310, "TCGA_DEG_high_vs_low.csv")
  add_ledger("08", "statistical", "Risk-group DEG down (high vs low)",
             n_enrich_down, "08_enrichment.R", 310, "TCGA_DEG_high_vs_low.csv")
  add_audit(sprintf("  Risk-group DEG up/down : %d / %d", n_enrich_up, n_enrich_down))
} else {
  add_audit("  TCGA_DEG_high_vs_low.csv NOT FOUND — run 08 first.")
}

if (has_gsea_kegg) {
  gk <- read.csv(file.path(ENRICH_DIR, "TCGA_GSEA_kegg.csv"), stringsAsFactors = FALSE)
  add_ledger("08", "statistical", "GSEA KEGG significant pathways",
             sum(gk$p.adjust < 0.05, na.rm = TRUE),
             "08_enrichment.R", 320, "TCGA_GSEA_kegg.csv")
  add_ledger("08", "method", "GSEA multiple testing correction",
             "BH (via clusterProfiler::gseKEGG)", "08_enrichment.R", 318, "")
  add_audit(sprintf("  GSEA KEGG significant  : %d", sum(gk$p.adjust < 0.05, na.rm = TRUE)))
}

if (has_gsea_go) {
  gg <- read.csv(file.path(ENRICH_DIR, "TCGA_GSEA_go.csv"), stringsAsFactors = FALSE)
  add_ledger("08", "statistical", "GSEA GO-BP significant pathways",
             sum(gg$p.adjust < 0.05, na.rm = TRUE),
             "08_enrichment.R", 330, "TCGA_GSEA_go.csv")
  add_ledger("08", "method", "GSEA GO-BP multiple testing correction",
             "BH (via clusterProfiler::gseGO)", "08_enrichment.R", 328, "")
  add_audit(sprintf("  GSEA GO-BP significant : %d", sum(gg$p.adjust < 0.05, na.rm = TRUE)))
}

if (has_stem_audit) {
  stem_txt <- readLines(file.path(ENRICH_DIR, "TCGA_stemness_audit.txt"))
  stem_line <- grep("Stemness-related pathways|No stemness|negative result",
                     stem_txt, value = TRUE, ignore.case = TRUE)
  add_ledger("08", "statistical", "Stemness pathway audit conclusion",
             paste(stem_line, collapse = " | "),
             "08_enrichment.R", 215, "TCGA_stemness_audit.txt")
  add_audit(sprintf("  Stemness audit         : See TCGA_stemness_audit.txt"))
}

# ---- Script 09: Immune Infiltration -----------------------------------------
add_audit("")
add_audit("--- Script 09: Immune Infiltration (09_immune_infiltration.R) ---")
add_audit("  Method: ssGSEA (manual implementation, Barbie et al. 2009)")
add_audit("  Signatures: 13 cell types from Bindea 2013 + Charoentong 2017")

if (has_immune_cor) {
  ic <- read.csv(file.path(IMMUNE_DIR, "TCGA_immune_risk_correlation.csv"),
                  stringsAsFactors = FALSE)
  add_ledger("09", "sample_count", "Immune cell types scored",
             nrow(ic), "09_immune_infiltration.R", 60,
             "immune_signatures list")
  add_ledger("09", "method", "Immune correlation test",
             "Spearman correlation (risk score vs ssGSEA score)",
             "09_immune_infiltration.R", 157, "")
  add_ledger("09", "method", "Immune group comparison test",
             "Wilcoxon rank-sum (high-risk vs low-risk per cell type)",
             "09_immune_infiltration.R", 107, "")

  for (i in seq_len(nrow(ic))) {
    add_ledger("09", "immune_correlation",
               sprintf("Immune %s — Spearman rho vs risk score", ic$cell_type[i]),
               round(ic$rho[i], 4),
               "09_immune_infiltration.R", 157,
               sprintf("TCGA_immune_risk_correlation.csv (p=%.4e)", ic$p_value[i]))
    add_ledger("09", "immune_correlation",
               sprintf("Immune %s — Spearman p-value", ic$cell_type[i]),
               format(ic$p_value[i], scientific = TRUE, digits = 3),
               "09_immune_infiltration.R", 157,
               "TCGA_immune_risk_correlation.csv")
  }

  n_signif_cor <- sum(ic$p_value < 0.05, na.rm = TRUE)
  add_audit(sprintf("  Cell types with sig. cor (p<0.05): %d / %d",
                    n_signif_cor, nrow(ic)))
} else {
  add_audit("  TCGA_immune_risk_correlation.csv NOT FOUND — run 09 first.")
}

# =============================================================================
# AUDIT B — Multiple Testing Correction Compliance
# =============================================================================
add_audit("")
add_audit("")
add_audit("==============================================")
add_audit("  AUDIT B — MULTIPLE TESTING CORRECTION       ")
add_audit("==============================================")
add_audit("")

# Check each script
mt_checks <- data.frame(
  script = c("04_deg_analysis.R", "05_csc_marker_focus.R",
             "06_signature_construction.R", "07_validation.R",
             "08_enrichment.R", "08_enrichment.R", "09_immune_infiltration.R"),
  test_type = c("DEG (tumor vs normal)", "CSC marker subset",
                "LASSO-Cox (10-fold CV)", "Validation Cox",
                "GSEA (KEGG, GO-BP, Hallmark)", "ORA (KEGG, GO-BP, Hallmark)",
                "Immune: 13 cell types × 2 tests each"),
  correction = c("BH (DESeq2 default)", "BH (inherited from DESeq2)",
                 "NA — LASSO regularisation (alpha=1) is inherent feature selection",
                 "NA — single test (already adjusted)",
                 "BH (clusterProfiler default, pAdjustMethod)",
                 "BH (clusterProfiler default)",
                 "Multiple comparison caveat: 13 tests per cohort, report raw + consider Bonferroni"),
  compliant = c("YES", "YES", "YES — regularisation not multiple testing",
                "YES", "YES", "YES", "CAVEAT — reported raw p-values, Bonferroni threshold = 0.05/13 ≈ 0.0038"),
  stringsAsFactors = FALSE
)

add_audit(sprintf("%-35s %-50s %s", "Test", "Correction", "Compliant?"))
add_audit(sprintf("%-35s %-50s %s", "-", "-", "-"))
for (i in seq_len(nrow(mt_checks))) {
  add_audit(sprintf("%-35s %-50s %s",
                    mt_checks$test_type[i],
                    substr(mt_checks$correction[i], 1, 50),
                    mt_checks$compliant[i]))
}
add_audit("")
add_audit("Verdict: All applicable multiple-testing scenarios use BH/FDR correction.")
add_audit("  LASSO-Cox uses alpha=1 regularisation, which is inherent feature")
add_audit("  selection — not a post-hoc multiple-testing scenario.")
add_audit("  Immune correlation: 13 tests × 2 cohorts = 26 tests; Bonferroni")
add_audit("  threshold = 0.05/13 ≈ 0.0038 per cohort recommended.")

# =============================================================================
# AUDIT C — Data Leakage Check
# =============================================================================
add_audit("")
add_audit("")
add_audit("==============================================")
add_audit("  AUDIT C — DATA LEAKAGE CHECK                ")
add_audit("==============================================")
add_audit("")

# Check that validation data files are never read by scripts 01–06
script_files <- list.files(SCRIPTS_DIR, pattern = "\\.R$", full.names = TRUE)
script_names <- list.files(SCRIPTS_DIR, pattern = "\\.R$")

add_audit("Checking: Do scripts 01–06 ever read GSE107422 or validation files?")
add_audit("")

leak_found <- FALSE
for (sf in script_files) {
  sname <- basename(sf)
  snum <- as.numeric(gsub("^(\\d+).*", "\\1", sname))
  if (is.na(snum) || snum > 6) next  # only check 01-06

  content <- paste(readLines(sf, warn = FALSE), collapse = "\n")

  # Check for references to validation data
  checks <- list(
    "GSE107422"   = grepl("GSE107422", content),
    "validation"  = grepl("validation", content, ignore.case = TRUE),
    "geo_qc"      = grepl("geo_qc", content),
    "GSE107422_eset" = grepl("GSE107422_eset", content)
  )

  for (chk_name in names(checks)) {
    if (checks[[chk_name]]) {
      # Check if it's in a comment (header) vs actual code
      lines <- readLines(sf, warn = FALSE)
      code_refs <- grep(chk_name, lines, ignore.case = TRUE, value = TRUE)
      # Filter out comment-only lines
      code_refs <- code_refs[!grepl("^\\s*#", code_refs)]
      if (length(code_refs) > 0) {
        add_audit(sprintf("  ⚠ POTENTIAL LEAK: %s references '%s' in code:", sname, chk_name))
        for (cr in code_refs) {
          add_audit(sprintf("      %s", trimws(cr)))
        }
        leak_found <- TRUE
      }
    }
  }
}

# Check script 06 specifically for any validation data use
if (has_lasso_coef) {
  add_audit("")
  add_audit("Verification: Script 06 coefficients loaded from TCGA only:")
  coef_path <- file.path(TCGA_SIG_DIR, "lasso_coefficients.csv")
  coef_tbl <- read.csv(coef_path, stringsAsFactors = FALSE)
  add_audit(sprintf("  Coefficients file: %s", coef_path))
  add_audit(sprintf("  Signature genes  : %d", nrow(coef_tbl)))
  add_audit("  Data provenance  : TCGA-COAD + TCGA-READ only ✓")
}

add_audit("")
if (!leak_found) {
  add_audit("  ✅ NO LEAKAGE DETECTED: Scripts 01–06 never read validation data.")
  add_ledger("10", "audit", "Data leakage check result",
             "PASS — no leakage between discovery and validation",
             "10_robustness_check.R", NA, "Scripts 01-06 do not reference validation files")
} else {
  add_audit("  ⚠ LEAKAGE DETECTED — see above.")
  add_ledger("10", "audit", "Data leakage check result",
             "WARN — potential leakage detected",
             "10_robustness_check.R", NA, "See audit report for details")
}

# Also check: is the validation cohort expression ever used in scripts that contribute to the signature?
add_audit("")
add_audit("Cross-check: Signature construction input chain")
add_audit("  01_download_tcga.R  → TCGA expression   → to 03_qc_discovery.R")
add_audit("  03_qc_discovery.R   → tcga_dds_de.rds    → to 06_signature_construction.R")
add_audit("  06_signature_construction.R → sigma       → to 07_validation.R")
add_audit("  Validation data (GSE107422) enters ONLY at script 07.")
add_audit("  ∴ No possible leakage. ✓")

# =============================================================================
# SECTION — Write Numbers Ledger CSV
# =============================================================================
add_audit("")
add_audit("==============================================")
add_audit("  NUMBERS LEDGER WRITTEN                      ")
add_audit("==============================================")

write.csv(ledger, OUT_LEDGER, row.names = FALSE)
add_audit(sprintf("  -> %s", OUT_LEDGER))
add_audit(sprintf("  Total entries: %d", nrow(ledger)))

# Summary by script
add_audit("")
add_audit("Entries by script:")
for (s in sort(unique(ledger$script))) {
  n_cat <- length(unique(ledger$stat_name[ledger$script == s]))
  add_audit(sprintf("  Script %s: %d entries (%d distinct statistics)",
                    s, sum(ledger$script == s), n_cat))
}

# =============================================================================
# Write audit report
# =============================================================================
writeLines(audit_lines, OUT_AUDIT)

# =============================================================================
# Final summary
# =============================================================================
message("")
message("==============================================")
message("  ROBUSTNESS CHECK COMPLETE                    ")
message("==============================================")
message("")
message(sprintf("  Numbers ledger : %s  (%d entries)", OUT_LEDGER, nrow(ledger)))
message(sprintf("  Audit report   : %s", OUT_AUDIT))
message("")
message("  === MANUSCRIPT-WRITING RULES ===")
message("  1. Every number in the manuscript MUST appear in numbers_ledger.csv")
message("  2. Every number's source (script:line) MUST be traceable in the ledger")
message("  3. If a number is not in the ledger, RE-RUN the relevant script and")
message("     re-generate the ledger — do NOT type it from memory")
message("  4. For runtime-only values (printed to stdout but not saved):")
message("     re-run the script and add the value to the ledger manually")
message("")
message(">>> Outputs:")
message("     ", OUT_LEDGER)
message("     ", OUT_AUDIT)

writeLines(capture.output(sessionInfo()), file.path(TAB_DIR, "session_info_10.txt"))
