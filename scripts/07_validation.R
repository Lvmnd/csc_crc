#!/usr/bin/env Rscript
# =============================================================================
# 07_validation.R — Independent validation of TCGA CSC signature in GSE107422
#
# HARD RULES (pre-registered):
#   - NEVER re-fit the model. Apply the exact TCGA-derived coefficients.
#   - Report validation outcome honestly: success, partial, or fail.
#   - Do NOT suppress or reframe unfavourable results.
#
# Workflow:
#   1. Load TCGA signature coefficients (fixed, from script 06)
#   2. Load TCGA expression data (to compute z-score reference params)
#   3. Load GSE107422 expression + clinical data (post-QC)
#   4. Map signature genes by symbol across platforms
#   5. Cross-platform normalisation: z-score each gene within each cohort
#   6. Compute risk scores using fixed coefficients
#   7. Stratify: primary = TCGA median; sensitivity = validation median
#   8. KM + log-rank for recurrence-free survival
#   9. Cox PH + C-index (both cohorts)
#  10. Validation outcome classification
#  11. Sensitivity: tertile-based stratification
#
# Inputs:
#   data/processed/tcga_sig/lasso_coefficients.csv   — fixed signature
#   data/processed/tcga_qc/tcga_dds_de.rds             — TCGA norm counts
#   data/processed/geo_qc/GSE107422/*.rds              — validation expression
#   data/raw/geo/GSE107422/GSE107422_eset.rds          — validation clinical
#
# Outputs:
#   data/processed/validation/
#   ├── risk_scores_validation.csv         — per-patient risk scores
#   ├── cross_platform_params.rds          — z-score ref params from TCGA
#   └── tcga_risk_scores_ref.rds           — TCGA risk scores for C-index
#
#   results/tables/
#   ├── Validation_cox_summary.txt         — HR + 95% CI + p (validation)
#   ├── Validation_outcome.txt             — (success/partial/fail) + evidence
#   ├── Validation_Cindex_comparison.txt   — C-index TCGA vs validation
#   └── Validation_tertile_sensitivity.txt — tertile-based results
#
#   results/figures/validation/
#   ├── 01_km_curve_median.png
#   ├── 02_km_curve_tertiles.png
#   ├── 03_forest_plot.png
#   ├── 04_cox_zph_plot.png
#   └── 05_risk_score_distribution.png
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(survival)
  library(survminer)
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(pROC)
  library(tidyr)
  library(org.Hs.eg.db)
})

# -- Directories -------------------------------------------------------------
TCGA_SIG_DIR <- file.path("data", "processed", "tcga_sig")
TCGA_QC_DIR  <- file.path("data", "processed", "tcga_qc")
GEO_QC_DIR   <- file.path("data", "processed", "geo_qc", "GSE107422")
GEO_RAW_DIR  <- file.path("data", "raw", "geo", "GSE107422")
OUT_DIR      <- file.path("data", "processed", "validation")
TAB_DIR      <- file.path("results", "tables")
FIG_DIR      <- file.path("results", "figures", "validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load TCGA signature coefficients
# =============================================================================
message(">>> Loading TCGA signature coefficients ...")

coef_file <- file.path(TCGA_SIG_DIR, "lasso_coefficients.csv")
if (!file.exists(coef_file)) {
  stop("Coefficient file not found: ", coef_file,
       "\n  Run scripts/06_signature_construction.R first.")
}

coef_tbl <- read.csv(coef_file, stringsAsFactors = FALSE)
message(sprintf("   Signature genes: %d", nrow(coef_tbl)))
message(sprintf("   Coefficients loaded from: %s", coef_file))

# Also load TCGA risk scores and survival data (for C-index comparison)
tcga_risk_file <- file.path(TCGA_SIG_DIR, "risk_scores.csv")
tcga_risk <- read.csv(tcga_risk_file, stringsAsFactors = FALSE)
message(sprintf("   TCGA risk scores loaded: %d patients", nrow(tcga_risk)))

# =============================================================================
# 2. Load TCGA expression data (reference for z-score params)
# =============================================================================
message(">>> Loading TCGA expression data for z-score reference ...")

dds_file <- file.path(TCGA_QC_DIR, "tcga_dds_de.rds")
if (!file.exists(dds_file)) {
  stop("TCGA dds not found: ", dds_file)
}
dds <- readRDS(dds_file)
tcga_norm <- counts(dds, normalized = TRUE)

# Log2-transform (same as in script 06)
tcga_log <- log2(tcga_norm + 1)

# Gene symbol mapping (ENSEMBL → symbol via org.Hs.eg.db)
ensembls <- gsub("\\..*$", "", rownames(dds))
symbol_map <- suppressMessages(mapIds(org.Hs.eg.db, keys = ensembls, keytype = "ENSEMBL",
                     column = "SYMBOL", multiVals = "first"))
gene_map <- data.frame(
  gene_id   = rownames(dds),
  ensembl   = ensembls,
  gene_name = ifelse(is.na(symbol_map), rownames(dds), symbol_map),
  stringsAsFactors = FALSE
)

message(sprintf("   TCGA expression matrix: %d genes x %d samples",
                nrow(tcga_log), ncol(tcga_log)))

# =============================================================================
# 3. Load GSE107422 expression data (post-QC)
# =============================================================================
message(">>> Loading GSE107422 expression data ...")

# Determine what QC output exists
counts_file <- file.path(GEO_QC_DIR, "counts_filtered.rds")
tpm_file    <- file.path(GEO_QC_DIR, "tpm_filtered.rds")
log_file    <- file.path(GEO_QC_DIR, "log_normalized.rds")
log2_file   <- file.path(GEO_QC_DIR, "log2_tpm.rds")

if (file.exists(counts_file)) {
  message("   Found count data from GSE107422")
  val_counts <- readRDS(counts_file)
  # Normalise + log2 — same pipeline as TCGA
  dds_val <- DESeqDataSetFromMatrix(
    countData = val_counts,
    colData = data.frame(row.names = colnames(val_counts), dummy = 1),
    design = ~ 1
  )
  dds_val <- estimateSizeFactors(dds_val)
  val_norm <- counts(dds_val, normalized = TRUE)
  val_log <- log2(val_norm + 1)
  data_type <- "counts_deseq2_normalized"
  rm(dds_val)
} else if (file.exists(tpm_file)) {
  message("   Found TPM data from GSE107422")
  val_tpm <- readRDS(tpm_file)
  val_log <- log2(val_tpm + 1)
  data_type <- "tpm_log2"
} else if (file.exists(log_file)) {
  message("   Found log-normalized data from GSE107422")
  val_log <- readRDS(log_file)
  data_type <- "pre_normalized_log"
} else if (file.exists(log2_file)) {
  message("   Found log2-TPM data from GSE107422")
  val_log <- readRDS(log2_file)
  data_type <- "log2_tpm"
} else {
  stop("No processed expression data found for GSE107422 in ", GEO_QC_DIR,
       "\n  Run scripts/03_qc_validation.R first.")
}

message(sprintf("   Validation expression: %d genes x %d samples (type: %s)",
                nrow(val_log), ncol(val_log), data_type))

# =============================================================================
# 4. Gene matching between TCGA and validation
# =============================================================================
message(">>> Matching signature genes between cohorts ...")

# -- Map validation gene IDs to symbols ---------------------------------------
# GEO expression data may have various row name formats:
#   - Gene symbols (e.g., "PROM1")
#   - Ensembl IDs (e.g., "ENSG00000000003")
#   - Probe IDs (e.g., "ILMN_12345")
# We try to detect the format and map to symbols.

val_gene_ids <- rownames(val_log)

# Heuristic: if most start with ENSG, they're Ensembl IDs
is_ensembl <- mean(grepl("^ENSG", val_gene_ids)) > 0.5
# If most are all-caps alphanumeric (symbols), treat as symbols
is_symbol  <- mean(grepl("^[A-Z0-9]+$", val_gene_ids) & !grepl("^ENSG", val_gene_ids)) > 0.5

if (is_ensembl) {
  message("   Validation gene IDs detected as: Ensembl")
  val_genes_df <- data.frame(
    val_id   = val_gene_ids,
    val_name = val_gene_ids,  # fallback
    stringsAsFactors = FALSE
  )
  # Try to match via TCGA gene_map (Ensembl ID → symbol)
  tcga_map <- gene_map[!duplicated(gene_map$gene_id), ]
  rownames(tcga_map) <- tcga_map$gene_id
  matched <- intersect(val_gene_ids, tcga_map$gene_id)
  val_genes_df$tcga_symbol <- tcga_map[matched, "gene_name"]
  message(sprintf("   Ensembl IDs matched to symbols: %d", length(matched)))
} else if (is_symbol) {
  message("   Validation gene IDs detected as: Symbols")
  val_genes_df <- data.frame(
    val_id      = val_gene_ids,
    val_name    = val_gene_ids,
    tcga_symbol = val_gene_ids,
    stringsAsFactors = FALSE
  )
} else {
  message("   Validation gene IDs format: other (non-Ensembl, non-standard symbols)")
  message("   Attempting RefSeq → Symbol mapping ...")
  # Extract RefSeq accessions from gi|...|ref|NM_...| format
  refseq_ids <- gsub("^.*\\|ref\\|([^|]+)\\|.*$", "\\1", val_gene_ids)
  refseq_ids <- gsub("\\..*$", "", refseq_ids)  # strip version
  refseq_map <- suppressMessages(mapIds(org.Hs.eg.db, keys = refseq_ids, keytype = "REFSEQ",
                       column = "SYMBOL", multiVals = "first"))
  n_mapped <- sum(!is.na(refseq_map))
  message(sprintf("   RefSeq → Symbol mapped: %d / %d", n_mapped, length(refseq_ids)))
  val_genes_df <- data.frame(
    val_id      = val_gene_ids,
    val_name    = ifelse(is.na(refseq_map), val_gene_ids, refseq_map),
    tcga_symbol = ifelse(is.na(refseq_map), val_gene_ids, refseq_map),
    stringsAsFactors = FALSE
  )
}

# For the signature genes, find matching rows in validation
sig_symbols <- coef_tbl$gene_symbol
sig_symbols_upper <- toupper(sig_symbols)

val_symbols_upper <- toupper(val_genes_df$tcga_symbol)

# For each signature gene, find the best match
gene_match <- data.frame(
  sig_gene      = sig_symbols,
  sig_coefficient = coef_tbl$coefficient,
  val_row       = NA_integer_,
  val_id        = NA_character_,
  matched       = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(gene_match))) {
  hit <- which(val_symbols_upper == sig_symbols_upper[i])
  if (length(hit) > 0) {
    gene_match$val_row[i]  <- hit[1]  # take first match
    gene_match$val_id[i]   <- val_genes_df$val_id[hit[1]]
    gene_match$matched[i]  <- TRUE
  }
}

n_matched <- sum(gene_match$matched)
message(sprintf("   Signature genes matched: %d / %d", n_matched, nrow(gene_match)))

if (n_matched == 0) {
  stop("No signature genes could be matched between TCGA and validation cohort.")
}

if (n_matched < nrow(gene_match) * 0.5) {
  warning(sprintf("  Only %d / %d signature genes matched. Signature may be unstable.",
                  n_matched, nrow(gene_match)))
}

# =============================================================================
# 5. Cross-platform normalisation: z-score per gene within each cohort
# =============================================================================
message(">>> Cross-platform normalisation (z-score per gene within cohort) ...")
message("")
message("   Method: For each signature gene, compute z-score within each")
message("   cohort independently: z = (x - mean_cohort) / sd_cohort.")
message("   This renders each gene distribution to mean=0, SD=1 in both")
message("   cohorts, making the TCGA-estimated coefficients β directly")
message("   transferable across platforms.")
message("")
message("   Limitation (transparently reported): The coefficients were")
message("   originally estimated on log2(norm_counts+1) data in TCGA.")
message("   Z-score transformation preserves relative ordering but alters")
message("   the absolute scale, which may affect the coefficient magnitude.")
message("   Alternative raw-scale results are reported for comparison.")

# Extract matched signature expression from TCGA
matched_ids <- gene_map$gene_id[match(sig_symbols, gene_map$gene_name)]
matched_ids <- matched_ids[!is.na(matched_ids)]

# TCGA: log2 expression for matched signature genes
tcga_sig_log <- tcga_log[intersect(matched_ids, rownames(tcga_log)), , drop = FALSE]
# Also try matching by symbol directly
tcga_sig_log2 <- tcga_log[which(gene_map$gene_name %in% sig_symbols), , drop = FALSE]

if (nrow(tcga_sig_log2) == 0) {
  stop("No signature genes found in TCGA expression matrix.")
}

# Reorder to match coefficient order
tcga_sig_log2 <- tcga_sig_log2[
  intersect(rownames(tcga_sig_log2),
            gene_map$gene_id[gene_map$gene_name %in% sig_symbols]),
  , drop = FALSE
]
rownames(tcga_sig_log2) <- gene_map$gene_name[match(rownames(tcga_sig_log2), gene_map$gene_id)]

# Validation: log2 expression for matched signature genes
val_sig_log <- val_log[gene_match$val_row[gene_match$matched], , drop = FALSE]
rownames(val_sig_log) <- sig_symbols[gene_match$matched]

# Intersect genes available in both
common_genes <- intersect(rownames(tcga_sig_log2), rownames(val_sig_log))
message(sprintf("   Common signature genes with expression in both cohorts: %d",
                length(common_genes)))

if (length(common_genes) < 2) {
  stop("   < 2 common signature genes. Validation not possible.")
}

# Subset both matrices to common genes
tcga_sig_expr <- tcga_sig_log2[common_genes, , drop = FALSE]
val_sig_expr  <- val_sig_log[common_genes, , drop = FALSE]
common_coefs  <- coef_tbl$coefficient[match(common_genes, coef_tbl$gene_symbol)]

# Z-score each gene within TCGA
tcga_mean <- apply(tcga_sig_expr, 1, mean, na.rm = TRUE)
tcga_sd   <- apply(tcga_sig_expr, 1, sd, na.rm = TRUE)
tcga_sd[tcga_sd == 0] <- 1  # avoid division by zero for invariant genes

tcga_z <- (tcga_sig_expr - tcga_mean) / tcga_sd

# Z-score each gene within validation (using validation's own mean/SD)
val_mean <- apply(val_sig_expr, 1, mean, na.rm = TRUE)
val_sd   <- apply(val_sig_expr, 1, sd, na.rm = TRUE)
val_sd[val_sd == 0] <- 1

val_z <- (val_sig_expr - val_mean) / val_sd

# Save z-score reference parameters
z_params <- list(
  tcga_mean = tcga_mean,
  tcga_sd   = tcga_sd,
  val_mean  = val_mean,
  val_sd    = val_sd,
  common_genes = common_genes,
  coefficients = common_coefs,
  method   = "z-score per gene within each cohort independently",
  limitation = "Coefficients were estimated on log2(norm_counts+1) in TCGA, not on z-scored data"
)
saveRDS(z_params, file.path(OUT_DIR, "cross_platform_params.rds"))

# =============================================================================
# 6. Compute risk scores — primary (z-scored)
# =============================================================================
message(">>> Computing risk scores (z-score method) ...")

# Risk score = sum(gene_z * coefficient) for each patient
# TCGA
tcga_risk_z <- as.numeric(t(common_coefs) %*% tcga_z)
names(tcga_risk_z) <- colnames(tcga_z)

# Validation
val_risk_z <- as.numeric(t(common_coefs) %*% val_z)
names(val_risk_z) <- colnames(val_z)

message(sprintf("   TCGA risk scores: n=%d, median=%.4f, range=[%.4f, %.4f]",
                length(tcga_risk_z), median(tcga_risk_z),
                min(tcga_risk_z), max(tcga_risk_z)))
message(sprintf("   Validation risk scores: n=%d, median=%.4f, range=[%.4f, %.4f]",
                length(val_risk_z), median(val_risk_z),
                min(val_risk_z), max(val_risk_z)))

# =============================================================================
# 7. Compute risk scores — alternative (raw log2 scale, no z-score)
# =============================================================================
message(">>> Computing risk scores (raw log2 scale — alternative) ...")

# Directly apply coefficients to log2 expression (no z-score)
tcga_risk_raw <- as.numeric(t(common_coefs) %*% tcga_sig_expr)
names(tcga_risk_raw) <- colnames(tcga_sig_expr)

val_risk_raw <- as.numeric(t(common_coefs) %*% val_sig_expr)
names(val_risk_raw) <- colnames(val_sig_expr)

# =============================================================================
# 8. Load GSE107422 clinical data
# =============================================================================
message(">>> Loading GSE107422 clinical data ...")

eset_file <- file.path(GEO_RAW_DIR, "GSE107422_eset.rds")
if (!file.exists(eset_file)) {
  stop("GSE107422 eset not found: ", eset_file)
}
eset <- readRDS(eset_file)
pheno <- pData(eset)

message(sprintf("   PhenoData columns: %d", ncol(pheno)))

# -- Identify recurrence-related columns -------------------------------------
# GEO column names often end with ":ch1"
recur_cols <- grep("recur|relapse|progression|event|rfs|dfs|survival|follow.up|outcome",
                   names(pheno), ignore.case = TRUE, value = TRUE)
message("   Potential recurrence columns found:")
for (c in recur_cols) {
  vals <- unique(pheno[[c]])
  message(sprintf("     %-45s : %d unique values (e.g., %s)",
                  c, length(vals),
                  paste(head(na.omit(vals), 3), collapse = ", ")))
}

# -- Extract recurrence status and time --------------------------------------
# This is GEO-specific and column names are dataset-specific.
# Try common patterns:

# Pattern 1: "recurrence:ch1" or "recurrence status:ch1" — binary
recur_status_col <- grep("recurrence.*status|recurrence:ch1|recur$|recurrence_event",
                         names(pheno), ignore.case = TRUE, value = TRUE)

# Pattern 2: time to recurrence / recurrence-free survival time
recur_time_col <- grep("recurrence.*free|rfs|dfs|recurrence.*time|time.to.recur|recurrence.*days|recurrence.*month",
                       names(pheno), ignore.case = TRUE, value = TRUE)

# Pattern 3: follow-up time (if no specific RFS column)
fu_time_col <- grep("follow.up|last.follow|months.to|days.to",
                    names(pheno), ignore.case = TRUE, value = TRUE)

message("")
message("   Clinical data mapping:")
message(sprintf("     Recurrence status column(s) : %s",
                if (length(recur_status_col) > 0) paste(recur_status_col, collapse = ", ") else "NOT FOUND"))
message(sprintf("     Recurrence time column(s)   : %s",
                if (length(recur_time_col) > 0) paste(recur_time_col, collapse = ", ") else "NOT FOUND"))
message(sprintf("     Follow-up time column(s)    : %s",
                if (length(fu_time_col) > 0) paste(fu_time_col, collapse = ", ") else "NOT FOUND"))

# -- Build survival data frame -----------------------------------------------
samples_val <- colnames(val_log)

# Match samples to phenotype
# GEO sample IDs are typically in the first column or "geo_accession" column
if ("geo_accession" %in% names(pheno)) {
  pheno$sample_id <- pheno$geo_accession
} else if (length(intersect(rownames(pheno), samples_val)) > 0) {
  pheno$sample_id <- rownames(pheno)
} else {
  # Try matching by title or description
  pheno$sample_id <- make.names(pheno[[1]])  # best guess
}

# Build survival data
val_clin <- data.frame(
  sample_id = samples_val,
  stringsAsFactors = FALSE
)

# Merge phenotype data
pheno_sub <- pheno[, intersect(names(pheno),
                                c("sample_id", recur_status_col, recur_time_col,
                                  fu_time_col, "title:ch1", "characteristics_ch1")),
                    drop = FALSE]
val_clin <- merge(val_clin, pheno_sub, by = "sample_id", all.x = TRUE,
                  all = FALSE, sort = FALSE)

# -- Parse recurrence status -------------------------------------------------
if (length(recur_status_col) > 0) {
  # Try to extract binary status from the column
  raw_status <- val_clin[[recur_status_col[1]]]
  if (is.character(raw_status) || is.factor(raw_status)) {
    # Look for common yes/no, 1/0, recurrence/no_recurrence patterns
    status_numeric <- rep(NA, length(raw_status))
    status_numeric[grepl("^1$|yes|recurrence|recurrent|event|dead|progression",
                         raw_status, ignore.case = TRUE)] <- 1
    status_numeric[grepl("^0$|no|no_recurrence|no_recur|non.recur|censor|alive|no_event|no_recurrence|norecurrence|recurrence_free",
                         raw_status, ignore.case = TRUE)] <- 0
    val_clin$recur_event <- status_numeric
  } else if (is.numeric(raw_status)) {
    val_clin$recur_event <- as.integer(raw_status)
  }
}

# -- Parse recurrence time ---------------------------------------------------
if (length(recur_time_col) > 0) {
  val_clin$recur_time <- as.numeric(val_clin[[recur_time_col[1]]])
} else if (length(fu_time_col) > 0) {
  # Use follow-up time as proxy
  val_clin$recur_time <- as.numeric(val_clin[[fu_time_col[1]]])
  message("   WARNING: No explicit RFS column found. Using follow-up time as proxy.")
}

# Report clinical data quality
n_with_status <- sum(!is.na(val_clin$recur_event))
n_with_time   <- sum(!is.na(val_clin$recur_time))
message(sprintf("   Patients with recurrence status: %d / %d", n_with_status, nrow(val_clin)))
message(sprintf("   Patients with recurrence time  : %d / %d", n_with_time, nrow(val_clin)))

if (n_with_status == 0) {
  stop("Could not parse recurrence status from GSE107422 clinical annotation.")
}

logistic_only <- (n_with_time == 0)
if (logistic_only) {
  message("   NOTE: No recurrence time data available. Using logistic regression")
  message("   (recurrence status ~ risk group) as fallback instead of Cox PH.")
  message("   AUC will be reported instead of C-index.")
  # Create a dummy time column (all 1s) so KM/Cox code doesn't crash downstream
  val_clin$recur_time <- 1
  val_clin$recur_event[is.na(val_clin$recur_event)] <- 0
}

# =============================================================================
# 9. Stratify — primary: apply TCGA median cutoff
# =============================================================================
message(">>> Stratifying patients (primary: TCGA median cutoff) ...")

tcga_median_z <- median(tcga_risk_z, na.rm = TRUE)

val_clin$risk_score_z  <- val_risk_z
val_clin$risk_group_z  <- ifelse(val_risk_z > tcga_median_z, "High", "Low")
val_clin$risk_group_z  <- factor(val_clin$risk_group_z, levels = c("Low", "High"))
val_clin$risk_score_raw <- val_risk_raw
val_clin$risk_group_raw <- ifelse(val_risk_raw > median(tcga_risk_raw), "High", "Low")
val_clin$risk_group_raw <- factor(val_clin$risk_group_raw, levels = c("Low", "High"))

n_high <- sum(val_clin$risk_group_z == "High", na.rm = TRUE)
n_low  <- sum(val_clin$risk_group_z == "Low", na.rm = TRUE)
message(sprintf("   TCGA median cutoff applied: %.4f", tcga_median_z))
message(sprintf("   Validation: High-risk n=%d, Low-risk n=%d", n_high, n_low))

# If median split is highly imbalanced, report and also compute validation median
if (min(n_high, n_low) < nrow(val_clin) * 0.3) {
  message("   NOTE: Median split is imbalanced in validation cohort.")
  message("   Also computing validation-specific median cutoff for comparison.")
  val_median_z <- median(val_risk_z, na.rm = TRUE)
  val_clin$risk_group_val_z <- ifelse(val_risk_z > val_median_z, "High", "Low")
  val_clin$risk_group_val_z <- factor(val_clin$risk_group_val_z, levels = c("Low", "High"))
  message(sprintf("   Validation median cutoff: %.4f", val_median_z))
  message(sprintf("   Validation re-stratified: High n=%d, Low n=%d",
                  sum(val_clin$risk_group_val_z == "High"),
                  sum(val_clin$risk_group_val_z == "Low")))
}

# Save risk scores
risk_out <- val_clin %>%
  dplyr::select(sample_id, risk_score_z, risk_group_z, risk_score_raw, recur_time, recur_event)
write.csv(risk_out, file.path(OUT_DIR, "risk_scores_validation.csv"), row.names = FALSE)

# =============================================================================
# 10. Survival / logistic regression (validation, z-score method, TCGA median)
# =============================================================================
if (!logistic_only) {
  message(">>> Kaplan-Meier survival analysis (validation cohort) ...")
} else {
  message(">>> Logistic regression analysis (validation cohort, no time data) ...")
}

# Remove NAs
val_valid <- val_clin %>% dplyr::filter(!is.na(recur_time), !is.na(recur_event),
                                  !is.na(risk_group_z))
message(sprintf("   Patients in analysis: %d (%d events)",
                nrow(val_valid), sum(val_valid$recur_event)))

if (nrow(val_valid) < 20) {
  warning("Very few patients with complete data. Results may be unreliable.")
}

if (!logistic_only) {
  # KM with TCGA median cutoff
  km_fit_z <- survfit(Surv(recur_time, recur_event) ~ risk_group_z, data = val_valid)

  km_plot_z <- ggsurvplot(
    km_fit_z,
    data          = val_valid,
    pval          = TRUE,
    pval.method   = TRUE,
    conf.int      = TRUE,
    risk.table    = TRUE,
    risk.table.col = "strata",
    xlab          = "Time (days)",
    ylab          = "Recurrence-free survival probability",
    title         = "GSE107422 — CSC signature validation",
    subtitle      = sprintf("TCGA median cutoff (Low n=%d, High n=%d)",
                            sum(val_valid$risk_group_z == "Low"),
                            sum(val_valid$risk_group_z == "High")),
    palette       = c("#377EB8", "#E41A1C"),
    legend.title  = "Risk group",
    ggtheme       = theme_minimal()
  )

  png(file.path(FIG_DIR, "01_km_curve_median.png"), width = 9, height = 8, units = "in", res = 150)
  print(km_plot_z)
  dev.off()

  # Log-rank test
  lr_z <- survdiff(Surv(recur_time, recur_event) ~ risk_group_z, data = val_valid)
  lr_p_z <- 1 - pchisq(lr_z$chisq, df = 1)
  message(sprintf("   Log-rank test p-value (TCGA median cutoff): %.4f", lr_p_z))
} else {
  # Logistic regression fallback (no recurrence time data)
  logit_val <- glm(recur_event ~ risk_group_z, data = val_valid, family = binomial)
  logit_or <- exp(coef(logit_val)["risk_group_zHigh"])
  logit_ci <- exp(confint.default(logit_val, "risk_group_zHigh"))
  logit_p <- summary(logit_val)$coefficients["risk_group_zHigh", "Pr(>|z|)"]
  val_valid$pred_prob <- predict(logit_val, type = "response")
  roc_obj <- roc(val_valid$recur_event, val_valid$pred_prob)
  logit_auc <- as.numeric(auc(roc_obj))
  message(sprintf("   Logistic regression OR (High vs Low): %.4f (%.4f-%.4f)",
                  logit_or, logit_ci[1], logit_ci[2]))
  message(sprintf("   Logistic regression p-value: %.4f", logit_p))
  message(sprintf("   AUC: %.4f", logit_auc))
  lr_p_z <- logit_p
}

# =============================================================================
# 11. Cox proportional hazards + C-index / logistic regression
# =============================================================================
if (!logistic_only) {
  message(">>> Cox proportional hazards (validation cohort) ...")

  # Unadjusted
  cox_val <- coxph(Surv(recur_time, recur_event) ~ risk_group_z, data = val_valid)

  # Print summary
  sink(file.path(TAB_DIR, "Validation_cox_summary.txt"))
  cat("GSE107422 — Validation of TCGA CSC Signature\n")
  cat("===============================================\n\n")
  cat(sprintf("Cutoff: TCGA median (z-score method)\n"))
  cat(sprintf("Samples: %d, Events: %d\n\n", nrow(val_valid), sum(val_valid$recur_event)))

  cat("--- Unadjusted Cox ---\n")
  print(summary(cox_val))
  cat("\n")

  # C-index (Harrell's)
  cindex_val <- survConcordance(Surv(recur_time, recur_event) ~ predict(cox_val), data = val_valid)
  cat(sprintf("\nC-index (validation): %.4f (SE = %.4f)\n",
              cindex_val$concordance, sqrt(cindex_val$var)))

  # C-index for TCGA (using z-scored risk scores on TCGA data)
  # Build TCGA survival data and compute C-index
  tcga_surv <- data.frame(
    risk_score = tcga_risk_z[names(tcga_risk_z) %in% tcga_risk$sample_id],
    stringsAsFactors = FALSE
  )
  tcga_surv$os_time  <- tcga_risk$os_time[match(rownames(tcga_surv), tcga_risk$sample_id)]
  tcga_surv$os_event <- tcga_risk$os_event[match(rownames(tcga_surv), tcga_risk$sample_id)]
  tcga_surv <- tcga_surv[complete.cases(tcga_surv), ]

  cox_tcga_z <- coxph(Surv(os_time, os_event) ~ risk_score, data = tcga_surv)
  cindex_tcga <- survConcordance(Surv(os_time, os_event) ~ predict(cox_tcga_z), data = tcga_surv)

  cat("--- C-index Comparison ---\n")
  cat(sprintf("  TCGA (discovery)  : %.4f (SE = %.4f)\n",
              cindex_tcga$concordance, sqrt(cindex_tcga$var)))
  cat(sprintf("  GSE107422 (validation): %.4f (SE = %.4f)\n",
              cindex_val$concordance, sqrt(cindex_val$var)))
  cat(sprintf("  Delta (val - TCGA): %+.4f\n",
              cindex_val$concordance - cindex_tcga$concordance))

  # PH assumption test
  zph_val <- cox.zph(cox_val)
  cat("\n--- Proportional Hazards Assumption ---\n")
  print(zph_val)

  sink()

  message("   Cox results → ", file.path(TAB_DIR, "Validation_cox_summary.txt"))

  hr_val <- summary(cox_val)$conf.int[1]
  hr_p   <- summary(cox_val)$coefficients[1, "Pr(>|z|)"]
  hr_lower <- summary(cox_val)$conf.int[3]
  hr_upper <- summary(cox_val)$conf.int[4]
  cindex_val_val <- cindex_val$concordance
} else {
  message(">>> Logistic regression (validation cohort, no time data) ...")

  sink(file.path(TAB_DIR, "Validation_logistic_summary.txt"))
  cat("GSE107422 — Validation of TCGA CSC Signature\n")
  cat("(Logistic regression fallback — no recurrence time data)\n")
  cat("============================================================\n\n")
  cat(sprintf("Cutoff: TCGA median (z-score method)\n"))
  cat(sprintf("Samples: %d, Events: %d\n\n", nrow(val_valid), sum(val_valid$recur_event)))

  cat("--- Logistic Regression ---\n")
  print(summary(logit_val))
  cat("\n")
  cat(sprintf("Odds Ratio (High vs Low): %.4f (%.4f - %.4f)\n",
              logit_or, logit_ci[1], logit_ci[2]))
  cat(sprintf("p-value: %.4f\n", logit_p))
  cat(sprintf("AUC: %.4f\n", logit_auc))
  sink()

  message("   Logistic results → ", file.path(TAB_DIR, "Validation_logistic_summary.txt"))

  hr_val <- logit_or
  hr_p <- logit_p
  hr_lower <- logit_ci[1]
  hr_upper <- logit_ci[2]
  cindex_val_val <- logit_auc
}

# =============================================================================
# 12. Validation outcome classification (HARD RULE: honest reporting)
# =============================================================================
message(">>> Validation outcome classification ...")

# Criteria for outcome classification:
# SUCCESS:      p < 0.05 AND effect > 1 AND discrimination > 0.5
# PARTIAL:      p < 0.10 OR (effect > 1 AND discrimination > 0.5)
# FAIL:         effect <= 1 OR discrimination <= 0.5

effect_label <- if (logistic_only) "OR" else "HR"
cindex_label <- if (logistic_only) "AUC" else "C-index"

hr_consistent <- hr_val > 1

if (hr_p < 0.05 && hr_consistent && cindex_val_val > 0.5) {
  outcome <- "SUCCESS"
  outcome_detail <- sprintf(
    "Validation successful: %s=%.3f (%.3f-%.3f), p=%.4f, %s=%.4f. Significant discrimination in expected direction.",
    effect_label, hr_val, hr_lower, hr_upper, hr_p, cindex_label, cindex_val_val)
} else if ((hr_p < 0.10 || cindex_val_val > 0.55) && hr_consistent) {
  outcome <- "PARTIAL"
  outcome_detail <- sprintf(
    "Partial validation: %s=%.3f (%.3f-%.3f), p=%.4f, %s=%.4f. Direction consistent but %s.",
    effect_label, hr_val, hr_lower, hr_upper, hr_p, cindex_label, cindex_val_val,
    ifelse(hr_p >= 0.05, "not significant at alpha=0.05",
           sprintf("%s only moderately > 0.5", cindex_label)))
} else if (hr_consistent && cindex_val_val > 0.5) {
  outcome <- "PARTIAL"
  outcome_detail <- sprintf(
    "Partial validation: %s=%.3f (%.3f-%.3f), p=%.4f, %s=%.4f. Direction consistent with TCGA but not significant.",
    effect_label, hr_val, hr_lower, hr_upper, hr_p, cindex_label, cindex_val_val)
} else {
  outcome <- "FAIL"
  outcome_detail <- sprintf(
    "Validation failed: %s=%.3f (%.3f-%.3f), p=%.4f, %s=%.4f. %s",
    effect_label, hr_val, hr_lower, hr_upper, hr_p, cindex_label, cindex_val_val,
    ifelse(!hr_consistent,
           sprintf("%s direction opposite to TCGA expectation (High risk shows lower hazard).",
                   effect_label),
           sprintf("No evidence of discrimination (%s <= 0.5).", cindex_label)))
}

message(sprintf("   Validation outcome: %s", outcome))
message(sprintf("   Detail: %s", outcome_detail))

# Save outcome
sink(file.path(TAB_DIR, "Validation_outcome.txt"))
cat("Validation Outcome — TCGA CSC Signature in GSE107422\n")
cat("======================================================\n\n")
cat(sprintf("Outcome: %s\n\n", outcome))
cat(sprintf("Detail: %s\n\n", outcome_detail))
cat("--- Evidence ---\n")
cat(sprintf("  %s (High vs Low)  : %.4f (%.4f - %.4f)\n", effect_label, hr_val, hr_lower, hr_upper))
cat(sprintf("  p-value           : %.4f\n", hr_p))
cat(sprintf("  %s               : %.4f\n", cindex_label, cindex_val_val))
if (!logistic_only) {
  cat(sprintf("  TCGA C-index      : %.4f\n", cindex_tcga$concordance))
}
cat(sprintf("  N (validation)    : %d\n", nrow(val_valid)))
cat(sprintf("  Events            : %d\n", sum(val_valid$recur_event)))
cat("\n")
cat("Classification criteria:\n")
cat(sprintf("  SUCCESS : p < 0.05 AND %s > 1 AND %s > 0.5\n", effect_label, cindex_label))
cat(sprintf("  PARTIAL : p < 0.10 OR (%s > 1 AND %s > 0.5)\n", effect_label, cindex_label))
cat(sprintf("  FAIL    : %s <= 1 OR %s <= 0.5\n", effect_label, cindex_label))
cat("\n")
cat("Note: The signature was developed for overall survival in TCGA and\n")
cat("validated for recurrence-free survival in GSE107422. These are related\n")
cat("but distinct endpoints, which should be considered when interpreting\n")
cat("the validation result.\n")
if (logistic_only) {
  cat("\n")
  cat("Note: Logistic regression was used instead of Cox PH because\n")
  cat("GSE107422 clinical annotation contains recurrence status but no\n")
  cat("recurrence time / follow-up information.\n")
}
sink()

message("   Outcome → ", file.path(TAB_DIR, "Validation_outcome.txt"))

# =============================================================================
# 13. Sensitivity analysis: tertile-based stratification
# =============================================================================
if (!logistic_only) {
  message(">>> Sensitivity analysis: tertile-based stratification ...")
} else {
  message(">>> Sensitivity analysis: tertile-based stratification (logistic) ...")
}

# Tertile split within validation cohort
tertiles_val <- quantile(val_risk_z, probs = c(1/3, 2/3), na.rm = TRUE)
val_clin$risk_tertile <- cut(val_risk_z,
                              breaks = c(-Inf, tertiles_val[1], tertiles_val[2], Inf),
                              labels = c("Low", "Medium", "High"),
                              include.lowest = TRUE)

message(sprintf("   Validation tertile cutoffs: %.4f, %.4f",
                tertiles_val[1], tertiles_val[2]))
message(sprintf("   Low: %d, Medium: %d, High: %d",
                sum(val_clin$risk_tertile == "Low", na.rm = TRUE),
                sum(val_clin$risk_tertile == "Medium", na.rm = TRUE),
                sum(val_clin$risk_tertile == "High", na.rm = TRUE)))

# Filter for analysis
val_tert <- val_clin %>% dplyr::filter(!is.na(recur_time), !is.na(recur_event),
                                 !is.na(risk_tertile))

if (!logistic_only) {
  # KM by tertiles
  km_tert <- survfit(Surv(recur_time, recur_event) ~ risk_tertile, data = val_tert)

  km_plot_tert <- ggsurvplot(
    km_tert,
    data          = val_tert,
    pval          = TRUE,
    pval.method   = TRUE,
    conf.int      = FALSE,
    risk.table    = TRUE,
    risk.table.col = "strata",
    xlab          = "Time (days)",
    ylab          = "Recurrence-free survival probability",
    title         = "GSE107422 — Sensitivity: tertile split",
    subtitle      = sprintf("Low n=%d, Medium n=%d, High n=%d",
                            sum(val_tert$risk_tertile == "Low"),
                            sum(val_tert$risk_tertile == "Medium"),
                            sum(val_tert$risk_tertile == "High")),
    palette       = c("#377EB8", "#4DAF4A", "#E41A1C"),
    legend.title  = "Risk group",
    ggtheme       = theme_minimal()
  )

  png(file.path(FIG_DIR, "02_km_curve_tertiles.png"), width = 9, height = 8, units = "in", res = 150)
  print(km_plot_tert)
  dev.off()

  # Cox with tertiles (treat as ordered or categorical)
  cox_tert <- coxph(Surv(recur_time, recur_event) ~ risk_tertile, data = val_tert)

  # Pairwise: High vs Low
  cox_high_low <- coxph(Surv(recur_time, recur_event) ~ I(risk_tertile == "High"),
                         data = val_tert)

  # Log-rank trend test (score test for trend ~ numeric tertile)
  val_tert$tertile_num <- as.numeric(val_tert$risk_tertile)
  cox_trend <- coxph(Surv(recur_time, recur_event) ~ tertile_num, data = val_tert)
} else {
  # Logistic regression by tertiles
  logit_tert <- glm(recur_event ~ risk_tertile, data = val_tert, family = binomial)
  logit_tert_high_vs_low <- glm(recur_event ~ I(risk_tertile == "High"), data = val_tert, family = binomial)
  val_tert$tertile_num <- as.numeric(val_tert$risk_tertile)
  logit_tert_trend <- glm(recur_event ~ tertile_num, data = val_tert, family = binomial)
}

# Save sensitivity results
sink(file.path(TAB_DIR, "Validation_tertile_sensitivity.txt"))
cat("Sensitivity Analysis — Tertile-Based Stratification\n")
cat("====================================================\n\n")
cat(sprintf("Validation dataset: GSE107422\n"))
cat(sprintf("Tertile cutoffs: %.4f, %.4f\n\n", tertiles_val[1], tertiles_val[2]))

cat("--- Tertile distribution ---\n")
print(table(val_tert$risk_tertile, useNA = "ifany"))
cat("\n")

if (!logistic_only) {
  cat("--- Cox by tertile ---\n")
  print(summary(cox_tert))
  cat("\n")

  cat("--- Cox: High vs Low (pairwise) ---\n")
  print(summary(cox_high_low))
  cat("\n")

  cat("--- Trend test (tertile as numeric) ---\n")
  print(summary(cox_trend))
  cat("\n")

  cat("--- C-index (tertile model) ---\n")
  cindex_tert <- survConcordance(Surv(recur_time, recur_event) ~ predict(cox_tert), data = val_tert)
  cat(sprintf("  C-index (tertile): %.4f (SE = %.4f)\n",
              cindex_tert$concordance, sqrt(cindex_tert$var)))
  cat("\n")

  cat("--- Comparison: median vs tertile ---\n")
  cat(sprintf("  Median split C-index : %.4f\n", cindex_val$concordance))
  cat(sprintf("  Tertile split C-index: %.4f\n", cindex_tert$concordance))
} else {
  cat("--- Logistic regression by tertile ---\n")
  print(summary(logit_tert))
  cat("\n")

  cat("--- Logistic regression: High vs Low (pairwise) ---\n")
  print(summary(logit_tert_high_vs_low))
  cat("\n")

  cat("--- Trend test (tertile as numeric) ---\n")
  print(summary(logit_tert_trend))
  cat("\n")

  # AUC for tertile model
  val_tert$pred_prob_tert <- predict(logit_tert, type = "response")
  roc_tert <- roc(val_tert$recur_event, val_tert$pred_prob_tert)
  auc_tert <- as.numeric(auc(roc_tert))
  cat(sprintf("  AUC (tertile model): %.4f\n", auc_tert))
  cat("\n")

  cat("--- Comparison: median vs tertile ---\n")
  cat(sprintf("  Median split AUC : %.4f\n", cindex_val_val))
  cat(sprintf("  Tertile split AUC: %.4f\n", auc_tert))
}
sink()

message("   Sensitivity → ", file.path(TAB_DIR, "Validation_tertile_sensitivity.txt"))

# =============================================================================
# 14. Additional figures
# =============================================================================
message(">>> Additional figures ...")

# Risk score distribution
p_dist <- ggplot(val_clin, aes(x = risk_score_z)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, colour = "black") +
  geom_vline(xintercept = tcga_median_z, linetype = "dashed", colour = "red", size = 1) +
  geom_vline(xintercept = tertiles_val, linetype = "dotted", colour = "darkgreen", size = 0.8) +
  annotate("text", x = tcga_median_z, y = 0, label = "TCGA median",
           hjust = -0.1, vjust = -0.5, colour = "red", size = 3) +
  labs(title = "GSE107422 — Risk score distribution",
       subtitle = "Red dashed = TCGA median; Green dotted = validation tertiles",
       x = "Risk score (z-score method)", y = "Count") +
  theme_minimal()

ggsave(file.path(FIG_DIR, "05_risk_score_distribution.png"), p_dist, width = 8, height = 5)

# Forest plot
hr_val_df <- data.frame(
  variable = "Risk group (High vs Low)",
  HR = hr_val,
  HR_lower = hr_lower,
  HR_upper = hr_upper,
  p = hr_p,
  stringsAsFactors = FALSE
)

cindex_label_plot <- if (logistic_only) "AUC" else "C-index"
hr_label <- if (logistic_only) "Odds Ratio" else "Hazard Ratio"

p_forest <- ggplot(hr_val_df, aes(x = HR, y = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2, size = 1) +
  geom_point(size = 4, colour = "steelblue") +
  scale_x_log10() +
  labs(title = "GSE107422 — Validation regression",
       subtitle = sprintf("p = %.4f, %s = %.4f", hr_p, cindex_label_plot, cindex_val_val),
       x = sprintf("%s (95%% CI)", hr_label), y = "") +
  theme_minimal()

ggsave(file.path(FIG_DIR, "03_forest_plot.png"), p_forest, width = 6, height = 3)

if (!logistic_only) {
  # PH assumption plot
  png(file.path(FIG_DIR, "04_cox_zph_plot.png"), width = 8, height = 6, units = "in", res = 150)
  par(mfrow = c(1, 1))
  plot(zph_val)
  dev.off()
}

# =============================================================================
# 15. Final summary
# =============================================================================
message("")
message("==============================================")
message("   VALIDATION COMPLETE                       ")
message("==============================================")
message("")
message(sprintf("  Signature genes          : %d", length(common_genes)))
message(sprintf("  Validation cohort        : GSE107422 (Korean CRC)"))
message(sprintf("  Normalisation            : %s", data_type))
message(sprintf("  Cross-platform method    : z-score per gene per cohort"))
message(sprintf("  Risk score cutoff        : TCGA median (%.4f)", tcga_median_z))
message(sprintf("  Validation outcome       : %s", outcome))
message(sprintf("  p-value                  : %.4f", hr_p))

effect_label_summary <- if (logistic_only) "OR" else "HR"
cindex_label_summary <- if (logistic_only) "AUC" else "C-index"
message(sprintf("  %s (High vs Low)     : %.4f (%.4f - %.4f)", effect_label_summary, hr_val, hr_lower, hr_upper))
if (logistic_only) {
  message(sprintf("  AUC (validation)         : %.4f", cindex_val_val))
} else {
  message(sprintf("  C-index (validation)     : %.4f", cindex_val_val))
  message(sprintf("  C-index (TCGA discovery) : %.4f", cindex_tcga$concordance))
}

message("")
message(">>> Outputs:")
message("     ", file.path(OUT_DIR, "risk_scores_validation.csv"))
message("     ", file.path(OUT_DIR, "cross_platform_params.rds"))
if (logistic_only) {
  message("     ", file.path(TAB_DIR, "Validation_logistic_summary.txt"))
} else {
  message("     ", file.path(TAB_DIR, "Validation_cox_summary.txt"))
}
message("     ", file.path(TAB_DIR, "Validation_outcome.txt"))
message("     ", file.path(TAB_DIR, "Validation_tertile_sensitivity.txt"))
message("     ", FIG_DIR, "/ (5 figures)")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info_07.txt"))
