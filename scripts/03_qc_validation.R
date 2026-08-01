#!/usr/bin/env Rscript
# =============================================================================
# 03_qc_validation.R — Quality control for GEO validation cohorts
#
# Purpose:
#   Apply identical QC rules to both Asian validation cohorts (GSE107422 and
#   GSE92921, GSE71187, GSE64857). Every exclusion is logged with stated rule and written to
#   Table S2.
#
# Inputs:
#   data/raw/geo/GSE107422/GSE107422_eset.rds
#   (GSE220148 removed — expression-only reference cohort excluded from study)
#
# Outputs:
#   data/processed/geo_qc/
#   ├── GSE107422_dds.rds / _norm.rds
#   (GSE220148_dds.rds / _norm.rds removed)
#   ├── geo_excluded_samples.csv          — Table S2
#   └── geo_qc_summary.txt
#
#   results/figures/qc_geo/
#   ├── GSE107422_* (6 figures)
#   (GSE220148 QC figures removed)
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

# -- Directories -------------------------------------------------------------
GEO_DIR    <- file.path("data", "raw", "geo")
OUT_DIR    <- file.path("data", "processed", "geo_qc")
FIG_DIR    <- file.path("results", "figures", "qc_geo")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# -- Shared QC constants (identical to discovery script) ---------------------
MIN_COUNTS      <- 10
MIN_SAMPLE_PROP <- 0.1
PCA_OUTLIER_SD  <- 4

# =============================================================================
# Helper: initialise and return exclusion table
# =============================================================================
init_exclusions <- function() {
  data.frame(
    sample_id      = character(),
    cohort         = character(),
    exclusion_step = character(),
    reason         = character(),
    rule           = character(),
    stringsAsFactors = FALSE
  )
}

add_exclusion <- function(excl, sid, cohort, step, reason, rule) {
  rbind(excl, data.frame(
    sample_id = sid, cohort = cohort,
    exclusion_step = step, reason = reason, rule = rule,
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# Helper: PCA outlier detection + plot
# =============================================================================
qc_pca <- function(expr_mat, sample_meta, title, fig_path, excl, cohort_name,
                   pca_outlier_sd = PCA_OUTLIER_SD) {

  if (ncol(expr_mat) < 4) {
    message("     Too few samples for PCA, skipping.")
    return(list(plot = NULL, excl = excl, n_outliers = 0))
  }

  # Remove rows and columns with NA/Inf values
  pca_mat <- expr_mat
  # Remove rows (genes) with any NA
  pca_mat <- pca_mat[complete.cases(pca_mat), , drop = FALSE]
  # Remove columns (samples) with any NA
  pca_mat <- pca_mat[, colSums(is.na(pca_mat)) == 0, drop = FALSE]
  # Replace any remaining Inf with max finite
  pca_mat[is.infinite(pca_mat)] <- NA
  pca_mat <- pca_mat[complete.cases(pca_mat), , drop = FALSE]

  if (ncol(pca_mat) < 4 || nrow(pca_mat) < 10) {
    message("     Too few features after NA removal for PCA, skipping.")
    return(list(plot = NULL, excl = excl, n_outliers = 0))
  }
  message(sprintf("     PCA on %d genes x %d samples", nrow(pca_mat), ncol(pca_mat)))

  pca <- prcomp(t(pca_mat), center = TRUE, scale. = FALSE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2], label = rownames(pca$x)
  )
  # Merge available metadata
  if (!is.null(sample_meta) && nrow(sample_meta) > 0) {
    df <- cbind(df, sample_meta[rownames(df), , drop = FALSE])
  }

  center_x <- mean(df$PC1, na.rm = TRUE)
  center_y <- mean(df$PC2, na.rm = TRUE)
  sd_x <- sd(df$PC1, na.rm = TRUE)
  sd_y <- sd(df$PC2, na.rm = TRUE)

  df$outlier <- abs(df$PC1 - center_x) > pca_outlier_sd * sd_x |
                abs(df$PC2 - center_y) > pca_outlier_sd * sd_y

  p <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(size = 2.5, alpha = 0.8) +
    ggrepel::geom_text_repel(data = subset(df, outlier),
                             aes(label = label), size = 2.5, max.overlaps = 20) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
    theme_minimal()

  if (!is.null(fig_path)) ggsave(fig_path, p, width = 8, height = 6)

  # Flag outliers
  outlier_ids <- rownames(df[df$outlier, ])
  for (s in outlier_ids) {
    excl <- add_exclusion(excl, s, cohort_name, "pca_outlier",
      sprintf("> %d SD from centroid on PC1 or PC2", pca_outlier_sd),
      sprintf("|PC| <= %d SD from mean", pca_outlier_sd))
  }

  list(plot = p, excl = excl, n_outliers = length(outlier_ids))
}

# =============================================================================
# Helper: sample correlation heatmap
# =============================================================================
qc_heatmap <- function(expr_mat, annotation, fig_path, title, max_samples = 100) {
  if (ncol(expr_mat) < 3) {
    message("     Too few samples for heatmap, skipping.")
    return(invisible(NULL))
  }

  set.seed(42)
  if (ncol(expr_mat) > max_samples) {
    idx <- sample(ncol(expr_mat), max_samples)
    hm_mat <- expr_mat[, idx]
    hm_annot <- if (!is.null(annotation)) annotation[idx, , drop = FALSE] else NULL
  } else {
    hm_mat <- expr_mat
    hm_annot <- annotation
  }

  # Remove rows/cols with NA/Inf in expression data
  hm_clean <- hm_mat
  hm_clean[is.na(hm_clean) | is.infinite(hm_clean)] <- 0

  cor_mat <- cor(hm_clean)
  # Replace any NA correlations (zero-variance samples) with 0
  cor_mat[is.na(cor_mat)] <- 0

  if (!is.null(fig_path)) {
    png(fig_path, width = 10, height = 9, units = "in", res = 150)
    pheatmap(cor_mat,
             annotation_col = hm_annot,
             show_colnames = FALSE, show_rownames = FALSE,
             main = title,
             color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean")
    dev.off()
  }
}

# =============================================================================
# Helper: library size plot
# =============================================================================
qc_libsize <- function(lib_sizes, sample_meta, fig_path, title) {
  df <- data.frame(
    sample = names(lib_sizes),
    lib_size = lib_sizes
  )
  if (!is.null(sample_meta) && nrow(sample_meta) > 0) {
    common <- intersect(names(lib_sizes), rownames(sample_meta))
    df <- cbind(df, sample_meta[common, , drop = FALSE])
  }

  p <- ggplot(df, aes(x = reorder(sample, lib_size), y = lib_size)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.9) +
    labs(title = title,
         x = "Sample", y = "Library size (total counts / expression units)") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank())
  if (!is.null(fig_path)) ggsave(fig_path, p, width = 12, height = 4.5)
  p
}

# =============================================================================
# Main QC function for one GEO dataset
# =============================================================================
qc_one_geo_dataset <- function(gse_id, excl) {

  message("\n==============================================")
  message("   QC for ", gse_id)
  message("==============================================")

  series_dir <- file.path(GEO_DIR, gse_id)
  eset_file <- file.path(series_dir, paste0(gse_id, "_eset.rds"))
  mat_file <- file.path(series_dir, paste0(gse_id, "_expression_matrix.rds"))

  if (!file.exists(eset_file) && !file.exists(mat_file)) {
    warning("   No data found for ", gse_id, ". Run 02_download_geo_validation.R first.")
    return(excl)
  }

  # Load metadata (pheno) from ExpressionSet
  pheno <- NULL
  if (file.exists(eset_file)) {
    eset <- readRDS(eset_file)
    pheno <- pData(eset)
    # Try expression data from ExpressionSet
    expr_mat <- tryCatch(exprs(eset), error = function(e) NULL)
  } else {
    expr_mat <- NULL
  }

  # Prefer pre-parsed expression matrix (from 02b_parse_geo_supplementary.R)
  if (file.exists(mat_file)) {
    expr_mat <- readRDS(mat_file)
    message("   Loaded pre-parsed expression matrix: ", nrow(expr_mat), " x ", ncol(expr_mat))
  }

  if (is.null(expr_mat) || ncol(expr_mat) < 3) {
    warning("   No usable expression matrix for ", gse_id, ". Reporting metadata only.")
    return(excl)
  }

  # -- 2. Determine data type (counts vs TPM/FPKM) --------------------------
  is_counts <- all(expr_mat == floor(expr_mat), na.rm = TRUE) && max(expr_mat, na.rm = TRUE) > 100
  message(sprintf("   Data type: %s", ifelse(is_counts, "integer counts", "continuous (TPM/FPKM)")))

  # -- 3. Library size -------------------------------------------------------
  message(">>> Library size distribution ...")
  if (is_counts) {
    lib_sizes <- colSums(expr_mat, na.rm = TRUE)
  } else {
    lib_sizes <- colSums(expr_mat, na.rm = TRUE)  # still informative for TPM
  }

  qc_libsize(lib_sizes, pheno, file.path(FIG_DIR, paste0(gse_id, "_01_library_size.png")),
             paste0(gse_id, " — Library size per sample"))

  message(sprintf("   Median: %s, Min: %s, Max: %s",
                  format(median(lib_sizes, na.rm = TRUE), big.mark = ","),
                  format(min(lib_sizes, na.rm = TRUE), big.mark = ","),
                  format(max(lib_sizes, na.rm = TRUE), big.mark = ",")))

  # -- 4. Gene filtering (only for raw counts) ------------------------------
  if (is_counts) {
    message(">>> Gene-level filtering ...")
    min_samples <- ceiling(ncol(expr_mat) * MIN_SAMPLE_PROP)
    keep <- rowSums(expr_mat >= MIN_COUNTS, na.rm = TRUE) >= min_samples
    message(sprintf("   Kept %d / %d genes (>=%d counts in >=%.0f%% samples)",
                    sum(keep), nrow(expr_mat), MIN_COUNTS, MIN_SAMPLE_PROP * 100))
    expr_filt <- expr_mat[keep, ]
  } else {
    message(">>> Continuous data (TPM/FPKM) — keeping all features with non-zero median > 0")
    # For TPM, keep genes with median > 0 (expressed in at least half the samples)
    keep <- apply(expr_mat, 1, function(x) {
      med <- median(x, na.rm = TRUE)
      !is.na(med) && med > 0 && is.finite(med)
    })
    expr_filt <- expr_mat[keep, ]
    message(sprintf("   Kept %d / %d genes (median expression > 0)", sum(keep), nrow(expr_mat)))
  }

  # -- 5. Transformation for PCA/heatmap (log2, or rlog for counts) ---------
  message(">>> Transformation for visualisation ...")
  if (is_counts && ncol(expr_filt) >= 4) {
    # Try DESeq2 rlog
    tryCatch({
      dds_tmp <- DESeqDataSetFromMatrix(
        countData = expr_filt,
        colData = data.frame(row.names = colnames(expr_filt), dummy = 1),
        design = ~ 1
      )
      dds_tmp <- estimateSizeFactors(dds_tmp)
      rld_tmp <- rlog(dds_tmp, blind = TRUE)
      viz_mat <- assay(rld_tmp)
      message("     Used DESeq2 rlog transformation")
    }, error = function(e) {
      # Fallback to log2-CPM via edgeR or simple log2
      message("     rlog failed, using log2-CPM: ", e$message)
      if (is_counts) {
        viz_mat <<- cpm(expr_filt, log = TRUE, prior.count = 1)
      } else {
        viz_mat <<- log2(expr_filt + 1)
      }
    })
  } else {
    viz_mat <- log2(expr_filt + 1)
    message("     Used log2(x + 1) transformation")
  }

  # Clean NA/Inf values for downstream analyses
  viz_mat[is.na(viz_mat) | is.infinite(viz_mat)] <- 0

  # -- 6. PCA ----------------------------------------------------------------
  message(">>> PCA for outlier detection ...")
  pca_res <- qc_pca(viz_mat, pheno,
                    title = paste0(gse_id, " — PCA (post-filter)"),
                    fig_path = file.path(FIG_DIR, paste0(gse_id, "_02_pca.png")),
                    excl = excl, cohort_name = gse_id)
  excl <- pca_res$excl

  # -- 7. Sample correlation heatmap ----------------------------------------
  message(">>> Sample correlation heatmap ...")
  # Build annotation from pheno (pick first few useful columns)
  useful_cols <- setdiff(names(pheno), c("title", "geo_accession", "status", "submission_date",
                                          "last_update_date", "type", "channel_count"))
  annot <- pheno[, intersect(useful_cols, names(pheno)), drop = FALSE]
  annot <- annot[, sapply(annot, function(x) length(unique(x))) < nrow(annot) * 0.9, drop = FALSE]
  annot <- annot[, sapply(annot, function(x) length(unique(x))) > 1, drop = FALSE]

  if (ncol(annot) > 5) annot <- annot[, 1:5]  # cap at 5 annotation columns

  qc_heatmap(viz_mat, annot,
             file.path(FIG_DIR, paste0(gse_id, "_03_heatmap.png")),
             paste0(gse_id, " — Sample correlation"))
  
  # -- 8. Missing clinical data ---------------------------------------------
  message(">>> Missing clinical data ...")
  clin_missing <- data.frame(
    variable = names(pheno),
    n_missing = sapply(pheno, function(x) sum(is.na(x) | x == "")),
    pct_missing = round(100 * sapply(pheno, function(x) sum(is.na(x) | x == "")) / nrow(pheno), 1)
  ) %>% arrange(desc(n_missing)) %>% head(20)

  for (i in seq_len(min(10, nrow(clin_missing)))) {
    r <- clin_missing[i, ]
    message(sprintf("     %-35s : %d / %d missing (%.0f%%)",
                    r$variable, r$n_missing, nrow(pheno), r$pct_missing))
  }

  # -- 9. Save processed data -----------------------------------------------
  message(">>> Saving processed data ...")
  dataset_dir <- file.path(OUT_DIR, gse_id)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

  if (is_counts) {
    saveRDS(expr_filt, file.path(dataset_dir, "counts_filtered.rds"))
    saveRDS(viz_mat, file.path(dataset_dir, "log_normalized.rds"))
  } else {
    saveRDS(expr_filt, file.path(dataset_dir, "tpm_filtered.rds"))
    saveRDS(viz_mat, file.path(dataset_dir, "log2_tpm.rds"))
  }

  message(">>> QC for ", gse_id, " complete.")

  excl
}

# =============================================================================
# Run QC on both GEO datasets
# =============================================================================
message(">>> Starting QC for validation cohorts ...")

all_excl <- init_exclusions()

all_excl <- qc_one_geo_dataset("GSE107422", all_excl)
# GSE220148 removed — expression-only reference cohort excluded from study


# -- Write exclusion table (Table S2) ----------------------------------------
write.csv(all_excl, file.path(OUT_DIR, "geo_excluded_samples.csv"), row.names = FALSE)

# -- Write summary ------------------------------------------------------------
sink(file.path(OUT_DIR, "geo_qc_summary.txt"))
cat("GEO Validation Cohorts — QC Summary\n")
cat("=====================================\n\n")
for (gse in c("GSE107422")) {
  cat(sprintf("--- %s ---\n", gse))
  eset_file <- file.path(GEO_DIR, gse, paste0(gse, "_eset.rds"))
  if (file.exists(eset_file)) {
    eset <- readRDS(eset_file)
    cat(sprintf("  Samples loaded : %d\n", ncol(eset)))
  }
  excl_sub <- all_excl[all_excl$cohort == gse, ]
  if (nrow(excl_sub) > 0) {
    cat(sprintf("  Excluded       : %d\n", nrow(excl_sub)))
    for (i in seq_len(nrow(excl_sub))) {
      cat(sprintf("    - %s: %s\n", excl_sub$sample_id[i], excl_sub$reason[i]))
    }
  } else {
    cat("  Excluded       : 0\n")
  }
  cat("\n")
}
cat("— End of QC summary —\n")
sink()

message("")
message(">>> QC complete for both validation cohorts.")
message("     Exclusions: ", file.path(OUT_DIR, "geo_excluded_samples.csv"))
message("     Summary   : ", file.path(OUT_DIR, "geo_qc_summary.txt"))
message("     Figures   : ", FIG_DIR, "/")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))
