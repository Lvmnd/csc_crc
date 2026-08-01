#!/usr/bin/env Rscript
# =============================================================================
# 03_qc_discovery.R — Quality control & filtering for TCGA discovery cohort
#
# Purpose:
#   Apply consistent, documented QC filters to TCGA-COAD + TCGA-READ raw
#   STAR-Counts before downstream analysis. Every exclusion is logged with
#   a stated rule and written to Table S2.
#
# Inputs:
#   data/raw/tcga/tcga_expression_SE.rds      — SummarizedExperiment
#   data/raw/tcga/clinical_COAD.rds
#   data/raw/tcga/clinical_READ.rds
#
# Outputs:
#   data/processed/tcga_qc/
#   ├── tcga_dds.rds                          — DESeqDataSet (filtered, pre-normalization)
#   ├── tcga_vsd.rds                          — DESeqDataSet (after vst)
#   ├── tcga_excluded_samples.csv             — Table S2 skeleton
#   └── tcga_qc_summary.txt                   — numerics for manuscript
#
#   results/figures/qc_tcga/
#   ├── 01_library_size_boxplot.png
#   ├── 02_pca_pre_filter.png
#   ├── 03_pca_post_filter.png
#   ├── 04_sample_correlation_heatmap.png
#   ├── 05_gene_filtering_diagnostic.png
#   └── 06_missing_clinical_data.png
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(matrixStats)
})

# Helper: extract TCGA sample type from barcode (chars 14-15)
# "01" = Primary Solid Tumor, "11" = Solid Tissue Normal
tcga_sample_type <- function(barcodes) {
  code <- substr(barcodes, 14, 15)
  ifelse(code == "01", "Tumor", ifelse(code == "11", "Normal", "Other"))
}
tcga_is_tumor <- function(barcodes) substr(barcodes, 14, 15) == "01"
tcga_is_normal <- function(barcodes) substr(barcodes, 14, 15) == "11"

# -- Directories -------------------------------------------------------------
RAW_DIR    <- file.path("data", "raw", "tcga")
OUT_DIR    <- file.path("data", "processed", "tcga_qc")
FIG_DIR    <- file.path("results", "figures", "qc_tcga")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# -- Constants (documented exclusion rules — see logs/decisions.md) ----------
MIN_COUNTS      <- 10    # minimum raw count
MIN_SAMPLE_PROP <- 0.1   # gene must have >= MIN_COUNTS in >= this fraction of samples
                         # Equivalent to: keep gene if counts >= 10 in >= 10% of samples
                         # (Removes genes that are essentially unexpressed)
MAX_NA_CLINICAL <- 0.5   # exclude clinical variable if >50% missing
PCA_OUTLIER_SD  <- 4     # flag samples >4 SD from mean on PC1 or PC2

# =============================================================================
# 1. Load data
# =============================================================================
message(">>> Loading TCGA expression and clinical data ...")

se <- readRDS(file.path(RAW_DIR, "tcga_expression_SE.rds"))
clin_COAD <- readRDS(file.path(RAW_DIR, "clinical_COAD.rds"))
clin_READ <- readRDS(file.path(RAW_DIR, "clinical_READ.rds"))

# Raw count matrix
counts_raw <- assay(se, "unstranded")
message("   Raw dimensions: ", nrow(counts_raw), " genes x ", ncol(counts_raw), " samples")

# Column data
col_data <- as.data.frame(colData(se))
col_data$sample_barcode <- colnames(se)
col_data$sample_type <- tcga_sample_type(colnames(se))

# =============================================================================
# 2. Initialise exclusion tracking (Table S2)
# =============================================================================
exclusions <- data.frame(
  sample_id     = character(),
  cohort        = character(),
  exclusion_step = character(),
  reason         = character(),
  rule           = character(),
  stringsAsFactors = FALSE
)

add_exclusion <- function(sid, cohort, step, reason, rule) {
  exclusions <<- rbind(exclusions, data.frame(
    sample_id = sid, cohort = cohort,
    exclusion_step = step, reason = reason, rule = rule,
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# 3. Library size distribution
# =============================================================================
message(">>> Library size distribution ...")

lib_sizes <- colSums(counts_raw)
lib_df <- data.frame(
  sample = colnames(counts_raw),
  lib_size = lib_sizes,
  project = col_data$project_id,
  type = col_data$sample_type
)

p <- ggplot(lib_df, aes(x = reorder(sample, lib_size), y = lib_size, fill = project)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "TCGA — Library size per sample",
       subtitle = paste0("Total samples = ", ncol(se), ", Total genes = ", nrow(counts_raw)),
       x = "Sample", y = "Total counts") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())
ggsave(file.path(FIG_DIR, "01_library_size_boxplot.png"), p, width = 14, height = 5)

message(sprintf("   Library sizes — median: %s, min: %s, max: %s",
                format(median(lib_sizes), big.mark = ","),
                format(min(lib_sizes), big.mark = ","),
                format(max(lib_sizes), big.mark = ",")))

# Flag very low library size (< 1M counts)
low_lib_cutoff <- 1e6
low_lib_samples <- names(lib_sizes[lib_sizes < low_lib_cutoff])
for (s in low_lib_samples) {
  add_exclusion(s, col_data[s, "project_id"], "library_size",
                sprintf("Library size = %d (< %d)", lib_sizes[s], low_lib_cutoff),
                "library_size >= 1e6 total counts")
}
message(sprintf("   Samples with library size < %s: %d",
                format(low_lib_cutoff, big.mark = ","), length(low_lib_samples)))

# =============================================================================
# 4. Gene filtering — remove low-count genes
# =============================================================================
message(">>> Gene-level filtering ...")

# Define filter: gene must have >= MIN_COUNTS in >= MIN_SAMPLE_PROP of samples
min_samples <- ceiling(ncol(counts_raw) * MIN_SAMPLE_PROP)
keep_genes <- rowSums(counts_raw >= MIN_COUNTS) >= min_samples
n_removed <- sum(!keep_genes)
n_kept <- sum(keep_genes)

message(sprintf("   Genes with >= %d counts in >= %d samples (%.0f%%): %d / %d (%.1f%% kept)",
                MIN_COUNTS, min_samples, MIN_SAMPLE_PROP * 100,
                n_kept, nrow(counts_raw), 100 * n_kept / nrow(counts_raw)))

# Diagnostic plot (use matrixStats for speed)
gene_diag <- data.frame(
  max_count = rowMaxs(counts_raw),
  n_nonzero = rowSums(counts_raw > 0),
  kept = keep_genes
)
p_gf <- ggplot(gene_diag, aes(x = log10(max_count + 1), y = n_nonzero, colour = kept)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = min_samples, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = log10(MIN_COUNTS + 1), linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "steelblue")) +
  labs(title = "TCGA — Gene filtering diagnostic",
       subtitle = sprintf("Threshold: >= %d counts in >= %d samples. Removed: %d genes",
                          MIN_COUNTS, min_samples, n_removed),
       x = "log10(max count + 1)", y = "Number of samples with count > 0") +
  theme_minimal()
ggsave(file.path(FIG_DIR, "05_gene_filtering_diagnostic.png"), p_gf, width = 8, height = 6)

# Apply filter
counts_filt <- counts_raw[keep_genes, ]
se_filt <- se[keep_genes, ]

# =============================================================================
# 5. Build DESeqDataSet and normalise
# =============================================================================
message(">>> Building DESeqDataSet ...")

# Coldata for DESeq2
cd <- col_data[, c("sample_barcode", "project_id", "sample_type"), drop = FALSE]
cd$project_id <- factor(cd$project_id)
cd$sample_type <- factor(cd$sample_type)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = cd,
  design = ~ sample_type + project_id   # adjust for project (COAD vs READ)
)

# Median-of-ratios normalisation (DESeq2 default)
message(">>> Estimating size factors (median-of-ratios) ...")
dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)
message(sprintf("   Size factor range: %.3f – %.3f", min(sizeFactors(dds)), max(sizeFactors(dds))))

# vst transformation for PCA / clustering (blind — much faster than rlog)
message(">>> vst transformation for visualisation ...")
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

# =============================================================================
# 6. PCA — pre-filter outlier detection
# =============================================================================
message(">>> PCA for outlier detection ...")

pca_plot <- function(mat, cd, pc_x = 1, pc_y = 2, title = "") {
  pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df <- data.frame(
    PC1 = pca$x[, pc_x],
    PC2 = pca$x[, pc_y],
    project = cd$project_id,
    type = cd$sample_type
  )
  df$label <- rownames(df)

  # Outlier flags based on distance
  center_x <- mean(df$PC1)
  center_y <- mean(df$PC2)
  sd_x <- sd(df$PC1)
  sd_y <- sd(df$PC2)

  # Detect samples > PCA_OUTLIER_SD from mean on either PC
  df$outlier <- abs(df$PC1 - center_x) > PCA_OUTLIER_SD * sd_x |
                abs(df$PC2 - center_y) > PCA_OUTLIER_SD * sd_y

  p <- ggplot(df, aes(x = PC1, y = PC2, colour = type, shape = project)) +
    geom_point(size = 2.5, alpha = 0.8) +
    ggrepel::geom_text_repel(data = subset(df, outlier),
                             aes(label = label), size = 2.5, max.overlaps = 20) +
    stat_ellipse(aes(group = type), level = 0.95, linetype = "dashed", alpha = 0.4) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = title,
         x = sprintf("PC%d (%.1f%%)", pc_x, var_exp[pc_x]),
         y = sprintf("PC%d (%.1f%%)", pc_y, var_exp[pc_y])) +
    theme_minimal()

  list(plot = p, outlier_df = df)
}

# PCA on vst-transformed counts
pca_res <- pca_plot(vst_mat, cd, title = "TCGA — PCA on vst counts (all samples)")
ggsave(file.path(FIG_DIR, "02_pca_pre_filter.png"), pca_res$plot, width = 9, height = 7)

# =============================================================================
# 7. Identify and exclude PCA outliers
# =============================================================================
outlier_samples <- pca_res$outlier_df$label[pca_res$outlier_df$outlier]
message(sprintf("   PCA outliers (> %d SD on PC1/PC2): %d", PCA_OUTLIER_SD, length(outlier_samples)))

for (s in outlier_samples) {
  add_exclusion(s, cd[s, "project_id"], "pca_outlier",
                sprintf("> %d SD from centroid on PC1 or PC2 (vst-transformed)", PCA_OUTLIER_SD),
                sprintf("|PC| <= %d SD from mean", PCA_OUTLIER_SD))
}

# =============================================================================
# 8. Sample correlation heatmap
# =============================================================================
message(">>> Sample correlation heatmap ...")

# Subsample for visualisation if too many samples
max_heatmap_samples <- 100
if (ncol(vst_mat) > max_heatmap_samples) {
  set.seed(42)
  heatmap_idx <- sample(ncol(vst_mat), max_heatmap_samples)
  heatmap_mat <- vst_mat[, heatmap_idx]
  heatmap_cd <- cd[heatmap_idx, , drop = FALSE]
} else {
  heatmap_mat <- vst_mat
  heatmap_cd <- cd
}

cor_mat <- cor(heatmap_mat)
annotation_col <- data.frame(
  Project = heatmap_cd$project_id,
  Type = heatmap_cd$sample_type,
  row.names = colnames(heatmap_mat)
)
ann_colors <- list(
  Project = c("TCGA-COAD" = "#E41A1C", "TCGA-READ" = "#377EB8"),
  Type = c("Tumor" = "black", "Normal" = "grey90", "Other" = "white")
)

png(file.path(FIG_DIR, "04_sample_correlation_heatmap.png"),
    width = 10, height = 9, units = "in", res = 150)
pheatmap(cor_mat,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE, show_rownames = FALSE,
         main = "TCGA — Sample correlation heatmap (vst)",
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
dev.off()

# =============================================================================
# 9. Missing clinical data summary
# =============================================================================
message(">>> Missing clinical data ...")

clin_all <- bind_rows(
  mutate(clin_COAD, project = "TCGA-COAD"),
  mutate(clin_READ, project = "TCGA-READ")
)

clin_missing <- data.frame(
  variable = colnames(clin_all),
  n_total = nrow(clin_all),
  n_missing = colSums(is.na(clin_all) | clin_all == ""),
  pct_missing = round(100 * colSums(is.na(clin_all) | clin_all == "") / nrow(clin_all), 1)
) %>% arrange(desc(n_missing))

message("   Top 10 most-missing clinical variables:")
for (i in seq_len(min(10, nrow(clin_missing)))) {
  r <- clin_missing[i, ]
  message(sprintf("     %-35s : %3d / %d missing (%.0f%%)",
                  r$variable, r$n_missing, r$n_total, r$pct_missing))
}

# Flag columns with > MAX_NA_CLINICAL missing
drop_cols <- clin_missing$variable[clin_missing$pct_missing > MAX_NA_CLINICAL * 100]
if (length(drop_cols) > 0) {
  message(sprintf("   Columns with > %.0f%% missing (excluded from analysis): %s",
                  MAX_NA_CLINICAL * 100, paste(drop_cols, collapse = ", ")))
}

# Plot
clin_missing_plot <- clin_missing %>%
  arrange(pct_missing) %>%
  slice_tail(n = 30)  # top 30 most-missing

p_clin <- ggplot(clin_missing_plot, aes(x = reorder(variable, pct_missing), y = pct_missing)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "TCGA — Missing clinical data",
       subtitle = sprintf("%d clinical variables (top 30 shown)", nrow(clin_missing)),
       x = "", y = "% missing") +
  theme_minimal()
ggsave(file.path(FIG_DIR, "06_missing_clinical_data.png"), p_clin, width = 8, height = 7)

# =============================================================================
# 10. Apply exclusions and rebuild DDS
# =============================================================================
keep_samples <- setdiff(colnames(dds), exclusions$sample_id)
message(sprintf(">>> Samples retained after QC: %d / %d", length(keep_samples), ncol(dds)))

dds_filt <- dds[, keep_samples]
dds_filt <- estimateSizeFactors(dds_filt)

# vst on filtered set (much faster than rlog)
vsd_filt <- vst(dds_filt, blind = TRUE)
vst_filt_mat <- assay(vsd_filt)

# Post-filter PCA
pca_post <- pca_plot(vst_filt_mat, as.data.frame(colData(dds_filt)),
                     title = sprintf("TCGA — PCA after QC (%d samples)", length(keep_samples)))
ggsave(file.path(FIG_DIR, "03_pca_post_filter.png"), pca_post$plot, width = 9, height = 7)

# =============================================================================
# 11. Save outputs
# =============================================================================
message(">>> Saving outputs ...")

saveRDS(dds_filt, file.path(OUT_DIR, "tcga_dds.rds"))
saveRDS(vsd_filt, file.path(OUT_DIR, "tcga_vsd.rds"))
saveRDS(dds, file.path(OUT_DIR, "tcga_dds_pre_filter.rds"))  # pre-filter for reference

# Write exclusion table (Table S2)
write.csv(exclusions, file.path(OUT_DIR, "tcga_excluded_samples.csv"), row.names = FALSE)

# Write numeric summary
sink(file.path(OUT_DIR, "tcga_qc_summary.txt"))
cat("TCGA Discovery Cohort — QC Summary\n")
cat("==================================\n\n")
cat(sprintf("Raw samples              : %d\n", ncol(counts_raw)))
cat(sprintf("Raw genes                : %d\n", nrow(counts_raw)))
cat(sprintf("  Tumor (TP)             : %d\n", sum(col_data$sample_type == "Tumor")))
cat(sprintf("  Normal (NT)            : %d\n", sum(col_data$sample_type == "Normal")))
cat(sprintf("  Other                  : %d\n", sum(col_data$sample_type == "Other")))
cat(sprintf("  COAD                   : %d\n", sum(col_data$project_id == "TCGA-COAD")))
cat(sprintf("  READ                   : %d\n", sum(col_data$project_id == "TCGA-READ")))
cat(sprintf("\nExcluded samples         : %d\n", nrow(exclusions)))
for (i in seq_len(nrow(exclusions))) {
  cat(sprintf("  - %s: %s\n", exclusions$sample_id[i], exclusions$reason[i]))
}
cat(sprintf("\nRetained samples         : %d\n", ncol(dds_filt)))
cat(sprintf("Retained genes           : %d (filter: >=%d counts in >=%.0f%% samples)\n",
            nrow(dds_filt), MIN_COUNTS, MIN_SAMPLE_PROP * 100))
cat(sprintf("Median library size      : %s\n", format(median(lib_sizes), big.mark = ",")))
cat(sprintf("Size factor range        : %.3f – %.3f\n",
            min(sizeFactors(dds_filt)), max(sizeFactors(dds_filt))))
cat(sprintf("PCA outliers detected    : %d\n", length(outlier_samples)))
cat(sprintf("PCA outlier threshold    : > %d SD from centroid\n", PCA_OUTLIER_SD))
cat(sprintf("Clinical variables dropped: %d (>%.0f%% missing)\n",
            length(drop_cols), MAX_NA_CLINICAL * 100))
cat("\n— End of QC summary —\n")
sink()

message("")
message(">>> QC complete. Outputs:")
message("     ", file.path(OUT_DIR, "tcga_dds.rds"))
message("     ", file.path(OUT_DIR, "tcga_vsd.rds"))
message("     ", file.path(OUT_DIR, "tcga_excluded_samples.csv"))
message("     ", file.path(OUT_DIR, "tcga_qc_summary.txt"))
message("     ", FIG_DIR, "/ (6 figures)")

message("")
message(">>> Proceed to scripts/03_qc_validation.R")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))
