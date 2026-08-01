#!/usr/bin/env Rscript
# =============================================================================
# 04_deg_analysis.R — Differential expression: tumour vs normal (DESeq2)
#
# Purpose:
#   Identify differentially expressed genes between primary CRC tumours and
#   solid-tissue normal samples in TCGA-COAD + TCGA-READ using DESeq2.
#   Apply pre-registered thresholds — no post-hoc cherry-picking.
#
# Inputs:
#   data/processed/tcga_qc/tcga_dds.rds     — DESeqDataSet (filtered, size-factor est.)
#
# Outputs:
#   results/tables/
#   ├── DEG_full_TCGA.csv                   — full gene-level results (all genes)
#   ├── DEG_summary_TCGA.txt                — counts: up/down, total tested
#   └── DEG_sig_TCGA.csv                    — significant only (for convenience)
#   data/processed/tcga_qc/tcga_dds_de.rds  — after DESeq() (for downstream)
#
# Thresholds (pre-registered in logs/decisions.md):
#   |log2FoldChange| > 1
#   adjusted p-value (BH)  < 0.05
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

# -- Thresholds (stated explicitly — never changed post-hoc) -----------------
LOG2FC_CUTOFF <- 1
PADJ_CUTOFF   <- 0.05

# -- Directories -------------------------------------------------------------
IN_DIR  <- file.path("data", "processed", "tcga_qc")
OUT_DIR <- file.path("results", "tables")
FIG_DIR <- file.path("results", "figures", "deg_tcga")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load filtered DESeqDataSet
# =============================================================================
message(">>> Loading filtered DESeqDataSet ...")
dds <- readRDS(file.path(IN_DIR, "tcga_dds.rds"))

# Only tumour vs normal samples for DEG analysis
valid_types <- c("Tumor", "Normal")
dds <- dds[, colData(dds)$sample_type %in% valid_types]
# Re-level so Normal is the reference
colData(dds)$sample_type <- relevel(factor(colData(dds)$sample_type), ref = "Normal")
design(dds) <- ~ sample_type + project_id

message(sprintf("   Samples: %d tumour, %d normal",
                sum(colData(dds)$sample_type == "Tumor"),
                sum(colData(dds)$sample_type == "Normal")))
message(sprintf("   Genes:   %d", nrow(dds)))

# =============================================================================
# 2. DESeq2
# =============================================================================
message(">>> Running DESeq2 ...")
dds <- DESeq(dds)

# =============================================================================
# 3. Extract results
# =============================================================================
message(">>> Extracting results for tumour vs Normal ...")
res <- results(dds, contrast = c("sample_type", "Tumor", "Normal"),
               alpha = PADJ_CUTOFF)

# Convert to data frame
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id")

# Add gene symbols from rowData if available
if ("gene_name" %in% colnames(rowData(dds))) {
  symbol_map <- data.frame(
    gene_id   = rownames(dds),
    gene_name = rowData(dds)$gene_name,
    stringsAsFactors = FALSE
  )
  res_df <- res_df %>%
    left_join(symbol_map, by = "gene_id") %>%
    select(gene_id, gene_name, everything())
} else {
  res_df <- res_df %>%
    mutate(gene_name = gene_id) %>%
    select(gene_id, gene_name, everything())
}

# Order by padj
res_df <- res_df %>% arrange(padj, pvalue)

message(sprintf("   Genes tested          : %d", nrow(res_df)))
message(sprintf("   Significant (padj<%s) : %d", PADJ_CUTOFF, sum(res_df$padj < PADJ_CUTOFF, na.rm = TRUE)))

# =============================================================================
# 4. Apply thresholds
# =============================================================================
res_df <- res_df %>%
  mutate(
    significant = ifelse(
      is.na(padj) | is.na(log2FoldChange), FALSE,
      abs(log2FoldChange) > LOG2FC_CUTOFF & padj < PADJ_CUTOFF
    ),
    direction = case_when(
      significant & log2FoldChange > 0  ~ "Up",
      significant & log2FoldChange < 0  ~ "Down",
      TRUE                              ~ "NS"
    )
  )

n_up   <- sum(res_df$direction == "Up", na.rm = TRUE)
n_down <- sum(res_df$direction == "Down", na.rm = TRUE)
cat(sprintf("   |log2FC| > %.1f  &  padj < %s\n", LOG2FC_CUTOFF, PADJ_CUTOFF))
cat(sprintf("   Up-regulated   : %d\n", n_up))
cat(sprintf("   Down-regulated : %d\n", n_down))
cat(sprintf("   Total DE genes : %d\n", n_up + n_down))

# =============================================================================
# 5. Volcano plot
# =============================================================================
message(">>> Volcano plot ...")
volcano_data <- res_df %>%
  mutate(log10padj = -log10(padj)) %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

p_volc <- ggplot(volcano_data, aes(x = log2FoldChange, y = log10padj, colour = direction)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey60")) +
  geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", alpha = 0.4) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", alpha = 0.4) +
  labs(title = "TCGA — Tumour vs Normal",
       subtitle = sprintf("|log2FC| > %.1f, padj < %s  |  ↑ %d  ↓ %d",
                          LOG2FC_CUTOFF, PADJ_CUTOFF, n_up, n_down),
       x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "01_volcano_plot.png"), p_volc, width = 8, height = 7)

# =============================================================================
# 6. Save outputs
# =============================================================================
message(">>> Saving results ...")

# Full table (all tested genes)
write.csv(res_df, file.path(OUT_DIR, "DEG_full_TCGA.csv"), row.names = FALSE)

# Significant only
sig_df <- res_df %>% filter(significant)
write.csv(sig_df, file.path(OUT_DIR, "DEG_sig_TCGA.csv"), row.names = FALSE)

# Summary text
sink(file.path(OUT_DIR, "DEG_summary_TCGA.txt"))
cat("Differential Expression Summary — TCGA-COAD + TCGA-READ\n")
cat("========================================================\n\n")
cat(sprintf("Comparison       : Tumour vs Normal\n"))
cat(sprintf("Method           : DESeq2 (Wald test)\n"))
cat(sprintf("Filter           : tumour (TP) + solid-tissue normal (NT) only\n"))
cat(sprintf("Genes tested     : %d\n", nrow(res_df)))
cat(sprintf("Significance     : |log2FC| > %.1f, padj < %s\n", LOG2FC_CUTOFF, PADJ_CUTOFF))
cat(sprintf("  Up-regulated   : %d (%.1f%%)\n", n_up, 100 * n_up / nrow(res_df)))
cat(sprintf("  Down-regulated : %d (%.1f%%)\n", n_down, 100 * n_down / nrow(res_df)))
cat(sprintf("  Not significant: %d (%.1f%%)\n",
            sum(res_df$direction == "NS", na.rm = TRUE),
            100 * sum(res_df$direction == "NS", na.rm = TRUE) / nrow(res_df)))
cat(sprintf("  NA (low count) : %d\n", sum(is.na(res_df$pvalue))))
cat("\nUpper quartile of |log2FC| among significant genes:\n")
up_fc <- sig_df$log2FoldChange[sig_df$direction == "Up"]
down_fc <- sig_df$log2FoldChange[sig_df$direction == "Down"]
cat(sprintf("  Up:   Q1 = %.2f, median = %.2f, Q3 = %.2f\n",
            quantile(up_fc, 0.25), median(up_fc), quantile(up_fc, 0.75)))
cat(sprintf("  Down: Q1 = %.2f, median = %.2f, Q3 = %.2f\n",
            quantile(down_fc, 0.25), median(down_fc), quantile(down_fc, 0.75)))
cat("\n— End of DEG summary —\n")
sink()

# Save updated dds for downstream scripts
saveRDS(dds, file.path(IN_DIR, "tcga_dds_de.rds"))

message("")
message(">>> Outputs:")
message("     ", file.path(OUT_DIR, "DEG_full_TCGA.csv"), "  — Full table (",
        nrow(res_df), " genes)")
message("     ", file.path(OUT_DIR, "DEG_sig_TCGA.csv"), "   — Significant only (",
        nrow(sig_df), " genes)")
message("     ", file.path(OUT_DIR, "DEG_summary_TCGA.txt"))
message("     ", file.path(FIG_DIR, "01_volcano_plot.png"))
message(">>> Proceed to scripts/05_csc_marker_focus.R")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info_04.txt"))
