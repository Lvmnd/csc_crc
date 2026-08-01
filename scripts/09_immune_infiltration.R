#!/usr/bin/env Rscript
# =============================================================================
# 09_immune_infiltration.R — Immune cell composition by ssGSEA
#
# Method:
#   Single-sample GSEA (ssGSEA) — rank-based enrichment score per sample per
#   immune cell type. Uses curated gene signatures from published literature
#   (Bindea et al. 2013 Immunity; Charoentong et al. 2017 Cell Reports).
#
#   No external deconvolution package required — the ssGSEA algorithm is
#   implemented manually for full reproducibility and transparency.
#
# Workflow:
#   1. Define immune cell type signatures (10 major types)
#   2. ssGSEA scoring function
#   3. Score TCGA tumour samples (rlog data)
#   4. Score validation cohort (log-normalised data)
#   5. Compare high- vs low-risk groups (boxplots + stats)
#   6. Correlate risk score with immune fractions (heatmap + scatter)
#   7. Save all figures and tables
#
# Inputs:
#   data/processed/tcga_qc/tcga_vsd.rds
#   data/processed/tcga_sig/risk_scores.csv
#   data/processed/geo_qc/GSE107422/*_filtered.rds  (or log2_tpm.rds)
#   data/processed/validation/risk_scores_validation.csv
#
# Outputs:
#   results/tables/immune_infiltration/
#   ├── TCGA_immune_scores.csv
#   ├── Validation_immune_scores.csv
#   ├── TCGA_immune_comparison_stats.txt      — Wilcoxon: high vs low per cell type
#   ├── Validation_immune_comparison_stats.txt
#   └── Immune_correlation_with_risk.txt      — Spearman rho + p per cell type
#
#   results/figures/immune_infiltration/
#   ├── TCGA_01_immune_boxplots.png
#   ├── TCGA_02_immune_heatmap.png
#   ├── TCGA_03_immune_correlation_heatmap.png
#   ├── TCGA_04_risk_vs_immune_scatter.png    (top 4 correlations)
#   ├── Validation_01_immune_boxplots.png
#   └── Validation_02_immune_heatmap.png
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(SummarizedExperiment)
  library(GSVA)
  library(org.Hs.eg.db)
})

# -- Directories -------------------------------------------------------------
TCGA_QC_DIR <- file.path("data", "processed", "tcga_qc")
TCGA_SIG_DIR <- file.path("data", "processed", "tcga_sig")
GEO_QC_DIR  <- file.path("data", "processed", "geo_qc", "GSE107422")
VAL_DIR     <- file.path("data", "processed", "validation")
TAB_DIR     <- file.path("results", "tables", "immune_infiltration")
FIG_DIR     <- file.path("results", "figures", "immune_infiltration")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Immune cell type signatures
# =============================================================================
# Sources:
#   Bindea et al. 2013 (Immunity) — 24 immune cell types
#   Charoentong et al. 2017 (Cell Reports) — immunophenogram
#   Angelova et al. 2015 (Genome Biology) — immune contexture
#
# Each signature contains 15-50 genes. All symbols in UPPERCASE.

immune_signatures <- list(

  `CD8 T cells` = c(
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "GZMA", "GZMB", "GZMK",
    "PRF1", "IFNG", "TNFSF14", "CCL5", "CXCL9", "CXCL10", "CXCL11",
    "EOMES", "TBX21", "CD27", "CD28", "IL2RB", "KLRK1", "NKG7"
  ),

  `CD4 T cells` = c(
    "CD4", "CD3D", "CD3E", "CD3G", "IL2RA", "CD40LG", "ICOS",
    "IL4", "IL5", "IL10", "IL13", "IL21", "CXCR5", "BCL6",
    "STAT6", "GATA3", "MAF", "IL17A", "RORC", "CCR6"
  ),

  `Tregs` = c(
    "FOXP3", "CD4", "IL2RA", "CTLA4", "TNFRSF18", "TNFRSF4",
    "IKZF2", "IKZF4", "CCR8", "ENTPD1", "LGALS3", "TIGIT",
    "LAG3", "PDCD1", "HAVCR2"
  ),

  `NK cells` = c(
    "NKG7", "GNLY", "KLRB1", "KLRD1", "KLRK1", "KIR2DL1",
    "KIR2DL3", "KIR3DL1", "KIR3DL2", "NCR1", "NCR3", "CD160",
    "CD244", "PRF1", "GZMA", "GZMB", "GZMM", "IFNG", "FCGR3A"
  ),

  `B cells` = c(
    "CD19", "CD20", "MS4A1", "CD22", "CD79A", "CD79B", "BLK",
    "BACH2", "PAX5", "EBF1", "IRF8", "POU2F2", "JCHAIN",
    "CXCR5", "CXCR4", "CR2", "FCER2", "SDC1", "MZB1", "TNFRSF17"
  ),

  `Macrophages` = c(
    "CD14", "CD68", "CD163", "FCGR1A", "FCGR2A", "FCGR3A",
    "ITGAM", "ITGAX", "CSF1R", "CSF1", "CCL2", "CCL7", "CCL8",
    "CXCL10", "CXCL11", "IL1B", "IL6", "TNF", "TLR2", "TLR4",
    "TLR8", "MARCO", "MSR1", "CD36"
  ),

  `M1 Macrophages` = c(
    "IL1A", "IL1B", "IL6", "IL12A", "IL12B", "IL23A", "TNF",
    "CXCL9", "CXCL10", "CXCL11", "CCL5", "CCR7", "NOS2",
    "IRF5", "IRF1", "STAT1", "TLR2", "TLR4", "CD80", "CD86",
    "FCGR1A", "FCGR1B"
  ),

  `M2 Macrophages` = c(
    "IL10", "IL4", "IL13", "CCL17", "CCL18", "CCL22", "CCL24",
    "CCL26", "CCR4", "CCR8", "CXCR1", "CXCR2", "CD163", "MRC1",
    "MSR1", "MMP9", "MMP12", "VEGFA", "ARG1", "ARG2", "TGFB1",
    "STAT6", "IRF4", "PPARG", "MYC", "MERTK"
  ),

  `Dendritic cells` = c(
    "CD1A", "CD1B", "CD1C", "CD1D", "CD207", "CD209", "CD80",
    "CD83", "CD86", "CLEC4C", "CLEC9A", "FCER1A", "FLT3",
    "ITGAX", "LAMP3", "LILRA4", "NRP1", "THBD", "XCR1",
    "BATF3", "IRF4", "IRF8", "ZBTB46", "TLR3", "TLR7", "TLR9"
  ),

  `Neutrophils` = c(
    "FCGR3B", "CEACAM8", "CEACAM6", "CSF3R", "CXCR1", "CXCR2",
    "FPR1", "FPR2", "ITGAM", "ITGB2", "MMP8", "MMP9", "MPO",
    "ELANE", "AZU1", "PRTN3", "BPI", "LACT", "LTF", "LCN2",
    "SLC11A1", "TLR4", "TLR8", "NCF1", "NCF2", "CYBB"
  ),

  `Monocytes` = c(
    "CD14", "CD16", "FCGR3A", "LYZ", "S100A8", "S100A9", "S100A12",
    "CCR2", "CX3CR1", "CD33", "CD86", "ITGAM", "CSF1R",
    "CLEC7A", "TLR2", "TLR4", "TLR8", "NLRP3", "IL1B", "TNF",
    "CCL2", "CCL3", "CCL4", "CCL5"
  ),

  `MDSCs` = c(
    "CD33", "CD11B", "ITGAM", "S100A8", "S100A9", "S100A12",
    "ARG1", "ARG2", "NOS2", "IL10", "TGFB1", "CXCR2", "CXCR4",
    "CCL2", "CCL4", "CCL5", "MMP9", "MMP2", "LOX", "COX2",
    "STAT3", "CEBPB", "CSF1R", "CSF3R", "LAG3", "PDCD1", "VEGFA"
  )
)

message(sprintf(">>> Immune signatures loaded: %d cell types", length(immune_signatures)))
for (nm in names(immune_signatures)) {
  message(sprintf("     %-20s : %d genes", nm, length(immune_signatures[[nm]])))
}

# =============================================================================
# 2. ssGSEA scoring function (manual implementation)
# =============================================================================
# Algorithm: Barbie et al. 2009 (Nature)
#   For each sample (expression vector) and each gene set:
#   1. Rank genes by expression (descending)
#   2. Compute enrichment score by walking down the ranked list
#   3. ES = sum(weight of genes in set) - sum(weight of genes not in set)
#
# Returns a normalised enrichment score (range approx -1 to 1)

ssgsea_score <- function(expr_vector, gene_set, all_genes, alpha = 0.25) {
  # expr_vector: named numeric vector of expression values
  # gene_set: character vector of gene symbols in the set
  # all_genes: character vector of all gene symbols (for ranking)
  # alpha: weight parameter (0.25 = default from Barbie et al.)

  # Rank genes by expression (descending)
  ranks <- rank(-expr_vector, ties.method = "average")

  N <- length(all_genes)
  # Indices of genes in the set
  hit_indices <- which(all_genes %in% gene_set)
  Nh <- length(hit_indices)

  if (Nh == 0) return(NA_real_)
  if (Nh < 3) return(NA_real_)  # too few genes for reliable score

  # Enrichment score computation
  step_up   <- sum(abs(ranks[hit_indices])^alpha)
  step_down <- 1 / (N - Nh)

  # Walk the ranked list
  es <- 0
  max_es <- 0
  min_es <- 0

  for (i in seq_len(N)) {
    if (i %in% hit_indices) {
      es <- es + (abs(ranks[i])^alpha) / step_up
    } else {
      es <- es - step_down
    }
    max_es <- max(max_es, es)
    min_es <- min(min_es, es)
  }

  # Normalise
  if (max_es > abs(min_es)) {
    return(max_es)
  } else {
    return(min_es)
  }
}

# Wrapper for matrix (genes x samples)
compute_ssgsea <- function(expr_mat, gene_sets) {
  # expr_mat: matrix with gene symbols as rownames, samples as colnames
  # gene_sets: named list of gene symbol vectors
  # Returns: matrix (cell types x samples) of enrichment scores

  # Ensure uppercase rownames for matching
  all_genes <- toupper(rownames(expr_mat))
  names(all_genes) <- rownames(expr_mat)

  scores <- matrix(NA, nrow = length(gene_sets), ncol = ncol(expr_mat))
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(expr_mat)

  for (ct in seq_along(gene_sets)) {
    gene_set <- toupper(gene_sets[[ct]])
    for (s in seq_len(ncol(expr_mat))) {
      expr_val <- expr_mat[, s]
      names(expr_val) <- all_genes
      # Remove NA
      expr_val <- expr_val[!is.na(expr_val) & names(expr_val) != ""]
      # Take max per gene (handle duplicates)
      if (any(duplicated(names(expr_val)))) {
        expr_val <- tapply(expr_val, names(expr_val), max)
      }
      scores[ct, s] <- ssgsea_score(expr_val, gene_set, names(expr_val))
    }
    if (ct %% 2 == 0) message(sprintf("     Scored: %s", names(gene_sets)[ct]))
  }
  message(sprintf("     Completed: %d cell types", length(gene_sets)))

  scores
}

# =============================================================================
# 3. Helper: immune boxplots with stats
# =============================================================================
plot_immune_boxplots <- function(score_mat, risk_df, cohort_name, fig_path) {
  # score_mat: cell types x samples
  # risk_df: must have sample_id and risk_group columns

  common_s <- intersect(colnames(score_mat), risk_df$sample_id)
  if (length(common_s) < 10) {
    message("     Too few common samples for immune boxplots.")
    return(invisible(NULL))
  }

  score_sub <- score_mat[, common_s, drop = FALSE]
  risk_sub <- risk_df[match(common_s, risk_df$sample_id), ]

  # Melt
  plot_df <- as.data.frame(t(score_sub))
  plot_df$sample_id <- rownames(plot_df)
  plot_df <- pivot_longer(plot_df, -sample_id, names_to = "cell_type", values_to = "enrichment")
  plot_df$risk_group <- risk_sub$risk_group[match(plot_df$sample_id, risk_sub$sample_id)]

  # Remove NA enrichment
  plot_df <- plot_df %>% dplyr::filter(!is.na(enrichment))

  # Compute stats for annotation
  stats_df <- plot_df %>%
    group_by(cell_type) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(enrichment ~ risk_group, data = cur_data())$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      significance = case_when(
        is.na(p_value)  ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE             ~ "ns"
      ),
      label = sprintf("p = %.3f", p_value)
    )

  # Write stats
  stats_df$cohort <- cohort_name
  write.csv(stats_df,
            file.path(TAB_DIR, paste0(cohort_name, "_immune_comparison_stats.csv")),
            row.names = FALSE)

  # Reorder cell types by median enrichment
  med_order <- plot_df %>%
    group_by(cell_type) %>%
    summarise(med = median(enrichment, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(cell_type)
  plot_df$cell_type <- factor(plot_df$cell_type, levels = med_order)

  p <- ggplot(plot_df, aes(x = cell_type, y = enrichment, fill = risk_group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    scale_fill_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C")) +
    labs(title = paste(cohort_name, "— Immune cell enrichment by risk group"),
         x = "", y = "ssGSEA enrichment score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(size = 11)) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5,
                 colour = "black", linewidth = 0.3, alpha = 0.3)

  # Add significance labels
  stats_df <- stats_df[match(med_order, stats_df$cell_type), ]
  y_max <- max(plot_df$enrichment, na.rm = TRUE)
  for (i in seq_len(nrow(stats_df))) {
    if (!is.na(stats_df$p_value[i]) && stats_df$p_value[i] < 0.05) {
      p <- p + annotate("text", x = i, y = y_max * (1 + 0.02 * i),
                         label = stats_df$significance[i], size = 3,
                         colour = "darkred")
    }
  }

  if (!is.null(fig_path)) ggsave(fig_path, p, width = 12, height = 6)
  p
}

# =============================================================================
# 4. Helper: immune heatmap (samples x cell types, annotated by risk group)
# =============================================================================
plot_immune_heatmap <- function(score_mat, risk_df, cohort_name, fig_path) {
  common_s <- intersect(colnames(score_mat), risk_df$sample_id)
  if (length(common_s) < 10) return(invisible(NULL))

  score_sub <- score_mat[, common_s, drop = FALSE]
  risk_sub <- risk_df[match(common_s, risk_df$sample_id), ]

  # Z-score across samples for each cell type
  score_z <- t(scale(t(score_sub)))
  # Cap at [-3, 3] for visualisation
  score_z[score_z > 3] <- 3
  score_z[score_z < -3] <- -3

  # Annotation
  annot_col <- data.frame(
    Risk = factor(risk_sub$risk_group, levels = c("Low", "High")),
    row.names = common_s
  )
  annot_colours <- list(Risk = c(Low = "#377EB8", High = "#E41A1C"))

  # Order by risk group then score
  ord <- order(risk_sub$risk_group, risk_sub$risk_score)
  score_z <- score_z[, ord, drop = FALSE]
  annot_col <- annot_col[ord, , drop = FALSE]

  pheatmap(score_z,
           annotation_col = annot_col,
           annotation_colors = annot_colours,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           main = paste(cohort_name, "— Immune cell enrichment"),
           fontsize_row = 9,
           filename = fig_path,
           width = 8, height = 5)
}

# =============================================================================
# 5. Helper: correlation heatmap (immune vs risk score)
# =============================================================================
plot_immune_correlation <- function(score_mat, risk_df, cohort_name,
                                     fig_path_cor, fig_path_scatter) {
  common_s <- intersect(colnames(score_mat), risk_df$sample_id)
  if (length(common_s) < 10) return(invisible(NULL))

  score_sub <- score_mat[, common_s, drop = FALSE]
  risk_sub <- risk_df[match(common_s, risk_df$sample_id), ]

  # Spearman correlation per cell type with risk score
  cor_df <- data.frame(
    cell_type = rownames(score_sub),
    rho = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(score_sub))) {
    test <- tryCatch(
      cor.test(score_sub[i, ], risk_sub$risk_score, method = "spearman"),
      error = function(e) NULL
    )
    if (!is.null(test)) {
      cor_df$rho[i]     <- test$estimate
      cor_df$p_value[i] <- test$p.value
    }
  }

  cor_df$cohort <- cohort_name
  write.csv(cor_df, file.path(TAB_DIR,
                               paste0(cohort_name, "_immune_risk_correlation.csv")),
            row.names = FALSE)

  # Correlation heatmap
  cor_mat <- as.matrix(cor_df$rho)
  rownames(cor_mat) <- cor_df$cell_type
  colnames(cor_mat) <- "Risk score"

  # P-value annotation
  pval_annot <- ifelse(cor_df$p_value < 0.001, "***",
                       ifelse(cor_df$p_value < 0.01, "**",
                              ifelse(cor_df$p_value < 0.05, "*", "")))

  # Colour by sign
  color_pal <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(50)
  breaks <- seq(-1, 1, length.out = 51)

  png(fig_path_cor, width = 5, height = max(5, nrow(cor_mat) * 0.4), units = "in", res = 150)
  par(mar = c(4, 10, 3, 4))
  image(t(cor_mat[nrow(cor_mat):1, , drop = FALSE]),
        axes = FALSE, col = color_pal, breaks = breaks,
        main = paste(cohort_name, "— Correlation with risk score"))
  axis(2, at = seq(0, 1, length.out = nrow(cor_mat)),
       labels = rev(cor_df$cell_type), las = 2, cex.axis = 0.8)
  axis(1, at = 0.5, labels = "Spearman rho", cex.axis = 0.8)
  # Add text annotations
  for (i in seq_len(nrow(cor_df))) {
    r_idx <- nrow(cor_mat) - i + 1
    text(0.5, (r_idx - 1) / (nrow(cor_mat) - 1),
         sprintf("%.2f %s", cor_df$rho[i], pval_annot[i]),
         cex = 0.7, col = ifelse(cor_df$rho[i] > 0.3, "white", "black"))
  }
  dev.off()

  # Scatter plot for top 4 most correlated cell types
  top_cor <- cor_df[order(abs(cor_df$rho), decreasing = TRUE), ]
  top_cor <- head(top_cor, 4)

  if (nrow(top_cor) > 0) {
    scatter_list <- list()
    for (i in seq_len(nrow(top_cor))) {
      ct <- top_cor$cell_type[i]
      df <- data.frame(
        risk_score     = risk_sub$risk_score,
        immune_score   = score_sub[ct, ],
        cell_type      = ct,
        rho            = sprintf("rho = %.3f", top_cor$rho[i]),
        pval           = sprintf("p = %.3f", top_cor$p_value[i]),
        stringsAsFactors = FALSE
      )
      p <- ggplot(df, aes(x = risk_score, y = immune_score)) +
        geom_point(aes(colour = risk_sub$risk_group), alpha = 0.6, size = 1.5) +
        geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.5) +
        scale_colour_manual(values = c("Low" = "#377EB8", "High" = "#E41A1C")) +
        labs(title = ct, subtitle = paste(df$rho[1], df$pval[1]),
             x = "Risk score", y = "ssGSEA enrichment") +
        theme_minimal() +
        theme(legend.position = "none", plot.title = element_text(size = 10))
      scatter_list[[i]] <- p
    }
    if (length(scatter_list) > 0) {
      p_combined <- gridExtra::arrangeGrob(grobs = scatter_list, ncol = 2,
                                             top = paste(cohort_name,
                                                         "— Immune correlation with risk score"))
      ggsave(fig_path_scatter, p_combined, width = 10, height = 8)
    }
  }

  cor_df
}

# =============================================================================
# 6. TCGA immune infiltration
# =============================================================================
message("")
message("==============================================")
message("  TCGA IMMUNE INFILTRATION                    ")
message("==============================================")
message("")

message(">>> Loading TCGA expression data for ssGSEA ...")

# Try VST data first (best for continuous scoring)
vsd_file <- file.path(TCGA_QC_DIR, "tcga_vsd.rds")
dds_file <- file.path(TCGA_QC_DIR, "tcga_dds_de.rds")
risk_file <- file.path(TCGA_SIG_DIR, "risk_scores.csv")

if (!file.exists(risk_file)) stop("Run 06_signature_construction.R first.")
tcga_risk <- read.csv(risk_file, stringsAsFactors = FALSE)

if (file.exists(vsd_file)) {
  vsd <- readRDS(vsd_file)
  tcga_expr <- assay(vsd)
  message(sprintf("   Using VST-transformed expression: %d x %d",
                  nrow(tcga_expr), ncol(tcga_expr)))
} else if (file.exists(dds_file)) {
  dds <- readRDS(dds_file)
  tcga_expr <- log2(counts(dds, normalized = TRUE) + 1)
  message(sprintf("   Using log2(norm_counts+1): %d x %d",
                  nrow(tcga_expr), ncol(tcga_expr)))
} else {
  stop("No TCGA expression data found. Run 03_qc_discovery.R first.")
}

# Map Ensembl IDs to gene symbols for immune signature matching
ensembl_ids <- gsub("\\..*$", "", rownames(tcga_expr))
gene_symbols <- suppressMessages(mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first"))
tcga_symbols <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)

# Take only tumour samples for immune analysis
tumour_idx <- if ("sample_type" %in% names(colData(vsd))) {
  colData(vsd)$sample_type == "Tumor"
} else {
  rep(TRUE, ncol(tcga_expr))  # all samples if no annotation
}

tcga_tumour <- tcga_expr[, tumour_idx, drop = FALSE]
rownames(tcga_tumour) <- tcga_symbols

# Remove NA and duplicate rownames
tcga_tumour <- tcga_tumour[!is.na(tcga_symbols) & tcga_symbols != "", , drop = FALSE]
# Keep unique gene symbols (keep max expression row)
dup_genes <- unique(rownames(tcga_tumour)[duplicated(rownames(tcga_tumour))])
for (g in dup_genes) {
  idx <- which(rownames(tcga_tumour) == g)
  max_row <- idx[which.max(rowMeans(tcga_tumour[idx, , drop = FALSE]))]
  tcga_tumour <- tcga_tumour[-setdiff(idx, max_row), ]
}
message(sprintf("   TCGA tumour expression: %d unique genes x %d samples",
                nrow(tcga_tumour), ncol(tcga_tumour)))

# Subset to samples with risk groups
common_t <- intersect(colnames(tcga_tumour), tcga_risk$sample_id)
tcga_tumour <- tcga_tumour[, common_t, drop = FALSE]
tcga_risk_sub <- tcga_risk[match(common_t, tcga_risk$sample_id), ]
message(sprintf("   TCGA samples with risk groups: %d", length(common_t)))

# ssGSEA
message(">>> Computing ssGSEA scores (TCGA) ...")
tcga_immune_scores <- compute_ssgsea(tcga_tumour, immune_signatures)

write.csv(as.data.frame(t(tcga_immune_scores)),
          file.path(TAB_DIR, "TCGA_immune_scores.csv"), row.names = TRUE)

# -- TCGA figures -----------------------------------------------------------
message(">>> TCGA — figures ...")

plot_immune_boxplots(tcga_immune_scores, tcga_risk_sub, "TCGA",
                     file.path(FIG_DIR, "TCGA_01_immune_boxplots.png"))

plot_immune_heatmap(tcga_immune_scores, tcga_risk_sub, "TCGA",
                    file.path(FIG_DIR, "TCGA_02_immune_heatmap.png"))

plot_immune_correlation(tcga_immune_scores, tcga_risk_sub, "TCGA",
                        file.path(FIG_DIR, "TCGA_03_immune_correlation_heatmap.png"),
                        file.path(FIG_DIR, "TCGA_04_risk_vs_immune_scatter.png"))

# =============================================================================
# 7. Validation cohort immune infiltration
# =============================================================================
message("")
message("==============================================")
message("  VALIDATION COHORT IMMUNE INFILTRATION       ")
message("==============================================")
message("")

val_risk_file <- file.path(VAL_DIR, "risk_scores_validation.csv")
if (!file.exists(val_risk_file)) {
  message("   Validation risk scores not found. Run 07_validation.R first.")
  message("   Validation immune infiltration skipped.")
} else {
  val_risk <- read.csv(val_risk_file, stringsAsFactors = FALSE)
  message(sprintf("   Validation risk scores loaded: %d patients", nrow(val_risk)))

  # Load validation expression data (prefer log-normalised)
  val_expr <- NULL
  log_file  <- file.path(GEO_QC_DIR, "log_normalized.rds")
  log2_file <- file.path(GEO_QC_DIR, "log2_tpm.rds")
  tpm_file  <- file.path(GEO_QC_DIR, "tpm_filtered.rds")

  if (file.exists(log_file)) {
    val_expr <- readRDS(log_file)
    message(sprintf("   Using log-normalised counts: %d x %d", nrow(val_expr), ncol(val_expr)))
  } else if (file.exists(log2_file)) {
    val_expr <- readRDS(log2_file)
    message(sprintf("   Using log2-TPM: %d x %d", nrow(val_expr), ncol(val_expr)))
  } else if (file.exists(tpm_file)) {
    val_expr <- log2(readRDS(tpm_file) + 1)
    message(sprintf("   Using log2(TPM+1): %d x %d", nrow(val_expr), ncol(val_expr)))
  } else {
    message("   No processed expression data for GSE107422.")
    message("   Validation immune infiltration skipped.")
  }

  if (!is.null(val_expr) && ncol(val_expr) >= 20) {
    # Normalise gene symbols (handle probe IDs etc.)
    val_symbols <- toupper(rownames(val_expr))
    val_expr <- val_expr[!is.na(val_symbols) & val_symbols != "", , drop = FALSE]
    rownames(val_expr) <- val_symbols[!is.na(val_symbols) & val_symbols != ""]

    # Deduplicate
    dup_v <- unique(rownames(val_expr)[duplicated(rownames(val_expr))])
    for (g in dup_v) {
      idx <- which(rownames(val_expr) == g)
      max_row <- idx[which.max(rowMeans(val_expr[idx, , drop = FALSE]))]
      val_expr <- val_expr[-setdiff(idx, max_row), ]
    }
    message(sprintf("   Validation expression (gene symbols): %d x %d",
                    nrow(val_expr), ncol(val_expr)))

    # Match samples to risk groups
    common_v <- intersect(colnames(val_expr), val_risk$sample_id)
    val_expr <- val_expr[, common_v, drop = FALSE]
    val_risk_sub <- val_risk[match(common_v, val_risk$sample_id), ]
    # Use the TCGA-median-based risk groups
    if ("risk_group_z" %in% names(val_risk_sub)) {
      names(val_risk_sub)[names(val_risk_sub) == "risk_group_z"] <- "risk_group"
    }
    message(sprintf("   Validation samples with risk groups: %d", nrow(val_risk_sub)))

    # ssGSEA
    message(">>> Computing ssGSEA scores (GSE107422) ...")
    val_immune_scores <- compute_ssgsea(val_expr, immune_signatures)

    write.csv(as.data.frame(t(val_immune_scores)),
              file.path(TAB_DIR, "Validation_immune_scores.csv"), row.names = TRUE)

    # -- Validation figures ---------------------------------------------------
    message(">>> GSE107422 — figures ...")

    plot_immune_boxplots(val_immune_scores, val_risk_sub, "GSE107422",
                         file.path(FIG_DIR, "Validation_01_immune_boxplots.png"))

    plot_immune_heatmap(val_immune_scores, val_risk_sub, "GSE107422",
                        file.path(FIG_DIR, "Validation_02_immune_heatmap.png"))

    plot_immune_correlation(val_immune_scores, val_risk_sub, "GSE107422",
                            file.path(FIG_DIR,
                                      "Validation_03_immune_correlation_heatmap.png"),
                            file.path(FIG_DIR,
                                      "Validation_04_risk_vs_immune_scatter.png"))
  }
}

# =============================================================================
# 8. Summary
# =============================================================================
message("")
message("==============================================")
message("   IMMUNE INFILTRATION COMPLETE               ")
message("==============================================")
message("")
message("   Method: ssGSEA (manual implementation, Barbie et al. 2009)")
message("   Signatures: curated from Bindea et al. 2013, Charoentong et al. 2017")
message(sprintf("   Cell types: %d", length(immune_signatures)))
message("")
message(">>> Outputs:")
message("     ", TAB_DIR, "/ (tables + correlation stats)")
message("     ", FIG_DIR, "/ (figures)")

writeLines(capture.output(sessionInfo()), file.path(TAB_DIR, "session_info_09.txt"))
