#!/usr/bin/env Rscript
# =============================================================================
# 08_enrichment.R — Functional enrichment (GSEA + ORA) of risk-group DEGs
#
# HARD RULES:
#   - Report honestly whether Wnt, Notch, Hedgehog, stemness pathways emerge.
#   - Do NOT cherry-pick significant pathways; report all tested.
#   - Negative result (no stemness enrichment) is scientifically informative.
#
# Workflow:
#   1. Load TCGA expression + risk groups
#   2. DESeq2: high-risk vs low-risk (tumour-only)
#   3. Pre-ranked GSEA (KEGG, GO-BP, Hallmark)
#   4. ORA on significant DEGs (KEGG, GO-BP, Hallmark)
#   5. Repeat for validation cohort (GSE107422) if possible
#   6. Honest stemness-pathway audit
#   7. Figures: dot plot, ridge plot, GSEA running-score, cnet
#
# Inputs:
#   data/processed/tcga_qc/tcga_dds_de.rds
#   data/processed/tcga_sig/risk_scores.csv
#   data/processed/geo_qc/GSE107422/[counts|tpm]_filtered.rds
#   data/processed/validation/risk_scores_validation.csv
#
# Outputs:
#   results/tables/enrichment/
#   ├── TCGA_DEG_high_vs_low.csv
#   ├── TCGA_GSEA_results.csv
#   ├── TCGA_ORA_results.csv
#   ├── TCGA_stemness_audit.txt        ← Honest report on Wnt/Notch/Hedgehog
#   └── Validation_stemness_audit.txt
#
#   results/figures/enrichment/
#   ├── TCGA_01_gsea_dotplot.png
#   ├── TCGA_02_gsea_ridgeplot.png
#   ├── TCGA_03_gsea_top_pathways.png   (individual running-score plots)
#   ├── TCGA_04_ora_dotplot.png
#   ├── TCGA_05_stemness_heatmap.png    (expression of core stemness genes)
#   └── Validation_* (same palette, if data available)
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# -- Optional: msigdbr for Hallmark gene sets ---------------------------------
has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
if (has_msigdbr) {
  library(msigdbr)
  message("   msigdbr available — Hallmark gene sets will be included.")
} else {
  message("   msigdbr not available — Hallmark analysis skipped.")
  message("   Install with: conda install -c conda-forge r-msigdbr")
}

# -- Directories -------------------------------------------------------------
TCGA_QC_DIR  <- file.path("data", "processed", "tcga_qc")
TCGA_SIG_DIR <- file.path("data", "processed", "tcga_sig")
GEO_QC_DIR   <- file.path("data", "processed", "geo_qc", "GSE107422")
VAL_DIR      <- file.path("data", "processed", "validation")
TAB_DIR      <- file.path("results", "tables", "enrichment")
FIG_DIR      <- file.path("results", "figures", "enrichment")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# -- Stemness-related pathway keywords (for honest audit) --------------------
stemness_keywords <- c("wnt", "notch", "hedgehog", "stem cell", "stemness",
                        "self-renewal", "embryonic", "pluripotency",
                        "hippo", "tgf.beta", "nanog", "oct4", "sox2",
                        "epithelial.mesenchymal", "emt")

# =============================================================================
# 1. Helper functions
# =============================================================================

# --- ssGSEA-like rank-based expression heatmap for a gene set ----------------
plot_stemness_heatmap <- function(expr_mat, gene_set, risk_df, title, fig_path) {
  # expr_mat: genes x samples (log2-normalised)
  # gene_set: vector of gene symbols
  # risk_df: data.frame with sample_id and risk_group

  common_genes <- intersect(toupper(rownames(expr_mat)), toupper(gene_set))
  if (length(common_genes) < 3) {
    message(sprintf("     Too few stemness genes matched (%d), skipping heatmap.",
                    length(common_genes)))
    return(invisible(NULL))
  }

  # Subset and order expression
  idx <- which(toupper(rownames(expr_mat)) %in% toupper(gene_set))
  sub_expr <- expr_mat[idx, , drop = FALSE]
  rownames(sub_expr) <- rownames(expr_mat)[idx]

  # Match to risk groups
  common_s <- intersect(colnames(sub_expr), risk_df$sample_id)
  sub_expr <- sub_expr[, common_s, drop = FALSE]
  risk_df <- risk_df[match(common_s, risk_df$sample_id), ]

  # Order by risk score
  ord <- order(risk_df$risk_score)
  sub_expr <- sub_expr[, ord, drop = FALSE]
  risk_df <- risk_df[ord, ]

  # Z-score per gene
  sub_expr <- t(scale(t(sub_expr)))

  # Melt for ggplot
  plot_df <- as.data.frame(sub_expr)
  plot_df$gene <- rownames(plot_df)
  plot_df <- pivot_longer(plot_df, -gene, names_to = "sample", values_to = "expression")

  plot_df$risk_group <- risk_df$risk_group[match(plot_df$sample, risk_df$sample_id)]

  p <- ggplot(plot_df, aes(x = sample, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                          midpoint = 0, name = "z-score") +
    labs(title = title, x = "Sample", y = "Gene") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank())

  if (!is.null(fig_path)) ggsave(fig_path, p, width = 10, height = max(4, nrow(sub_expr) * 0.35))
  p
}

# --- Honest stemness pathway audit ------------------------------------------
audit_stemness_pathways <- function(gsea_results, kegg_results, go_results,
                                     hallmark_results = NULL, prefix = "TCGA") {

  # Collect all pathways containing stemness keywords
  all_pathways <- data.frame()

  if (!is.null(gsea_results$kegg)) {
    kegg_df <- as.data.frame(gsea_results$kegg)
    if (nrow(kegg_df) > 0) {
      kegg_df$database <- "KEGG"
      all_pathways <- bind_rows(all_pathways, kegg_df)
    }
  }
  if (!is.null(gsea_results$go)) {
    go_df <- as.data.frame(gsea_results$go)
    if (nrow(go_df) > 0) {
      go_df$database <- "GO-BP"
      all_pathways <- bind_rows(all_pathways, go_df)
    }
  }
  if (!is.null(hallmark_results)) {
    hm_df <- as.data.frame(hallmark_results)
    if (nrow(hm_df) > 0) {
      hm_df$database <- "Hallmark"
      all_pathways <- bind_rows(all_pathways, hm_df)
    }
  }

  if (nrow(all_pathways) == 0) {
    return("No enrichment results available for audit.")
  }

  # Search for stemness keywords in pathway descriptions
  stemness_hits <- all_pathways %>%
    dplyr::filter(grepl(paste(stemness_keywords, collapse = "|"),
                 Description, ignore.case = TRUE))

  audit_lines <- character()
  audit_lines <- c(audit_lines,
    sprintf("=== Stemness Pathway Audit — %s ===\n", prefix),
    sprintf("Method: GSEA (pre-ranked by log2FC high-risk vs low-risk)\n"),
    sprintf("Databases searched: %s\n",
            paste(unique(all_pathways$database), collapse = ", ")),
    sprintf("Total pathways tested: %d\n", nrow(all_pathways)),
    sprintf("Total significant (padj < 0.05): %d\n",
            sum(all_pathways$p.adjust < 0.05, na.rm = TRUE)),
    "")

  if (nrow(stemness_hits) > 0) {
    stemness_sig <- stemness_hits %>% dplyr::filter(p.adjust < 0.05)
    audit_lines <- c(audit_lines,
      sprintf("Stemness-related pathways identified: %d\n", nrow(stemness_hits)),
      sprintf("Stemness-related pathways SIGNIFICANT (padj < 0.05): %d\n",
              nrow(stemness_sig)),
      "")
    if (nrow(stemness_sig) > 0) {
      audit_lines <- c(audit_lines, "Significant stemness pathways:")
      for (i in seq_len(nrow(stemness_sig))) {
        audit_lines <- c(audit_lines,
          sprintf("  [%s] %s — NES = %+.3f, padj = %.2e",
                  stemness_sig$database[i],
                  stemness_sig$Description[i],
                  stemness_sig$NES[i],
                  stemness_sig$p.adjust[i]))
      }
    } else {
      audit_lines <- c(audit_lines,
        "No stemness-related pathways reached statistical significance.")
      audit_lines <- c(audit_lines, "Non-significant stemness-related pathways detected:")
      for (i in seq_len(nrow(stemness_hits))) {
        audit_lines <- c(audit_lines,
          sprintf("  [%s] %s — NES = %+.3f, pval = %.2e",
                  stemness_hits$database[i],
                  stemness_hits$Description[i],
                  stemness_hits$NES[i],
                  stemness_hits$pvalue[i]))
      }
    }
  } else {
    audit_lines <- c(audit_lines,
      "No stemness-related pathways (Wnt, Notch, Hedgehog, EMT, pluripotency)")
    audit_lines <- c(audit_lines,
      "were detected among all enriched pathways.")
    audit_lines <- c(audit_lines,
      "This constitutes a NEGATIVE result for stemness enrichment.")
  }

  audit_lines <- c(audit_lines, "",
    "Interpretation:",
    "  If stemness pathways ARE enriched: supports the hypothesis that the",
    "  CSC signature captures stem-like tumour subpopulations.",
    "  If stemness pathways are NOT enriched: the signature may capture",
    "  non-stemness biological processes (proliferation, invasion, immune",
    "  evasion) or the stemness signal is not captured by bulk RNA-seq.",
    "",
    paste(rep("=", 60), collapse = ""))

  audit_lines
}

# --- GSEA dot plot wrapper --------------------------------------------------
plot_gsea_dot <- function(gsea_obj, title, fig_path, n_show = 20) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("     No significant results to plot.")
    return(invisible(NULL))
  }
  p <- dotplot(gsea_obj, showCategory = n_show,
                title = title, font.size = 10) +
    theme(plot.title = element_text(size = 12))
  if (!is.null(fig_path)) ggsave(fig_path, p, width = 10, height = max(6, n_show * 0.35))
  p
}

# --- GSEA ridge plot wrapper -------------------------------------------------
plot_gsea_ridge <- function(gsea_obj, title, fig_path, n_show = 20) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("     No significant results to plot.")
    return(invisible(NULL))
  }
  p <- ridgeplot(gsea_obj, showCategory = n_show,
                  title = title, font.size = 10) +
    theme(plot.title = element_text(size = 12))
  if (!is.null(fig_path)) ggsave(fig_path, p, width = 11, height = max(6, n_show * 0.35))
  p
}

# --- GSEA running-score plots for top pathways --------------------------------
plot_top_gsea_running <- function(gsea_obj, title, fig_path, n_show = 6) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("     No significant results to plot.")
    return(invisible(NULL))
  }
  p <- gseaplot2(gsea_obj, geneSetID = seq_len(min(n_show, length(gsea_obj@ID))),
                  title = title, subplots = 1:2, pvalue_table = TRUE,
                  ES_geom = "line")
  if (!is.null(fig_path)) ggsave(fig_path, p, width = 12, height = 9)
  p
}

# =============================================================================
# 2. Load TCGA data
# =============================================================================
message("==============================================")
message("  TCGA ENRICHMENT ANALYSIS                    ")
message("==============================================")
message("")

message(">>> Loading TCGA expression + risk groups ...")

dds_file <- file.path(TCGA_QC_DIR, "tcga_dds_de.rds")
if (!file.exists(dds_file)) stop("Run 03_qc_discovery.R first.")
dds <- readRDS(dds_file)

risk_file <- file.path(TCGA_SIG_DIR, "risk_scores.csv")
if (!file.exists(risk_file)) stop("Run 06_signature_construction.R first.")
tcga_risk <- read.csv(risk_file, stringsAsFactors = FALSE)

# Subset to tumour samples for risk-group comparison
tumour_samples <- colnames(dds)[colData(dds)$sample_type == "Tumor"]
risk_tumour <- tcga_risk %>% dplyr::filter(sample_id %in% tumour_samples)
message(sprintf("   Tumour samples with risk groups: %d", nrow(risk_tumour)))
message(sprintf("     High-risk: %d, Low-risk: %d",
                sum(risk_tumour$risk_group == "High"),
                sum(risk_tumour$risk_group == "Low")))

# Check balance
n_high <- sum(risk_tumour$risk_group == "High")
n_low  <- sum(risk_tumour$risk_group == "Low")
if (n_high < 5 || n_low < 5) {
  stop("Risk groups too small for DEG analysis.")
}

# =============================================================================
# 3. DESeq2: high-risk vs low-risk
# =============================================================================
message(">>> DESeq2: high-risk vs low-risk (tumour-only) ...")

# Subset dds to tumour samples with risk group annotation
dds_tumour <- dds[, colnames(dds) %in% risk_tumour$sample_id]

# Ensure order matches
dds_tumour <- dds_tumour[, risk_tumour$sample_id]
colData(dds_tumour)$risk_group <- factor(risk_tumour$risk_group,
                                          levels = c("Low", "High"))

design(dds_tumour) <- ~ risk_group

# Re-run DESeq on this subset (size factors already estimated)
dds_tumour <- DESeq(dds_tumour)

# Extract results
res <- results(dds_tumour, contrast = c("risk_group", "High", "Low"),
                alpha = 0.05)
res <- res[order(res$pvalue), ]

# DEG table
deg_df <- as.data.frame(res)
deg_df$gene_id <- rownames(deg_df)

# Add gene symbols
if ("gene_name" %in% colnames(rowData(dds_tumour))) {
  deg_df$gene_symbol <- rowData(dds_tumour)[rownames(deg_df), "gene_name"]
} else {
  deg_df$gene_symbol <- rownames(deg_df)
}

write.csv(deg_df, file.path(TAB_DIR, "TCGA_DEG_high_vs_low.csv"), row.names = FALSE)

n_up   <- sum(deg_df$padj < 0.05 & deg_df$log2FoldChange > 1, na.rm = TRUE)
n_down <- sum(deg_df$padj < 0.05 & deg_df$log2FoldChange < -1, na.rm = TRUE)
message(sprintf("   DEGs (|log2FC| > 1, padj < 0.05): %d up, %d down", n_up, n_down))

# =============================================================================
# 4. Pre-ranked GSEA
# =============================================================================
message(">>> Pre-ranked GSEA ...")

# Rank genes by: sign(log2FC) * -log10(pvalue)  — robust ranking metric
deg_df <- deg_df[!is.na(deg_df$pvalue), ]
rank_metric <- sign(deg_df$log2FoldChange) * -log10(deg_df$pvalue)
names(rank_metric) <- deg_df$gene_symbol

# Remove duplicates (keep max absolute value)
dup_idx <- duplicated(names(rank_metric))
if (any(dup_idx)) {
  dup_genes <- unique(names(rank_metric)[dup_idx])
  for (g in dup_genes) {
    pos <- which(names(rank_metric) == g)
    rank_metric[pos[1]] <- rank_metric[pos[which.max(abs(rank_metric[pos]))]]
  }
  rank_metric <- rank_metric[!duplicated(names(rank_metric))]
}

# Sort descending
rank_metric <- sort(rank_metric, decreasing = TRUE)
message(sprintf("   Ranked genes: %d", length(rank_metric)))

# -- 4a. GSEA — KEGG ---------------------------------------------------------
message("   GSEA — KEGG ...")
gsea_kegg <- tryCatch({
  gseKEGG(
    geneList     = rank_metric,
    organism     = "hsa",
    keyType      = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    nPermSimple  = 10000,
    seed         = 42,
    verbose      = FALSE
  )
}, error = function(e) {
  message("     KEGG GSEA failed: ", e$message)
  NULL
})

# -- 4b. GSEA — GO-BP --------------------------------------------------------
message("   GSEA — GO-BP ...")
gsea_go <- tryCatch({
  gseGO(
    geneList     = rank_metric,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    nPermSimple  = 10000,
    seed         = 42,
    verbose      = FALSE
  )
}, error = function(e) {
  message("     GO-BP GSEA failed: ", e$message)
  NULL
})

# -- 4c. GSEA — Hallmark (if msigdbr available) ------------------------------
gsea_hallmark <- NULL
if (has_msigdbr) {
  message("   GSEA — Hallmark ...")
  hm_df <- msigdbr(species = "Homo sapiens", category = "H")
  hallmark_list <- split(hm_df$gene_symbol, hm_df$gs_name)

  gsea_hallmark <- tryCatch({
    GSEA(
      geneList     = rank_metric,
      TERM2GENE    = data.frame(term = hm_df$gs_name, gene = hm_df$gene_symbol,
                                 stringsAsFactors = FALSE),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      nPermSimple  = 10000,
      seed         = 42,
      verbose      = FALSE
    )
  }, error = function(e) {
    message("     Hallmark GSEA failed: ", e$message)
    NULL
  })
}

# -- Save GSEA results -------------------------------------------------------
gsea_list <- list(kegg = gsea_kegg, go = gsea_go, hallmark = gsea_hallmark)

for (nm in names(gsea_list)) {
  if (!is.null(gsea_list[[nm]])) {
    df <- as.data.frame(gsea_list[[nm]])
    write.csv(df, file.path(TAB_DIR, paste0("TCGA_GSEA_", nm, ".csv")), row.names = FALSE)
  }
}

message(sprintf("   GSEA significant: KEGG=%d, GO-BP=%d, Hallmark=%d",
                if (!is.null(gsea_kegg)) nrow(as.data.frame(gsea_kegg)) else 0,
                if (!is.null(gsea_go)) nrow(as.data.frame(gsea_go)) else 0,
                if (!is.null(gsea_hallmark)) nrow(as.data.frame(gsea_hallmark)) else 0))

# =============================================================================
# 5. Over-representation analysis (ORA) on significant DEGs
# =============================================================================
message(">>> ORA on significant DEGs ...")

sig_up   <- deg_df$gene_symbol[which(deg_df$padj < 0.05 & deg_df$log2FoldChange > 0)]
sig_down <- deg_df$gene_symbol[which(deg_df$padj < 0.05 & deg_df$log2FoldChange < 0)]

# Run ORA on combined DEG set
sig_genes_all <- unique(c(sig_up, sig_down))
message(sprintf("   Significant DEGs (padj < 0.05, any FC): %d", length(sig_genes_all)))

if (length(sig_genes_all) >= 5) {
  # KEGG ORA
  ora_kegg <- tryCatch({
    enrichKEGG(gene = sig_genes_all, organism = "hsa",
               keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH")
  }, error = function(e) NULL)

  # GO-BP ORA
  ora_go <- tryCatch({
    enrichGO(gene = sig_genes_all, OrgDb = org.Hs.eg.db,
              keyType = "SYMBOL", ont = "BP",
              pvalueCutoff = 0.05, pAdjustMethod = "BH")
  }, error = function(e) NULL)

  # Hallmark ORA (if available)
  ora_hallmark <- NULL
  if (has_msigdbr) {
    ora_hallmark <- tryCatch({
      enricher(gene = sig_genes_all, TERM2GENE = hm_df[, c("gs_name", "gene_symbol")],
               pvalueCutoff = 0.05, pAdjustMethod = "BH")
    }, error = function(e) NULL)
  }
} else {
  message("   Too few DEGs for meaningful ORA. Skipping.")
  ora_kegg <- ora_go <- ora_hallmark <- NULL
}

# Save ORA
ora_list <- list(kegg = ora_kegg, go = ora_go, hallmark = ora_hallmark)
for (nm in names(ora_list)) {
  if (!is.null(ora_list[[nm]])) {
    df <- as.data.frame(ora_list[[nm]])
    if (nrow(df) > 0) {
      write.csv(df, file.path(TAB_DIR, paste0("TCGA_ORA_", nm, ".csv")),
                row.names = FALSE)
    }
  }
}

message(sprintf("   ORA significant: KEGG=%d, GO-BP=%d, Hallmark=%d",
                if (!is.null(ora_kegg)) nrow(as.data.frame(ora_kegg)) else 0,
                if (!is.null(ora_go)) nrow(as.data.frame(ora_go)) else 0,
                if (!is.null(ora_hallmark)) nrow(as.data.frame(ora_hallmark)) else 0))

# =============================================================================
# 6. TCGA — figures
# =============================================================================
message(">>> TCGA — figures ...")

# GSEA dot plots
plot_gsea_dot(gsea_kegg, "TCGA — KEGG GSEA: High-risk vs Low-risk",
              file.path(FIG_DIR, "TCGA_01a_gsea_dotplot_KEGG.png"))
plot_gsea_dot(gsea_go, "TCGA — GO-BP GSEA: High-risk vs Low-risk",
              file.path(FIG_DIR, "TCGA_01b_gsea_dotplot_GO.png"))
if (!is.null(gsea_hallmark)) {
  plot_gsea_dot(gsea_hallmark, "TCGA — Hallmark GSEA: High-risk vs Low-risk",
                file.path(FIG_DIR, "TCGA_01c_gsea_dotplot_Hallmark.png"))
}

# GSEA ridge plots
plot_gsea_ridge(gsea_kegg, "TCGA — KEGG GSEA: High-risk vs Low-risk",
                file.path(FIG_DIR, "TCGA_02a_gsea_ridgeplot_KEGG.png"))
plot_gsea_ridge(gsea_go, "TCGA — GO-BP GSEA: High-risk vs Low-risk",
                file.path(FIG_DIR, "TCGA_02b_gsea_ridgeplot_GO.png"))

# Top running-score plots (KEGG)
if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  plot_top_gsea_running(gsea_kegg, "TCGA — KEGG: Top pathways",
                        file.path(FIG_DIR, "TCGA_03a_gsea_running_KEGG.png"))
}
if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  plot_top_gsea_running(gsea_go, "TCGA — GO-BP: Top pathways",
                        file.path(FIG_DIR, "TCGA_03b_gsea_running_GO.png"))
}

# ORA dot plots
if (!is.null(ora_kegg) && nrow(as.data.frame(ora_kegg)) > 0) {
  p <- dotplot(ora_kegg, showCategory = 20,
               title = "TCGA — ORA KEGG: High-risk DEGs") +
    theme(plot.title = element_text(size = 12))
  ggsave(file.path(FIG_DIR, "TCGA_04a_ora_dotplot_KEGG.png"), p, width = 9, height = 7)
}
if (!is.null(ora_go) && nrow(as.data.frame(ora_go)) > 0) {
  p <- dotplot(ora_go, showCategory = 20,
               title = "TCGA — ORA GO-BP: High-risk DEGs") +
    theme(plot.title = element_text(size = 12))
  ggsave(file.path(FIG_DIR, "TCGA_04b_ora_dotplot_GO.png"), p, width = 9, height = 7)
}

# =============================================================================
# 7. Stemness pathway audit
# =============================================================================
message(">>> Stemness pathway audit (TCGA) ...")

tcga_audit <- audit_stemness_pathways(
  gsea_list, ora_kegg, ora_go, gsea_hallmark, prefix = "TCGA"
)
writeLines(tcga_audit, file.path(TAB_DIR, "TCGA_stemness_audit.txt"))

# Print summary to console
for (l in tcga_audit) message(l)

# Stemness heatmap
core_stemness_genes <- c("WNT1", "WNT2", "WNT3", "WNT3A", "WNT5A", "WNT5B",
                          "CTNNB1", "MYC", "CCND1", "AXIN2", "LEF1", "TCF7",
                          "NOTCH1", "NOTCH2", "NOTCH3", "JAG1", "DLL1", "HES1",
                          "DHH", "IHH", "SHH", "GLI1", "GLI2", "PTCH1", "SMO",
                          "NANOG", "POU5F1", "SOX2", "KLF4", "MYC", "ZFP42",
                          "CD44", "PROM1", "LGR5", "EPCAM", "ALDH1A1")

# Use TCGA vsd data (tumour-only) for heatmap
vsd_file <- file.path(TCGA_QC_DIR, "tcga_vsd.rds")
if (file.exists(vsd_file)) {
  vsd <- readRDS(vsd_file)
  vsd_tumour <- assay(vsd)[, colnames(vsd) %in% risk_tumour$sample_id, drop = FALSE]
  message(sprintf("   VST data loaded: %d x %d", nrow(vsd_tumour), ncol(vsd_tumour)))

  # Map ENSEMBL IDs to symbols via org.Hs.eg.db
  ensembls_vsd <- gsub("\\..*$", "", rownames(vsd_tumour))
  symbol_map_vsd <- suppressMessages(mapIds(org.Hs.eg.db, keys = ensembls_vsd, keytype = "ENSEMBL",
                       column = "SYMBOL", multiVals = "first"))
  vsd_symbols <- ifelse(is.na(symbol_map_vsd), rownames(vsd_tumour), symbol_map_vsd)
  vsd_df <- as.data.frame(vsd_tumour)
  vsd_df <- vsd_df[!is.na(vsd_symbols) & !duplicated(vsd_symbols), , drop = FALSE]
  rownames(vsd_df) <- vsd_symbols[!is.na(vsd_symbols) & !duplicated(vsd_symbols)]

  plot_stemness_heatmap(
    as.matrix(vsd_df), core_stemness_genes,
    risk_tumour,
    "TCGA — Core stemness gene expression (z-score)",
    file.path(FIG_DIR, "TCGA_05_stemness_heatmap.png")
  )
}

# =============================================================================
# 8. Validation cohort (GSE107422) — enrichment if possible
# =============================================================================
message("")
message("==============================================")
message("  VALIDATION COHORT ENRICHMENT                ")
message("==============================================")
message("")

val_risk_file <- file.path(VAL_DIR, "risk_scores_validation.csv")
if (file.exists(val_risk_file)) {
  val_risk <- read.csv(val_risk_file, stringsAsFactors = FALSE)
  message(sprintf("   Validation risk scores loaded: %d patients", nrow(val_risk)))

  # Try to load count data
  val_counts_file <- file.path(GEO_QC_DIR, "counts_filtered.rds")
  val_tpm_file    <- file.path(GEO_QC_DIR, "tpm_filtered.rds")

  val_expr <- NULL
  if (file.exists(val_counts_file)) {
    val_expr <- readRDS(val_counts_file)
    message(sprintf("   Validation count data: %d x %d", nrow(val_expr), ncol(val_expr)))
  } else if (file.exists(val_tpm_file)) {
    val_expr <- readRDS(val_tpm_file)
    message(sprintf("   Validation TPM data: %d x %d", nrow(val_expr), ncol(val_expr)))
    val_expr <- log2(val_expr + 1)
    message("   (log2-transformed for ranking)")
  }

  if (!is.null(val_expr) && ncol(val_expr) >= 20) {
    # Get risk groups
    val_risk <- val_risk %>% dplyr::filter(sample_id %in% colnames(val_expr))

    # For count data, try DESeq2; for TPM, use limma-style approach
    if (file.exists(val_counts_file)) {
      message(">>> DESeq2: high-risk vs low-risk (GSE107422) ...")
      dds_val <- DESeqDataSetFromMatrix(
        countData = val_expr[, val_risk$sample_id],
        colData = val_risk,
        design = ~ risk_group_z   # TCGA-median cutoff
      )
      dds_val <- DESeq(dds_val)
      res_val <- results(dds_val, contrast = c("risk_group_z", "High", "Low"),
                          alpha = 0.05)
      res_val <- res_val[order(res_val$pvalue), ]

      # Rank metric
      deg_val <- as.data.frame(res_val)
      deg_val <- deg_val[!is.na(deg_val$pvalue), ]
      rank_val <- sign(deg_val$log2FoldChange) * -log10(deg_val$pvalue)
      names(rank_val) <- rownames(deg_val)
      # Try to map to symbols
      names(rank_val) <- toupper(names(rank_val))
      rank_val <- sort(rank_val[!duplicated(names(rank_val))], decreasing = TRUE)

      if (length(rank_val) >= 50) {
        # GSEA — just KEGG (limited power expected)
        gsea_val_kegg <- tryCatch({
          gseKEGG(geneList = rank_val, organism = "hsa",
                   pvalueCutoff = 0.1, pAdjustMethod = "BH",
                   nPermSimple = 1000, seed = 42, verbose = FALSE)
        }, error = function(e) NULL)

        if (!is.null(gsea_val_kegg) && nrow(as.data.frame(gsea_val_kegg)) > 0) {
          write.csv(as.data.frame(gsea_val_kegg),
                    file.path(TAB_DIR, "Validation_GSEA_KEGG.csv"), row.names = FALSE)
          plot_gsea_dot(gsea_val_kegg, "Validation — KEGG GSEA: High vs Low risk",
                        file.path(FIG_DIR, "Validation_01_gsea_dotplot_KEGG.png"))
        } else {
          message("   No significant KEGG pathways in validation cohort.")
        }

        # Stemness audit for validation
        val_audit <- audit_stemness_pathways(
          list(kegg = gsea_val_kegg, go = NULL, hallmark = NULL),
          NULL, NULL, NULL, prefix = "GSE107422"
        )
        writeLines(val_audit, file.path(TAB_DIR, "Validation_stemness_audit.txt"))
      } else {
        message("   Too few genes for GSEA in validation cohort.")
      }
    }

    # Stemness heatmap for validation
    if (!is.null(val_expr)) {
      # Try to get gene symbols from rownames
      val_heatmap_expr <- val_expr[, val_risk$sample_id, drop = FALSE]
      rownames(val_heatmap_expr) <- toupper(rownames(val_heatmap_expr))

      # Subset to core stemness genes
      stem_in_val <- intersect(core_stemness_genes, rownames(val_heatmap_expr))
      if (length(stem_in_val) >= 3) {
        val_stem_expr <- val_heatmap_expr[stem_in_val, , drop = FALSE]
        # Reorder by risk score
        val_risk_ordered <- val_risk[match(colnames(val_stem_expr), val_risk$sample_id), ]
        ord <- order(val_risk_ordered$risk_score_z)
        val_stem_expr <- val_stem_expr[, ord, drop = FALSE]
        val_risk_ordered <- val_risk_ordered[ord, ]

        # Z-score per gene
        val_stem_z <- t(scale(t(val_stem_expr)))

        # Plot
        plot_df <- as.data.frame(val_stem_z)
        plot_df$gene <- rownames(plot_df)
        plot_df <- pivot_longer(plot_df, -gene, names_to = "sample", values_to = "expression")
        plot_df$risk_group <- val_risk_ordered$risk_group_z[match(plot_df$sample,
                                                                   val_risk_ordered$sample_id)]

        p <- ggplot(plot_df, aes(x = sample, y = gene, fill = expression)) +
          geom_tile() +
          scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                                midpoint = 0, name = "z-score") +
          labs(title = "GSE107422 — Core stemness gene expression",
               x = "Sample", y = "Gene") +
          theme_minimal() +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                panel.grid = element_blank())
        ggsave(file.path(FIG_DIR, "Validation_05_stemness_heatmap.png"),
               p, width = 10, height = max(4, length(stem_in_val) * 0.35))
      }
    }
  } else {
    message("   Validation cohort expression data unavailable or too small for enrichment.")
  }
} else {
  message("   Validation risk scores not found. Run 07_validation.R first.")
  message("   Validation enrichment skipped.")
}

# =============================================================================
# 9. Final summary
# =============================================================================
message("")
message("==============================================")
message("   ENRICHMENT ANALYSIS COMPLETE               ")
message("==============================================")
message("")
message(sprintf("  TCGA DEGs (|log2FC|>1, padj<0.05): %d up, %d down", n_up, n_down))
message(sprintf("  GSEA significant: KEGG=%d  GO-BP=%d  Hallmark=%d",
                if (!is.null(gsea_kegg)) nrow(as.data.frame(gsea_kegg)) else 0,
                if (!is.null(gsea_go)) nrow(as.data.frame(gsea_go)) else 0,
                if (!is.null(gsea_hallmark)) nrow(as.data.frame(gsea_hallmark)) else 0))
message("")
message(">>> Outputs:")
message("     ", TAB_DIR, "/ (tables + audit)")
message("     ", FIG_DIR, "/ (figures)")

writeLines(capture.output(sessionInfo()), file.path(TAB_DIR, "session_info_08.txt"))
