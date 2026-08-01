#!/usr/bin/env Rscript
# =============================================================================
# 05_csc_marker_focus.R — CSC marker dysregulation in TCGA
#
# Purpose:
#   Subset the full TCGA DEG table to the curated CRC CSC marker panel
#   and report honestly which markers are significantly dysregulated,
#   with actual fold-changes and adjusted p-values. Do NOT assume all
#   panel genes will be significant.
#
# Inputs:
#   results/tables/DEG_full_TCGA.csv    — full DESeq2 results
#
# Outputs:
#   results/tables/
#   ├── CSC_markers_TCGA.csv            — full panel with DE status
#   └── CSC_markers_summary_TCGA.txt    — narrative for manuscript
#
# CSC Panel (surface/functional markers with published evidence in CRC):
#   PROM1 (CD133), LGR5, CD44, CD24, ALCAM (CD166), ITGB1 (CD29),
#   EPCAM, ALDH1A1, CTNNB1 (β-catenin), CXCR4, CDCP1
#
# Source: PMC3645378, PMC3981033, PMC4657969, PMC10383310, PMC3506563
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(org.Hs.eg.db)
})

# Null-coalescing operator (used for alias display)
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# -- Thresholds (identical to script 04) -------------------------------------
LOG2FC_CUTOFF <- 1
PADJ_CUTOFF   <- 0.05

# -- Directories -------------------------------------------------------------
OUT_DIR <- file.path("results", "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. CSC marker panel
# =============================================================================
# Gene symbols accepted by the panel.  Some markers have multiple aliases;
# we search both the symbol and common alias.
csc_panel <- data.frame(
  symbol        = c("PROM1", "LGR5", "CD44", "CD24", "ALCAM",
                     "ITGB1", "EPCAM", "ALDH1A1", "CTNNB1", "CXCR4", "CDCP1"),
  alias         = c("CD133", "GPR49", NA, NA, "CD166",
                     "CD29", "ESA", "ALDH1", NA, "CD184", NA),
  full_name     = c("Prominin-1 (CD133)", "Leucine-rich repeat-containing G-protein coupled receptor 5",
                     "CD44 antigen", "CD24 antigen", "Activated leukocyte cell adhesion molecule (CD166)",
                     "Integrin beta-1 (CD29)", "Epithelial cell adhesion molecule",
                     "Aldehyde dehydrogenase 1 family member A1",
                     "Catenin beta-1 (β-catenin)", "C-X-C chemokine receptor type 4 (CD184)",
                     "CUB domain-containing protein 1"),
  stringsAsFactors = FALSE
)

# =============================================================================
# 2. Load DEG results
# =============================================================================
message(">>> Loading DEG results ...")
deg <- read.csv(file.path(OUT_DIR, "DEG_full_TCGA.csv"), stringsAsFactors = FALSE)

# =============================================================================
# 3. Map ENSEMBL IDs to gene symbols
# =============================================================================
message(">>> Mapping ENSEMBL IDs to gene symbols ...")
deg$ensembl <- gsub("\\..*$", "", deg$gene_name)  # strip version suffix
symbol_map <- mapIds(org.Hs.eg.db, keys = deg$ensembl, keytype = "ENSEMBL",
                     column = "SYMBOL", multiVals = "first")
deg$symbol <- symbol_map[deg$ensembl]
n_mapped <- sum(!is.na(deg$symbol))
message(sprintf("   Mapped %d / %d genes to symbols (%.1f%%)", n_mapped, nrow(deg),
                100 * n_mapped / nrow(deg)))

# =============================================================================
# 4. Match CSC markers in DEG table
# =============================================================================
message(">>> Matching CSC markers in DEG results ...")

# Build search vectors
search_terms <- unique(c(
  tolower(csc_panel$symbol),
  tolower(csc_panel$alias[!is.na(csc_panel$alias)])
))

# Match symbol column (case-insensitive)
deg$symbol_lower <- tolower(deg$symbol)
matched_idx <- which(deg$symbol_lower %in% search_terms)

if (length(matched_idx) == 0) {
  # Fallback: search in gene_name (ENSEMBL) 
  message("   No symbol matches — attempting ENSEMBL ID fallback ...")
  matched_idx <- which(
    sapply(tolower(deg$symbol), function(g) {
      any(sapply(search_terms, function(st) grepl(st, g, fixed = TRUE)))
    })
  )
}

csc_results <- deg[matched_idx, ] %>%
  mutate(
    csc_symbol   = toupper(symbol),
    significant  = ifelse(is.na(padj) | is.na(log2FoldChange), FALSE,
                          abs(log2FoldChange) > LOG2FC_CUTOFF & padj < PADJ_CUTOFF),
    direction    = case_when(
      significant & log2FoldChange > 0  ~ "Up",
      significant & log2FoldChange < 0  ~ "Down",
      TRUE                              ~ "NS"
    )
  ) %>%
  arrange(desc(significant), padj)

# Merge with panel metadata
# Join: csc_panel (left) + deg results (right) on symbol (case-insensitive)
csc_results <- csc_panel %>%
  mutate(search_lower = tolower(symbol)) %>%
  left_join(
    csc_results %>% mutate(search_lower = tolower(symbol)),
    by = "search_lower",
    suffix = c(".panel", ".deg")
  ) %>%
  dplyr::select(symbol = symbol.panel, alias, full_name,
         gene_id, gene_name,
         baseMean, log2FoldChange, lfcSE, pvalue, padj, significant, direction)

# Handle unmatched markers
csc_results$significant[is.na(csc_results$significant)] <- FALSE

# =============================================================================
# 4. Report
# =============================================================================
n_sig <- sum(csc_results$significant, na.rm = TRUE)
n_total <- nrow(csc_results)

message("")
message("==============================================")
message("   CSC MARKER DYSREGULATION IN TCGA           ")
message("==============================================")
for (i in seq_len(n_total)) {
  r <- csc_results[i, ]
  if (is.na(r$gene_id)) {
    msg <- sprintf("   %-8s (%-25s) %s", r$symbol, r$alias %||% "",
                   "⚠ NOT FOUND in expression data")
  } else if (r$significant) {
    msg <- sprintf("   %-8s (%-25s) %s  log2FC = %+.2f  padj = %.2e",
                   r$symbol, r$alias %||% "",
                   ifelse(r$direction == "Up", "↑ SIGNIFICANT", "↓ SIGNIFICANT"),
                   r$log2FoldChange, r$padj)
  } else if (!is.na(r$padj)) {
    msg <- sprintf("   %-8s (%-25s) %s  log2FC = %+.2f  padj = %.2e",
                   r$symbol, r$alias %||% "",
                   "— NOT significant",
                   r$log2FoldChange, r$padj)
  } else {
    msg <- sprintf("   %-8s (%-25s) %s", r$symbol, r$alias %||% "",
                   "— low count / filtered out")
  }
  message(msg)
}

message("")
message(sprintf("   %d / %d CSC markers significantly dysregulated", n_sig, n_total))

# =============================================================================
# 5. Save outputs
# =============================================================================
message(">>> Saving results ...")

write.csv(csc_results, file.path(OUT_DIR, "CSC_markers_TCGA.csv"), row.names = FALSE)

sink(file.path(OUT_DIR, "CSC_markers_summary_TCGA.txt"))
cat("CSC Marker Dysregulation in TCGA CRC\n")
cat("======================================\n\n")
cat(sprintf("Panel size               : %d markers\n", n_total))
cat(sprintf("Significantly dysregulated: %d (%.0f%%)\n", n_sig, 100 * n_sig / n_total))
cat(sprintf("  Up-regulated           : %d\n", sum(csc_results$direction == "Up", na.rm = TRUE)))
cat(sprintf("  Down-regulated         : %d\n", sum(csc_results$direction == "Down", na.rm = TRUE)))
cat(sprintf("  Not significant        : %d\n", sum(csc_results$direction == "NS" & !is.na(csc_results$padj), na.rm = TRUE)))
cat(sprintf("  Not detected           : %d\n", sum(is.na(csc_results$gene_id))))
cat("\nPer-marker detail:\n")
cat(sprintf("%-8s %-6s %8s %10s %8s\n", "Marker", "Status", "log2FC", "padj", "Sig?"))
cat(sprintf("%-8s %-6s %8s %10s %8s\n", "------", "------", "------", "----", "----"))
for (i in seq_len(n_total)) {
  r <- csc_results[i, ]
  sig_label <- if (isTRUE(r$significant)) "YES" else "no"
  fc <- if (!is.na(r$log2FoldChange)) sprintf("%+.2f", r$log2FoldChange) else "NA"
  p  <- if (!is.na(r$padj)) format(r$padj, scientific = TRUE, digits = 3) else "NA"
  status <- if (is.na(r$gene_id)) "missing" else if (r$direction == "Up") "UP" else if (r$direction == "Down") "DOWN" else "NS"
  cat(sprintf("%-8s %-6s %8s %10s %8s\n", r$symbol, status, fc, p, sig_label))
}
cat("\n— End of CSC marker summary —\n")
sink()

message("")
message(">>> Outputs:")
message("     ", file.path(OUT_DIR, "CSC_markers_TCGA.csv"))
message("     ", file.path(OUT_DIR, "CSC_markers_summary_TCGA.txt"))
message(">>> Proceed to scripts/06_signature_construction.R")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info_05.txt"))
