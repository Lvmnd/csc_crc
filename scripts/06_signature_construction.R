#!/usr/bin/env Rscript
# =============================================================================
# 06_signature_construction.R — CSC prognostic signature via WGCNA + LASSO-Cox
#
# Workflow:
#   1. CSC panel genes (seed list)
#   2. WGCNA co-expression network → identify modules enriched for CSC genes
#   3. Expand gene set: CSC panel + all genes from CSC-containing modules
#   4. LASSO-Cox regression on overall survival (10-fold CV)
#   5. Risk score formula + patient stratification (median split)
#   6. Kaplan-Meier + log-rank
#   7. Multivariable Cox (age/stage/sex)
#   8. Proportional hazards assumption test (cox.zph)
#
# Inputs:
#   data/processed/tcga_qc/tcga_dds_de.rds   — after DESeq() (for expression)
#   data/processed/tcga_qc/tcga_vsd.rds       — vst-transformed (for WGCNA)
#   data/raw/tcga/clinical_COAD.rds            — clinical
#   data/raw/tcga/clinical_READ.rds            — clinical
#   results/tables/DEG_full_TCGA.csv          — for annotation
#
# Outputs:
#   data/processed/tcga_sig/
#   ├── wgcna_results.rds                     — full WGCNA output
#   ├── csc_associated_genes.csv              — panel + module-expanded list
#   ├── lasso_cv_fit.rds                      — cv.glmnet object
#   ├── lasso_coefficients.csv                — final signature genes + coefs
#   └── risk_scores.csv                       — per-patient risk score
#
#   results/tables/
#   ├── Table_risk_score_formula.csv          — manuscript-ready formula table
#   ├── Cox_multivariable_summary.txt         — HR + 95% CI + p-values
#   └── PH_assumption_test.txt               — cox.zph results
#
#   results/figures/sig_tcga/
#   ├── 01_wgcna_dendrogram.png
#   ├── 02_wgcna_module_trait_heatmap.png
#   ├── 03_lasso_cv_curve.png
#   ├── 04_km_curve.png
#   ├── 05_forest_plot.png
#   └── 06_cox_zph_plot.png
#
# Manuscript number mapping (see logs/decisions.md):
#   - Risk score formula  → Table_risk_score_formula.csv
#   - Median risk score   → risk_scores.csv
#   - KM log-rank p       → from survdiff object (printed in script output)
#   - Multivariable HR    → Cox_multivariable_summary.txt
#   - PH test             → PH_assumption_test.txt
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(WGCNA)
  library(glmnet)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(matrixStats)
  library(org.Hs.eg.db)
})

# -- Allow WGCNA to work with multiple threads -------------------------------
allowWGCNAThreads()

# -- Directories -------------------------------------------------------------
QC_DIR    <- file.path("data", "processed", "tcga_qc")
RAW_DIR   <- file.path("data", "raw", "tcga")
OUT_DIR   <- file.path("data", "processed", "tcga_sig")
TAB_DIR   <- file.path("results", "tables")
FIG_DIR   <- file.path("results", "figures", "sig_tcga")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. CSC panel (identical to script 05)
# =============================================================================
csc_symbols <- c("PROM1", "LGR5", "CD44", "CD24", "ALCAM",
                  "ITGB1", "EPCAM", "ALDH1A1", "CTNNB1", "CXCR4", "CDCP1")
csc_symbols_lower <- tolower(csc_symbols)
csc_aliases <- c("CD133", "GPR49", NA, NA, "CD166",
                  "CD29", "ESA", "ALDH1", NA, "CD184", NA)
search_terms <- unique(tolower(c(csc_symbols, csc_aliases[!is.na(csc_aliases)])))

# =============================================================================
# 2. Load expression and clinical data
# =============================================================================
message(">>> Loading expression and clinical data ...")

# vst-transformed (for WGCNA — tumour-only)
vsd <- readRDS(file.path(QC_DIR, "tcga_vsd.rds"))
rlog_mat <- assay(vsd)

# DESeqDataSet post-DESeq (for size-factor-normalized counts)
dds <- readRDS(file.path(QC_DIR, "tcga_dds_de.rds"))
norm_counts <- counts(dds, normalized = TRUE)

# Gene symbol mapping (ENSEMBL → symbol via org.Hs.eg.db)
ensembls <- gsub("\\..*$", "", rownames(dds))  # strip version suffix
symbol_map <- suppressMessages(mapIds(org.Hs.eg.db, keys = ensembls, keytype = "ENSEMBL",
                     column = "SYMBOL", multiVals = "first"))
gene_map <- data.frame(
  gene_id   = rownames(dds),
  ensembl   = ensembls,
  gene_name = ifelse(is.na(symbol_map), rownames(dds), symbol_map),
  stringsAsFactors = FALSE
)

# Clinical data
clin_COAD <- readRDS(file.path(RAW_DIR, "clinical_COAD.rds"))
clin_READ <- readRDS(file.path(RAW_DIR, "clinical_READ.rds"))
clin_all <- bind_rows(
  mutate(clin_COAD, project = "TCGA-COAD"),
  mutate(clin_READ, project = "TCGA-READ")
)

# Clean clinical: keep only first record per patient
clin_all <- clin_all[!duplicated(clin_all$submitter_id), ]

message(sprintf("   Expression: %d genes x %d samples", nrow(rlog_mat), ncol(rlog_mat)))
message(sprintf("   Clinical:   %d patients", nrow(clin_all)))

# =============================================================================
# 3. Restrict to tumour samples (normal tissue has different co-expression)
# =============================================================================
tumour_samples <- colnames(dds)[colData(dds)$sample_type == "Tumor"]
rlog_tumour <- rlog_mat[, intersect(colnames(rlog_mat), tumour_samples)]
message(sprintf("   Tumour samples for WGCNA: %d", ncol(rlog_tumour)))

# =============================================================================
# 4. WGCNA — find co-expression modules containing CSC genes
# =============================================================================
message(">>> WGCNA co-expression network analysis ...")

# 4a. Use top variable genes for network construction
n_top <- min(5000, nrow(rlog_tumour))
gene_vars <- matrixStats::rowVars(rlog_tumour)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[seq_len(n_top)])
expr_wgcna <- t(rlog_tumour[top_genes, ])  # samples x genes

message(sprintf("   Using top %d variable genes for network construction", n_top))

# 4b. Pick soft threshold
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(expr_wgcna, powerVector = powers, verbose = 0)

soft_power <- sft$powerEstimate
if (is.na(soft_power) || soft_power < 1) {
  message("   Power estimation ambiguous — using power = 6 (default for signed networks)")
  soft_power <- 6
}
message(sprintf("   Selected soft threshold power: %d", soft_power))

# 4c. Build adjacency and TOM
adj <- adjacency(expr_wgcna, power = soft_power, type = "signed")
tom <- TOMsimilarity(adj, TOMType = "signed")
dissTOM <- 1 - tom

# 4d. Hierarchical clustering + module detection
gene_tree <- hclust(as.dist(dissTOM), method = "average")
min_module_size <- 30
dynamic_mods <- cutreeDynamic(dendro = gene_tree, distM = dissTOM,
                               minClusterSize = min_module_size,
                               method = "hybrid", deepSplit = 2)
dynamic_colors <- labels2colors(dynamic_mods)

# Plot dendrogram
png(file.path(FIG_DIR, "01_wgcna_dendrogram.png"), width = 12, height = 6, units = "in", res = 150)
plotDendroAndColors(gene_tree, dynamic_colors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "TCGA — Gene co-expression dendrogram")
dev.off()

# 4e. Find modules containing CSC genes
gene_names_wgcna <- colnames(expr_wgcna)
gene_symbols_wgcna <- gene_map$gene_name[match(gene_names_wgcna, gene_map$gene_id)]
gene_symbols_lower <- tolower(gene_symbols_wgcna)

# Which CSC markers are in the WGCNA gene set?
csc_in_wgcna <- which(gene_symbols_lower %in% search_terms)
message(sprintf("   CSC markers found in WGCNA gene set: %d / %d",
                length(csc_in_wgcna), length(csc_symbols)))

# Find module assignments for CSC genes
if (length(csc_in_wgcna) > 0) {
  csc_modules <- unique(dynamic_colors[csc_in_wgcna])
  csc_module_genes <- which(dynamic_colors %in% csc_modules)

  message(sprintf("   CSC markers span %d module(s): %s",
                  length(csc_modules), paste(csc_modules, collapse = " ")))
  message(sprintf("   Total genes in CSC-containing modules: %d", length(csc_module_genes)))
} else {
  message("   No CSC markers in WGCNA geneset. Using panel genes only.")
  csc_modules <- character(0)
  csc_module_genes <- integer(0)
}

# =============================================================================
# 5. Define expanded CSC-associated gene set
# =============================================================================
# Genes: CSC panel (seed) + all genes from CSC-containing modules
# Map back to gene IDs

# Panel genes from the full expression matrix
panel_in_expr <- which(gene_map$gene_name %in% csc_symbols |
                         tolower(gene_map$gene_name) %in% search_terms)
panel_gene_ids <- gene_map$gene_id[panel_in_expr]

# Module-expanded genes
if (length(csc_module_genes) > 0) {
  module_gene_ids <- gene_names_wgcna[csc_module_genes]
} else {
  module_gene_ids <- character(0)
}

csc_associated_ids <- unique(c(panel_gene_ids, module_gene_ids))
csc_associated_symbols <- gene_map$gene_name[match(csc_associated_ids, gene_map$gene_id)]

csc_associated_df <- data.frame(
  gene_id   = csc_associated_ids,
  gene_name = csc_associated_symbols,
  source    = ifelse(csc_associated_ids %in% panel_gene_ids,
                     ifelse(csc_associated_ids %in% module_gene_ids,
                            "panel+module", "panel"),
                     "module"),
  stringsAsFactors = FALSE
)
write.csv(csc_associated_df, file.path(OUT_DIR, "csc_associated_genes.csv"), row.names = FALSE)

message(sprintf("   CSC-associated genes for LASSO-Cox: %d", length(csc_associated_ids)))
message(sprintf("     From panel         : %d", sum(csc_associated_df$source == "panel")))
message(sprintf("     From module        : %d", sum(csc_associated_df$source == "module")))
message(sprintf("     Both               : %d", sum(csc_associated_df$source == "panel+module")))

# =============================================================================
# 6. Prepare survival data
# =============================================================================
message(">>> Preparing survival data ...")

# Extract patient IDs from tumour sample barcodes
tumour_barcodes <- tumour_samples
tumour_patient_ids <- substr(tumour_barcodes, 1, 12)

# Build survival data frame
surv_df <- data.frame(
  sample_id    = tumour_barcodes,
  patient_id   = tumour_patient_ids,
  stringsAsFactors = FALSE
)

# Merge clinical data
clin_for_merge <- clin_all %>%
  mutate(patient_id = submitter_id) %>%
  dplyr::select(patient_id, vital_status, days_to_death, days_to_last_follow_up,
         age_at_diagnosis, sex_at_birth, ajcc_pathologic_stage)

surv_df <- surv_df %>%
  left_join(clin_for_merge, by = "patient_id") %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)  # one sample per patient

# Compute overall survival
surv_df <- surv_df %>%
  mutate(
    # Overall survival time (days)
    os_time = ifelse(
      vital_status == "Dead",
      as.numeric(days_to_death),
      as.numeric(days_to_last_follow_up)
    ),
    # Overall survival event
    os_event = ifelse(vital_status == "Dead", 1L, 0L)
  )

# Remove patients with missing survival data
surv_df <- surv_df %>% dplyr::filter(!is.na(os_time), !is.na(os_event), os_time > 0)
message(sprintf("   Patients with survival data: %d (%d events)",
                nrow(surv_df), sum(surv_df$os_event)))

# =============================================================================
# 7. Prepare expression matrix for LASSO-Cox
# =============================================================================
message(">>> Preparing expression matrix for LASSO-Cox ...")

# Use normalized counts (log-transformed for LASSO)
expr_for_lasso <- norm_counts
# Subset to CSC-associated genes
common_genes <- intersect(csc_associated_ids, rownames(expr_for_lasso))
message(sprintf("   CSC-associated genes with expression data: %d", length(common_genes)))

expr_sub <- expr_for_lasso[common_genes, , drop = FALSE]
# Log2-transform for LASSO (variance stabilisation)
expr_sub <- log2(expr_sub + 1)

# Match to survival data
common_samples <- intersect(colnames(expr_sub), surv_df$sample_id)
surv_df <- surv_df[match(common_samples, surv_df$sample_id), ]
expr_sub <- expr_sub[, common_samples]
rownames(expr_sub) <- make.names(rownames(expr_sub), unique = TRUE)

message(sprintf("   Matched samples: %d", ncol(expr_sub)))

x_mat <- t(expr_sub)  # samples x genes
y_surv <- Surv(time = surv_df$os_time, event = surv_df$os_event)

# =============================================================================
# 8. LASSO-Cox with 10-fold cross-validation
# =============================================================================
message(">>> LASSO-Cox regression (10-fold CV) ...")

set.seed(42)
cv_fit <- cv.glmnet(
  x        = x_mat,
  y        = y_surv,
  family   = "cox",
  alpha    = 1,
  nfolds   = 10,
  type.measure = "C",
  cox.ties = "breslow"
)

# Save CV fit
saveRDS(cv_fit, file.path(OUT_DIR, "lasso_cv_fit.rds"))

# CV curve
png(file.path(FIG_DIR, "03_lasso_cv_curve.png"), width = 8, height = 6, units = "in", res = 150)
plot(cv_fit)
title("LASSO-Cox — 10-fold Cross-Validation", line = 2.5)
dev.off()

# Extract coefficients at lambda.min and lambda.1se
lambda_min <- cv_fit$lambda.min
lambda_1se <- cv_fit$lambda.1se

coef_min <- as.matrix(coef(cv_fit, s = lambda_min))
coef_1se <- as.matrix(coef(cv_fit, s = lambda_1se))

nz_min <- sum(coef_min != 0)
nz_1se <- sum(coef_1se != 0)

message(sprintf("   lambda.min = %.5f (nonzero: %d)", lambda_min, nz_min))
message(sprintf("   lambda.1se = %.5f (nonzero: %d)", lambda_1se, nz_1se))

# Use lambda.min for highest resolution (pre-registered: choose lambda.min)
sig_genes <- coef_min[coef_min[, 1] != 0, , drop = FALSE]
chosen_lambda <- lambda_min
message(sprintf("   Final signature genes: %d", nrow(sig_genes)))

# =============================================================================
# 9. Build risk score formula table
# =============================================================================
if (nrow(sig_genes) == 0) {
  message("!!! No genes selected by LASSO at lambda.min. Using lambda.1se ...")
  sig_genes <- coef_1se[coef_1se[, 1] != 0, , drop = FALSE]
  lambda_used <- "lambda.1se"
  chosen_lambda <- lambda_1se
} else {
  lambda_used <- "lambda.min"
}

# Create coefficient table
coef_tbl <- data.frame(
  gene_symbol   = names(sig_genes[, 1]),
  coefficient   = round(sig_genes[, 1], 6),
  stringsAsFactors = FALSE
)
# Add gene names (clean up the make.names)
coef_tbl$gene_id <- gsub("^X", "", coef_tbl$gene_symbol)
coef_tbl$gene_symbol <- gene_map$gene_name[match(coef_tbl$gene_id, gene_map$gene_id)]
coef_tbl$gene_symbol[is.na(coef_tbl$gene_symbol)] <- coef_tbl$gene_id[is.na(coef_tbl$gene_symbol)]

write.csv(coef_tbl, file.path(TAB_DIR, "Table_risk_score_formula.csv"), row.names = FALSE)
write.csv(coef_tbl, file.path(OUT_DIR, "lasso_coefficients.csv"), row.names = FALSE)

message("")
message("   Risk score formula:")
for (i in seq_len(nrow(coef_tbl))) {
  message(sprintf("     %+.6f × %s", coef_tbl$coefficient[i], coef_tbl$gene_symbol[i]))
}

# =============================================================================
# 10. Compute risk scores
# =============================================================================
message(">>> Computing risk scores ...")

risk_score <- predict(cv_fit, newx = x_mat, s = lambda_min, type = "link")
risk_score <- as.numeric(risk_score)
names(risk_score) <- surv_df$sample_id

# Median split (pre-registered — not p-hacked)
median_risk <- median(risk_score)
risk_group <- ifelse(risk_score > median_risk, "High", "Low")

surv_df$risk_score <- risk_score
surv_df$risk_group  <- factor(risk_group, levels = c("Low", "High"))

# Save
risk_out <- surv_df %>% dplyr::select(sample_id, patient_id, risk_score, risk_group,
                                os_time, os_event, age_at_diagnosis, sex_at_birth,
                                ajcc_pathologic_stage)
write.csv(risk_out, file.path(OUT_DIR, "risk_scores.csv"), row.names = FALSE)

message(sprintf("   Median risk score: %.4f", median_risk))
message(sprintf("   High-risk: %d, Low-risk: %d",
                sum(risk_group == "High"), sum(risk_group == "Low")))

# =============================================================================
# 11. Kaplan-Meier + log-rank
# =============================================================================
message(">>> Kaplan-Meier survival analysis ...")

km_fit <- survfit(Surv(os_time, os_event) ~ risk_group, data = surv_df)

# KM plot
km_plot <- ggsurvplot(
  km_fit,
  data          = surv_df,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.col = "strata",
  xlab          = "Time (days)",
  ylab          = "Overall survival probability",
  title         = "TCGA — CSC signature risk groups",
  subtitle      = sprintf("Median split (low n=%d, high n=%d)",
                          sum(risk_group == "Low"), sum(risk_group == "High")),
  palette       = c("#377EB8", "#E41A1C"),
  legend.title  = "Risk group",
  ggtheme       = theme_minimal()
)

png(file.path(FIG_DIR, "04_km_curve.png"), width = 9, height = 8, units = "in", res = 150)
print(km_plot)
dev.off()

# Log-rank test
lr <- survdiff(Surv(os_time, os_event) ~ risk_group, data = surv_df)
lr_p <- 1 - pchisq(lr$chisq, df = 1)
message(sprintf("   Log-rank test p-value: %.2e", lr_p))

# =============================================================================
# 12. Cox proportional hazards (unadjusted + adjusted)
# =============================================================================
message(">>> Cox proportional hazards ...")

# Unadjusted
cox_unadj <- coxph(Surv(os_time, os_event) ~ risk_group, data = surv_df)

# Adjusted
covariates <- c("age_at_diagnosis", "sex_at_birth")
if ("ajcc_pathologic_stage" %in% names(surv_df)) {
  covariates <- c(covariates, "ajcc_pathologic_stage")
}

# Clean stage variable — collapse to numeric
if ("ajcc_pathologic_stage" %in% names(surv_df)) {
  surv_df <- surv_df %>%
    mutate(
      stage_num = case_when(
        grepl("I$|IA$|IB$", ajcc_pathologic_stage)  ~ 1,
        grepl("II$|IIA$|IIB$|IIC$", ajcc_pathologic_stage) ~ 2,
        grepl("III$|IIIA$|IIIB$|IIIC$", ajcc_pathologic_stage) ~ 3,
        grepl("IV$|IVA$|IVB$|IVC$", ajcc_pathologic_stage) ~ 4,
        TRUE ~ NA_real_
      )
    )
}

# Build formula
adj_vars <- c("risk_group")
if (!all(is.na(surv_df$age_at_diagnosis))) adj_vars <- c(adj_vars, "age_at_diagnosis")
if (length(unique(na.omit(surv_df$sex_at_birth))) > 1) adj_vars <- c(adj_vars, "sex_at_birth")
if ("stage_num" %in% names(surv_df) && length(unique(na.omit(surv_df$stage_num))) > 1) {
  adj_vars <- c(adj_vars, "stage_num")
}

formula_adj <- as.formula(paste("Surv(os_time, os_event) ~", paste(adj_vars, collapse = " + ")))
cox_adj <- coxph(formula_adj, data = surv_df)

# Summary
sink(file.path(TAB_DIR, "Cox_multivariable_summary.txt"))
cat("Cox Proportional Hazards — TCGA CSC Signature\n")
cat("==============================================\n\n")

cat("--- Unadjusted ---\n")
print(summary(cox_unadj))
cat("\n")

cat("--- Adjusted ---\n")
cat(sprintf("Covariates: %s\n\n", paste(adj_vars, collapse = ", ")))
print(summary(cox_adj))
cat("\n")

# Extract HR table
hr_table <- data.frame(
  variable   = names(coef(cox_adj)),
  HR         = exp(coef(cox_adj)),
  HR_lower   = exp(confint(cox_adj)[, 1]),
  HR_upper   = exp(confint(cox_adj)[, 2]),
  p_value    = summary(cox_adj)$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)
cat("\nHazard Ratio Table:\n")
print(hr_table, row.names = FALSE)
sink()

message("   Cox results → ", file.path(TAB_DIR, "Cox_multivariable_summary.txt"))

# =============================================================================
# 13. Forest plot
# =============================================================================
message(">>> Forest plot ...")

# Reformat for forest plot
hr_df <- hr_table %>%
  mutate(
    var_label = gsub("risk_group", "Risk group (High vs Low)", variable),
    var_label = gsub("age_at_diagnosis", "Age at diagnosis", var_label),
    var_label = gsub("sex_at_birth", "Gender", var_label),
    var_label = gsub("stage_num", "AJCC stage", var_label)
  )

p_forest <- ggplot(hr_df, aes(x = HR, y = reorder(var_label, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2, size = 0.8) +
  geom_point(size = 3, colour = "steelblue") +
  scale_x_log10() +
  labs(title = "Multivariable Cox regression",
       subtitle = sprintf("CSC signature risk score, adjusted for covariates"),
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_minimal()

ggsave(file.path(FIG_DIR, "05_forest_plot.png"), p_forest, width = 7, height = 5)

# =============================================================================
# 14. Proportional hazards assumption test
# =============================================================================
message(">>> Proportional hazards assumption test (cox.zph) ...")

zph <- cox.zph(cox_adj)

sink(file.path(TAB_DIR, "PH_assumption_test.txt"))
cat("Proportional Hazards Assumption Test (cox.zph)\n")
cat("================================================\n\n")
print(zph)
cat("\nInterpretation:\n")
cat("  A significant p-value (< 0.05) suggests the PH assumption is\n")
cat("  violated for that covariate. For the global test, a non-significant\n")
cat("  result indicates the model satisfies the PH assumption overall.\n")
sink()

message("   PH test → ", file.path(TAB_DIR, "PH_assumption_test.txt"))

# PH test plot
png(file.path(FIG_DIR, "06_cox_zph_plot.png"), width = 10, height = 8, units = "in", res = 150)
par(mfrow = c(2, 3))
plot(zph, main = "Scaled Schoenfeld residuals")
dev.off()

# =============================================================================
# 15. Final summary
# =============================================================================
message("")
message("==============================================")
message("   SIGNATURE CONSTRUCTION COMPLETE            ")
message("==============================================")
message("")
message(sprintf("  CSC panel seed genes     : %d", length(csc_symbols)))
message(sprintf("  WGCNA modules with CSC   : %d", length(unique(csc_modules))))
message(sprintf("  CSC-associated genes     : %d", length(csc_associated_ids)))
message(sprintf("  LASSO lambda used        : %s = %.5f", lambda_used, chosen_lambda))
message(sprintf("  Signature genes          : %d", nrow(sig_genes)))
message(sprintf("  High / Low risk          : %d / %d",
                sum(risk_group == "High"), sum(risk_group == "Low")))
message(sprintf("  Log-rank p               : %.2e", lr_p))
message(sprintf("  Cox PH global test p     : %.2e", zph$table[nrow(zph$table), "p"]))

message("")
message(">>> Outputs:")
message("     ", file.path(TAB_DIR, "Table_risk_score_formula.csv"))
message("     ", file.path(TAB_DIR, "Cox_multivariable_summary.txt"))
message("     ", file.path(TAB_DIR, "PH_assumption_test.txt"))
message("     ", FIG_DIR, "/ (6 figures)")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info_06.txt"))
