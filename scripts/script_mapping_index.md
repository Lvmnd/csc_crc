# Script-to-Output Mapping Index

## Overview

This index maps each pipeline script to its generated tables, figures, and intermediate data files. All 17 scripts are included in `supplementary/scripts/` and zipped as `pipeline_scripts.zip`.

---

## Core Pipeline Scripts

### 00_install_packages.R
- **Purpose:** Install all required R/Bioconductor packages
- **Outputs:** None (side-effect: package installation)
- **Input dependencies:** conda environment (R 4.4.3)

### 00_set_seed.R
- **Purpose:** Set global random seed to 42 for reproducibility
- **Outputs:** None (sourced by all other scripts)
- **Tables/Figures:** None

### 01_download_tcga.R
- **Purpose:** Download TCGA-COAD and TCGA-READ STAR-Counts RNA-seq + clinical data
- **Output tables:** None (raw data download)
- **Figures:** None
- **Key outputs:** `data/raw/tcga/tcga_expression_SE.rds`, `data/raw/tcga/clinical_COAD.rds`, `data/raw/tcga/clinical_READ.rds`

### 02_download_geo_validation.R
- **Purpose:** Download GSE107422 validation datasets
- **Output tables:** None (raw data download)
- **Figures:** None
- **Key outputs:** `data/raw/geo/GSE107422/`

### 02b_parse_geo_supplementary.R
- **Purpose:** Parse supplementary files from GEO datasets
- **Output tables:** None
- **Figures:** None
- **Key outputs:** GSE107422 expression matrix

### 02c_align_sample_names.R
- **Purpose:** Harmonize sample naming conventions across GEO datasets
- **Output tables:** None
- **Figures:** None

### 03_qc_discovery.R
- **Purpose:** QC filtering for TCGA cohort (library size, gene abundance, PCA outliers, clinical missingness)
- **Output tables:**
  - `results/tables/Table_S2_excluded_samples.csv` (via tcga_excluded_samples.csv)
- **Figures:** Supplementary Figure S1–S6
  - `results/figures/qc_tcga/01_library_size_boxplot.png`
  - `results/figures/qc_tcga/02_pca_pre_filter.png`
  - `results/figures/qc_tcga/03_pca_post_filter.png`
  - `results/figures/qc_tcga/04_sample_correlation_heatmap.png`
  - `results/figures/qc_tcga/05_gene_filtering_diagnostic.png`
  - `results/figures/qc_tcga/06_missing_clinical_data.png`
- **Key outputs:** `data/processed/tcga_qc/tcga_dds.rds`, `data/processed/tcga_qc/tcga_vsd.rds`

### 03_qc_validation.R
- **Purpose:** QC filtering for GEO validation cohorts
- **Figures:** Supplementary Figure S7–S8
  - `results/figures/qc_geo/GSE107422_01_library_size.png`
  - `results/figures/qc_geo/GSE107422_02_pca.png`
  - `results/figures/qc_geo/GSE107422_03_heatmap.png`
  # (GSE220148 QC figures removed — cohort excluded from study)
- **Key outputs:** `data/processed/geo_qc/geo_excluded_samples.csv`

### 04_deg_analysis.R
- **Purpose:** Differential expression analysis (tumor vs normal) with DESeq2
- **Output tables:**
  - **Supplementary Table S3** — `results/tables/DEG_full_TCGA.csv` (23,429 genes × 10 columns)
  - `results/tables/DEG_sig_TCGA.csv` (significant DEGs only)
  - `results/tables/DEG_summary_TCGA.txt`
- **Figures:** Supplementary Figure S9
  - `results/figures/deg_tcga/01_volcano_plot.png`

### 05_csc_marker_focus.R
- **Purpose:** CSC marker expression analysis from DEG results
- **Output tables:**
  - `results/tables/CSC_markers_TCGA.csv`
  - `results/tables/CSC_markers_summary_TCGA.txt`
- **Figures:** None (but data feeds Table S4)

### 06_signature_construction.R
- **Purpose:** WGCNA + LASSO-Cox signature construction + survival analysis
- **Output tables:**
  - **Supplementary Table S4** — `results/tables/Table_risk_score_formula.csv` (187 genes + coefficients)
  - `results/tables/Cox_multivariable_summary.txt`
  - `results/tables/PH_assumption_test.txt`
- **Figures:** Supplementary Figure S10–S14
  - `results/figures/sig_tcga/01_wgcna_dendrogram.png`
  - `results/figures/sig_tcga/03_lasso_cv_curve.png`
  - `results/figures/sig_tcga/04_km_curve.png`
  - `results/figures/sig_tcga/05_forest_plot.png`
  - `results/figures/sig_tcga/06_cox_zph_plot.png`

### 07_validation.R
- **Purpose:** Independent validation in GSE107422 (logistic regression)
- **Output tables:**
  - **Supplementary Table S5** — validation results detail
  - `results/tables/Validation_outcome.txt`
  - `results/tables/Validation_logistic_summary.txt`
  - `results/tables/Validation_tertile_sensitivity.txt`
- **Figures:** Supplementary Figure S15–S16
  - `results/figures/validation/03_forest_plot.png`
  - `results/figures/validation/05_risk_score_distribution.png`

### 08_enrichment.R
- **Purpose:** Functional enrichment (GSEA/ORA) for risk-group DEGs
- **Output tables:**
  - `results/tables/enrichment/TCGA_DEG_high_vs_low.csv`
  - `results/tables/enrichment/TCGA_stemness_audit.txt`
- **Figures:** Supplementary Figure S22
  - `results/figures/enrichment/TCGA_05_stemness_heatmap.png`

### 09_immune_infiltration.R
- **Purpose:** ssGSEA immune infiltration scoring + correlation with risk score
- **Output tables:**
  - `results/tables/immune_infiltration/TCGA_immune_risk_correlation.csv`
  - `results/tables/immune_infiltration/TCGA_immune_scores.csv`
  - `results/tables/immune_infiltration/TCGA_immune_comparison_stats.csv`
  - `results/tables/immune_infiltration/Validation_immune_scores.csv`
- **Figures:** Supplementary Figure S17–S21
  - `results/figures/immune_infiltration/TCGA_01_immune_boxplots.png`
  - `results/figures/immune_infiltration/TCGA_02_immune_heatmap.png`
  - `results/figures/immune_infiltration/TCGA_03_immune_correlation_heatmap.png`
  - `results/figures/immune_infiltration/TCGA_04_risk_vs_immune_scatter.png`
  - `results/figures/immune_infiltration/Validation_01_immune_boxplots.png`

### 10_robustness_check.R
- **Purpose:** Self-audit, numbers ledger, data leakage check
- **Output tables:**
  - `results/tables/numbers_ledger.csv` (291 entries — every manuscript statistic)
  - `results/tables/robustness_audit.txt`
- **Figures:** None

### run_pipeline.ps1
- **Purpose:** PowerShell driver script to execute the entire pipeline sequentially
- **Outputs:** Console log of execution times and exit codes

---

## Supplementary Table Reference

| Table | File | Source Script | Description |
|---|---|---|---|
| Table S1 | `supplementary/tables/Table_S1_cohort_summary.md` | Compiled manually | Dataset and cohort characteristics |
| Table S2 | `supplementary/tables/Table_S2_excluded_samples.md` | 03_qc_discovery.R, 03_qc_validation.R | Samples excluded during QC with reasons |
| Table S3 | `results/tables/DEG_full_TCGA.csv` | 04_deg_analysis.R | Full differential expression results (23,429 genes) |
| Table S4 | `results/tables/Table_risk_score_formula.csv` | 06_signature_construction.R | 187 signature genes with LASSO coefficients |
| Table S5 | `supplementary/tables/Table_S5_validation_results.md` | 07_validation.R | Validation cohort results with sensitivity analysis |

## Supplementary Figure Reference

| Figure | Source File(s) | Source Script | Description |
|---|---|---|---|
| S1 | qc_tcga/01_library_size_boxplot.png | 03_qc_discovery.R | TCGA library size distribution |
| S2 | qc_tcga/02_pca_pre_filter.png | 03_qc_discovery.R | TCGA PCA before filtering |
| S3 | qc_tcga/03_pca_post_filter.png | 03_qc_discovery.R | TCGA PCA after filtering |
| S4 | qc_tcga/04_sample_correlation_heatmap.png | 03_qc_discovery.R | TCGA sample correlation heatmap |
| S5 | qc_tcga/05_gene_filtering_diagnostic.png | 03_qc_discovery.R | TCGA gene filtering diagnostics |
| S6 | qc_tcga/06_missing_clinical_data.png | 03_qc_discovery.R | TCGA missing clinical data |
| S7 | qc_geo/GSE107422_*.png | 03_qc_validation.R | GSE107422 QC (library size, PCA, heatmap) |
| S8 | (removed) | — | GSE220148 QC (removed — cohort excluded from study) |
| S9 | deg_tcga/01_volcano_plot.png | 04_deg_analysis.R | DEG volcano plot (tumor vs normal) |
| S10 | sig_tcga/01_wgcna_dendrogram.png | 06_signature_construction.R | WGCNA gene dendrogram and module colors |
| S11 | sig_tcga/03_lasso_cv_curve.png | 06_signature_construction.R | LASSO cross-validation deviance curve |
| S12 | sig_tcga/04_km_curve.png | 06_signature_construction.R | Kaplan-Meier survival curves (high vs low risk) |
| S13 | sig_tcga/05_forest_plot.png | 06_signature_construction.R | Forest plot — multivariable Cox regression |
| S14 | sig_tcga/06_cox_zph_plot.png | 06_signature_construction.R | Cox proportional hazards diagnostics (Schoenfeld residuals) |
| S15 | validation/03_forest_plot.png | 07_validation.R | Validation forest plot — logistic regression OR |
| S16 | validation/05_risk_score_distribution.png | 07_validation.R | Validation risk score distribution |
| S17 | immune_infiltration/TCGA_01_immune_boxplots.png | 09_immune_infiltration.R | TCGA immune cell scores by risk group |
| S18 | immune_infiltration/TCGA_02_immune_heatmap.png | 09_immune_infiltration.R | TCGA immune infiltration heatmap |
| S19 | immune_infiltration/TCGA_03_immune_correlation_heatmap.png | 09_immune_infiltration.R | Immune cell correlation heatmap |
| S20 | immune_infiltration/TCGA_04_risk_vs_immune_scatter.png | 09_immune_infiltration.R | Risk score vs immune score scatter plots |
| S21 | immune_infiltration/Validation_01_immune_boxplots.png | 09_immune_infiltration.R | GSE107422 immune cell scores by risk group |
| S22 | enrichment/TCGA_05_stemness_heatmap.png | 08_enrichment.R | Stemness-related pathway expression heatmap |
