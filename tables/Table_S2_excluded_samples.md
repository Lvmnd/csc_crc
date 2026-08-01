## Table S2 — Excluded Samples with Reasons

### Excluded TCGA Samples (n = 4)

All exclusions were pre-registered in `03_qc_discovery.R` with automatic application — no ad-hoc exclusions.

| Sample ID | Cohort | Step | Reason | Rule |
|---|---|---|---|---|
| TCGA-A6-5656-01B-02R-A277-07 | TCGA-COAD | PCA outlier | > 4 SD from centroid on PC1 or PC2 (VST-transformed) | \|PC\| ≤ 4 SD from mean |
| TCGA-A6-5659-01B-04R-A277-07 | TCGA-COAD | PCA outlier | > 4 SD from centroid on PC1 or PC2 (VST-transformed) | \|PC\| ≤ 4 SD from mean |
| TCGA-A6-2684-01C-08R-A277-07 | TCGA-COAD | PCA outlier | > 4 SD from centroid on PC1 or PC2 (VST-transformed) | \|PC\| ≤ 4 SD from mean |
| TCGA-A6-3810-01B-04R-A277-07 | TCGA-COAD | PCA outlier | > 4 SD from centroid on PC1 or PC2 (VST-transformed) | \|PC\| ≤ 4 SD from mean |

No samples were excluded at the library size or gene filtering steps.

### Excluded GEO Samples (n = 2 from GSE107422)

| Sample ID | Cohort | Step | Reason | Rule |
|---|---|---|---|---|
| GSM2866736 | GSE107422 | PCA outlier | > 4 SD from centroid on PC1 or PC2 | \|PC\| ≤ 4 SD from mean |
| GSM2866764 | GSE107422 | PCA outlier | > 4 SD from centroid on PC1 or PC2 | \|PC\| ≤ 4 SD from mean |

(GSE220148 removed — expression-only reference cohort excluded from study)

**Total: 6 samples excluded across all cohorts (4 TCGA + 2 GEO).**

*Source files: `data/processed/tcga_qc/tcga_excluded_samples.csv`, `data/processed/geo_qc/geo_excluded_samples.csv`*
