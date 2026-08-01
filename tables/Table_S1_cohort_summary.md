## Table S1 — Dataset and Cohort Summary

### Discovery Cohort

| Characteristic | Discovery (TCGA-COAD/READ) |
|---|---|
| **Total samples** | 698 (raw) → 694 (post-QC) |
| **Tumor samples** | 647 (TP) |
| **Normal samples** | 51 (NT) |
| **COAD** | 522 |
| **READ** | 176 |
| **Patients with clinical data** | COAD: 461, READ: 172 |
| **Patients with survival data** | COAD: 102, READ: 27 |
| **Patients with risk scores** | 609 |
| **High-risk group** | 304 |
| **Low-risk group** | 305 |
| **Platform** | Illumina HiSeq (STAR-Counts) |
| **Data type** | RNA-seq raw counts |
| **Genes pre-filter** | 60,660 |
| **Genes post-filter** | 23,429 |
| **Endpoint** | Overall survival |
| **Population** | Multi-ethnic (USA-based, predominantly European ancestry) |
| **Sex distribution** | 325 M / ... (varies by available clinical data) |
| **Median age (range)** | — |
| **TNM stages** | I–IV |
| **Data access** | TCGA dbGaP phs000178 via GDC |
| **Reference** | https://portal.gdc.cancer.gov/ |

### Validation Cohorts

| Characteristic | GSE107422 | GSE71187 | GSE64857 | GSE92921 | Zenodo #8333650 |
|---|---|---|---|---|---|---|
| **Population** | Korean (Asan Medical Center, Seoul) | Chinese (Shanghai) | Chinese (Fudan University, Shanghai) | Japanese | Korean (Asan: historical cohort) |
| **Total samples** | 112 (raw) → 110 (post-QC) | 52 | 75 | 59 | 176 |
| **Tumor samples** | 110 | 52 | 75 | 59 | 176 |
| **Normal samples** | 0 | 0 | 0 | 0 | 0 |
| **Survival endpoint** | Recurrence (binary) | OS (binary) | Recurrence (binary) | DFS (time-to-event) | Recurrence (binary) |
| **Survival events** | 38 | 22 | 26 recurred | 6 | 48 |
| **Patients with risk scores** | 110 | 52 | 75 | 59 | 176 |
| **High-risk group** | 58 | — | — | — | — |
| **Low-risk group** | 52 | — | — | — | — |
| **Matched signature genes** | 187 (full) | 139 / 187 | 147 / 187 | 148 / 187 | 117 / 187 |
| **Platform** | Illumina HiSeq 2000 | Affymetrix HG-U133 Plus 2.0 (GPL570) | Affymetrix HG-U133 Plus 2.0 (GPL570) | Affymetrix HG-U133 Plus 2.0 (GPL570) | Affymetrix (custom array) |
| **Data type** | Expression matrix (log2) | Expression matrix (log2) | Expression matrix (log2) | Expression matrix (log2) | Expression matrix (log2) |
| **Endpoint** | Recurrence (binary) | OS (binary status) | Recurrence (binary) | DFS (time-to-event) | Recurrence (binary) |
| **TNM stages** | I–IV | I–IV | I–IV | I–IV | I–IV |
| **Data access** | GEO GSE107422 | GEO GSE71187 | GEO GSE64857 | GEO GSE92921 | Zenodo 10.5281/zenodo.8333650 |

*Note: Detailed demographic data (age, sex, stage) were extracted from TCGA clinical annotations. Some cells intentionally left blank where data are unavailable or not published by the original study. See the numbers ledger (`results/tables/numbers_ledger.csv`) for exact counts and source script line numbers.*
