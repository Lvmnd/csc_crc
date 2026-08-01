#!/usr/bin/env Rscript
# =============================================================================
# 01_download_tcga.R — Download TCGA-COAD & TCGA-READ RNA-seq + clinical data
#
# Purpose:
#   Query and download STAR-Counts gene expression quantification for
#   TCGA colon adenocarcinoma (COAD) and rectum adenocarcinoma (READ)
#   together with matched clinical and survival data via TCGAbiolinks.
#
# Outputs:
#   data/raw/tcga/
#   ├── GDCquery_manifest.txt          # query results summary
#   ├── GDCdownload_log.txt             # download log
#   ├── clinical_COAD.csv               # clinical data (COAD)
#   ├── clinical_READ.csv               # clinical data (READ)
#   └── GDCdata/                        # downloaded expression files
#
# Cohort A (discovery): TCGA-COAD + TCGA-READ combined.
# =============================================================================

source("scripts/00_set_seed.R")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# -- Directories -------------------------------------------------------------
RAW_DIR   <- "data/raw/tcga"
MANIFEST  <- file.path(RAW_DIR, "GDCquery_manifest.txt")
LOG_FILE  <- file.path(RAW_DIR, "GDCdownload_log.txt")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

# -- 1. Query RNA-seq expression (STAR - Counts) ----------------------------
message(">>> Querying GDC for TCGA-COAD and TCGA-READ RNA-seq ...")

query_exp <- GDCquery(
  project       = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

# Save manifest
manifest_df <- getResults(query_exp)
write.table(manifest_df, file = MANIFEST, sep = "\t", quote = FALSE, row.names = FALSE)
message("   Manifest saved → ", MANIFEST)

# -- 2. Download ------------------------------------------------------------
message(">>> Downloading RNA-seq data (API method) ...")
sink(LOG_FILE, split = TRUE)
tryCatch(
  GDCdownload(
    query          = query_exp,
    method         = "api",
    directory      = RAW_DIR,
    files.per.chunk = 50
  ),
  error = function(e) {
    message("API download failed, trying GDC client ...")
    GDCdownload(
      query     = query_exp,
      method    = "client",
      directory = RAW_DIR
    )
  }
)
sink()
message("   Download log → ", LOG_FILE)

# -- 3. Prepare SummarizedExperiment ----------------------------------------
message(">>> Preparing expression data ...")
data_exp <- GDCprepare(
  query                 = query_exp,
  directory             = RAW_DIR,
  summarizedExperiment  = TRUE
)

# Raw count matrix
counts <- assay(data_exp, "unstranded")
message("   Expression matrix dims: ", nrow(counts), " genes x ", ncol(counts), " samples")

# -- 4. Download clinical data ----------------------------------------------
message(">>> Downloading clinical data ...")
clin_COAD <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
clin_READ <- GDCquery_clinic(project = "TCGA-READ", type = "clinical")
write.csv(clin_COAD, file.path(RAW_DIR, "clinical_COAD.csv"), row.names = FALSE)
write.csv(clin_READ, file.path(RAW_DIR, "clinical_READ.csv"), row.names = FALSE)

# -- 5. Sample summary -------------------------------------------------------
message("")
message("==============================================")
message("   SAMPLE SUMMARY — TCGA-COAD + TCGA-READ     ")
message("==============================================")

# Classify sample types from barcodes
sample_barcodes <- colnames(data_exp)
tumor_idx  <- TCGAquery_SampleTypes(sample_barcodes, typesample = "TP")
normal_idx <- TCGAquery_SampleTypes(sample_barcodes, typesample = "NT")

n_tumor  <- length(tumor_idx)
n_normal <- length(normal_idx)

message("   Total samples      : ", length(sample_barcodes))
message("   Primary Tumor (TP) : ", n_tumor)
message("   Solid Tissue Normal: ", n_normal)

# Per project
proj_vec <- colData(data_exp)$project_id
get_sample_type_code <- function(barcodes) {
  substr(barcodes, 14, 15)
}
is_tumor <- function(barcodes) get_sample_type_code(barcodes) == "01"
is_normal <- function(barcodes) get_sample_type_code(barcodes) == "11"

message("")
message("   By project:")
for (proj in unique(proj_vec)) {
  idx <- which(proj_vec == proj)
  bcs <- colnames(data_exp)[idx]
  t <- sum(is_tumor(bcs))
  n <- sum(is_normal(bcs))
  message(sprintf("     %-12s : %3d tumor + %d normal = %d", proj, t, n, length(idx)))
}

# Clinical summary
message("")
message("   Clinical data available:")
message(sprintf("     TCGA-COAD : %d patients (%d with survival)", 
                nrow(clin_COAD), sum(!is.na(clin_COAD$days_to_death))))
message(sprintf("     TCGA-READ : %d patients (%d with survival)", 
                nrow(clin_READ), sum(!is.na(clin_READ$days_to_death))))

# -- 6. Save prepared data for downstream -----------------------------------
message("")
message(">>> Saving prepared SummarizedExperiment and clinical data as RDS ...")

saveRDS(data_exp, file.path(RAW_DIR, "tcga_expression_SE.rds"))
saveRDS(clin_COAD, file.path(RAW_DIR, "clinical_COAD.rds"))
saveRDS(clin_READ, file.path(RAW_DIR, "clinical_READ.rds"))
saveRDS(query_exp, file.path(RAW_DIR, "GDCquery_object.rds"))

message("     tcga_expression_SE.rds    — SummarizedExperiment (counts + colData)")
message("     clinical_COAD.rds          — clinical data (COAD)")
message("     clinical_READ.rds          — clinical data (READ)")
message("     GDCquery_object.rds        — original GDCquery for reloading")

message("")
message(">>> TCGA download complete.")
message(">>> Proceed to scripts/02_download_geo_validation.R")

# -- Save session info for reproducibility ----------------------------------
writeLines(capture.output(sessionInfo()), file.path(RAW_DIR, "session_info.txt"))
