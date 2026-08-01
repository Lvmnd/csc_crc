# A Pre-Specified Cancer Stem Cell Gene Panel Shows No Independent Prognostic Value in Colorectal Cancer After Stage Adjustment: An Asian Populations Multi-Cohort Validation Analysis 

> ⚠️ **PLACEHOLDER NOTICE:** Sections marked `[CONFIRM BEFORE PUBLISHING]` contain
> values that are not yet locked/reconciled as of the last analysis run. Do not
> commit this README with those placeholders filled from memory — pull them
> directly from `results/tables/numbers_ledger.csv` only once finalized.

This repository contains the full analysis pipeline, scripts, and supplementary
materials for a study testing whether a pre-specified cancer stem cell (CSC)
gene panel independently predicts colorectal cancer (CRC) prognosis, using
TCGA as a discovery cohort and three independent East Asian validation
cohorts (Korea, Japan, China).

## Study summary

- **Discovery cohort:** TCGA-COAD/READ (RNA-seq, N=[CONFIRM])
- **Validation cohorts:**
  | Cohort | Population | Platform | Endpoint | N |
  |---|---|---|---|---|
  | GSE107422 | Korea | RNA-seq (HiSeq 2000) | Binary recurrence | 110 |
  | GSE71187 | China | GPL6480 (Agilent) | OS (binary) | 52 |
  | GSE92921 | Japan | GPL570 (Affymetrix) | DFS (continuous risk-score model) | 59 |
  | GSE64857 | China | Microarray | Binary recurrence | 75 |
  | Zenodo 10.5281/zenodo.8333650 | Korea | RNA-seq | Binary recurrence | [CONFIRM] |

- **Key finding:** the CSC panel's apparent unadjusted association with TCGA
  overall survival (HR=0.871, 95% CI 0.731–1.038) is fully attenuated after
  adjusting for AJCC stage (HR=1.012, LRT p=0.92), indicating no independent
  prognostic value once stage is accounted for. This pattern was tested for
  consistency across five independent validation cohorts (see Results in the
  manuscript for cohort-by-cohort findings — reported separately, never
  pooled).

- **Primary signature used for the core manuscript claim:**
  `[CONFIRM BEFORE PUBLISHING — 11/20-gene curated ssGSEA panel score, per
  locked Global Rule 13 of the project playbook, NOT the 187-gene LASSO
  signature from script 06, unless that has since been validated as stable
  and explicitly adopted as primary. Verify this against the current
  manuscript draft before publishing this README.]`

## Repository structure

```
.
├── docs/
│   └── master_playbook.md       # Full project design doc + phase-by-phase protocol
├── scripts/
│   ├── 00_set_seed.R            # Fixed random seed (42), sourced by all scripts
│   ├── 01_download_tcga.R       # TCGA-COAD/READ download via TCGAbiolinks
│   ├── 02_download_geo_validation.R
│   ├── 02b_parse_geo_supplementary.R
│   ├── 02c_align_sample_names.R
│   ├── 03_qc_discovery.R / 03_qc_validation.R
│   ├── 04_deg_tcga.R            # DESeq2 differential expression, tumor vs normal
│   ├── 05_csc_panel_focus.R     # Curated 11-gene CSC panel differential expression
│   ├── 06_signature_construction.R  # LASSO-Cox / WGCNA-expanded signature — see caveat above
│   ├── 07_validation.R / 07b_validation_extended.R
│   ├── 08_enrichment.R          # GSEA (KEGG, GO-BP, Hallmark)
│   ├── 09_immune_infiltration.R
│   ├── 10_robustness_check.R    # Numbers ledger + audit
│   ├── 11-12_*                  # Forest plots, publication/supplementary figures
├── results/
│   ├── tables/
│   │   ├── numbers_ledger.csv   # Single source of truth for every statistic
│   │   │                        # reported in the manuscript — every number
│   │   │                        # anywhere else must trace back to this file
│   │   ├── DEG_full_TCGA.csv
│   │   ├── Table_S1_cohort_summary.csv
│   │   ├── Table_S2_excluded_samples.csv
│   │   └── validation_extended/ # Per-cohort validation results
│   └── figures/
├── manuscript/
│   └── manuscript_draft.md / .docx
├── supplementary/
│   ├── supplementary_materials.pdf
│   ├── methods_extended.md
│   ├── citation_audit.csv       # Every reference independently verified
│   │                             # against its real PubMed/DOI record
│   └── script_mapping_index.md  # Maps each script to the table/figure it produces
├── presentation/
│   └── conference_presentation.pptx
├── logs/
│   └── decisions.md             # Timestamped audit trail of every analytical
│                                 # decision made during the project, including
│                                 # dead ends, corrections, and dataset
│                                 # verification steps
├── environment.yml              # Pinned R/Python/Bioconductor package versions
└── README.md
```

## Reproducing this analysis

```bash
# 1. Create the environment (pinned versions — see environment.yml)
conda env create -f environment.yml
conda activate crc-csc-multicohort

# 2. Run the pipeline in order (each script writes to results/ and logs
#    its decisions to logs/decisions.md)
Rscript scripts/01_download_tcga.R
Rscript scripts/02_download_geo_validation.R
Rscript scripts/02b_parse_geo_supplementary.R
Rscript scripts/02c_align_sample_names.R
Rscript scripts/03_qc_discovery.R
Rscript scripts/03_qc_validation.R
Rscript scripts/04_deg_tcga.R
Rscript scripts/05_csc_panel_focus.R
Rscript scripts/06_signature_construction.R
Rscript scripts/07_validation.R
Rscript scripts/07b_validation_extended.R
Rscript scripts/08_enrichment.R
Rscript scripts/09_immune_infiltration.R
Rscript scripts/10_robustness_check.R
Rscript scripts/11_forest_plot.R
Rscript scripts/12_publication_figures.R

# 3. Verify: every number in results/tables/numbers_ledger.csv should match
#    the manuscript exactly. A robustness audit is run automatically by
#    script 10 and logged to logs/decisions.md.
```

Random seed is fixed at 42 throughout (`scripts/00_set_seed.R`). No cloud/
subagent execution was used for any analytical step in this project — all
scripts were run as direct, foreground, human-supervised executions.

## Data availability

All data used in this study are publicly available:

- **TCGA-COAD/READ:** [GDC Data Portal](https://portal.gdc.cancer.gov/), accessed via `TCGAbiolinks`
- **GSE107422:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107422
- **GSE71187:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71187
- **GSE92921:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92921
- **GSE64857:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64857
- **Zenodo (Korean RNA-seq cohort):** https://doi.org/10.5281/zenodo.8333650

No new data were generated. This is a fully secondary-data analysis.

## Methodology notes worth knowing before reusing this pipeline

- **Cross-platform validation, not meta-analysis:** each validation cohort is
  normalized independently and tested separately with a locked model from
  discovery. Results are never pooled/combined into a single meta-analytic
  estimate — see `docs/master_playbook.md` Global Rule 12 for the rationale.
- **Stage-adjustment is central to the primary finding**, not a side
  robustness check — the unadjusted CSC panel association (HR=0.871)
  disappears after adjusting for AJCC stage (HR=1.012). Any reuse of this
  pipeline for a different marker panel should replicate this adjustment
  step, since stage confounding is common in this class of analysis.
- **GSE92921 is reported with a continuous risk-score model only**, not a
  binary high/low split — the binary split was found to be critically
  unstable (only 1 event in the reference group; leave-one-out sensitivity
  showed the HR could swing by several orders of magnitude by removing a
  single patient). See `logs/decisions.md` for the full instability audit.
- **Citation verification:** every reference in the manuscript was
  independently checked against its real PubMed/DOI record before
  inclusion (`supplementary/citation_audit.csv`) — this was added after an
  earlier draft was found to contain a misattributed citation.

## Software

R (>=4.3), Bioconductor (DESeq2, TCGAbiolinks, GEOquery, limma, WGCNA,
clusterProfiler, org.Hs.eg.db, glmnet, survival, survminer, fgsea), Python
(pandas, numpy, scikit-learn, python-pptx). Full pinned versions in
`environment.yml`.

## Citation

If you use this code or data, please cite:

```
[CONFIRM BEFORE PUBLISHING — full citation once manuscript is accepted/
published. Include a preprint DOI here if posted before formal publication.]
```

## License

`[Choose and confirm: e.g., MIT for code, CC-BY-4.0 for supplementary
materials/figures — check your target journal's data/code sharing policy
first, some require a specific license.]`

## Contact

`[Author name(s), affiliation, and email — withheld here as placeholder;
fill in before making the repository public if not already public.]`
