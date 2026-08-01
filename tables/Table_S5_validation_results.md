## Table S5 — Extended Validation Results Across Asian CRC Cohorts

### Summary Table

| Cohort | Population | N | Endpoint | Matched Genes | Method | Point Estimate | 95% CI | p-value | AUC | Outcome |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE107422 | Korean | 110 | Recurrence | 187/187 | Logistic regression | OR = 3.246 | 1.399 – 7.535 | 0.006 | 0.640 | **SUCCESS** |
| GSE71187 | Chinese | 52 | OS (binary) | 139/187 | Logistic regression | OR = 2.000 | 0.647 – 6.185 | 0.229 | 0.583 | **PARTIAL** |
| GSE64857 | Chinese | 75 | Recurrence (binary) | 147/187 | Logistic regression | OR = 1.333 | 0.513 – 3.463 | 0.555 | 0.536 | **PARTIAL** |
| Zenodo #8333650 | Korean | 176 | Recurrence (binary) | 117/187 | Logistic regression | OR = 1.889 | 0.957 – 3.727 | 0.067 | 0.578 | **PARTIAL** |

**Classification criteria:** SUCCESS (p < 0.05 AND OR > 1 AND AUC > 0.5); PARTIAL (p < 0.10 OR [OR > 1 AND AUC > 0.5]); FAIL (OR ≤ 1 OR AUC ≤ 0.5).

---

### Cohort-Level Details

#### 1. GSE107422 (Korean CRC, N = 110) — SUCCESS

**Primary Analysis: Median Risk Score Split**

| Metric | Value | 95% CI | p-value |
|---|---|---|---|
| Odds Ratio (High vs Low) | 3.246 | 1.399 – 7.535 | 0.0061 |
| AUC | 0.640 | — | — |

**Logistic regression model:**
```
glm(formula = recur_event ~ risk_group_z, family = binomial)
Coefficients:
                   Estimate Std. Error z value Pr(>|z|)
(Intercept)       -1.3157     0.3396  -3.875 0.000107
risk_group_zHigh   1.1775     0.4296   2.741 0.006131
```
- Null deviance: 141.81 on 109 df
- Residual deviance: 133.79 on 108 df
- AIC: 137.79
- High-risk: 58 patients, Low-risk: 52 patients, Events: 38

**Sensitivity Analysis: Tertile-Based Stratification**

Tertile cutoffs: −0.3224, 0.3270 (based on TCGA risk score distribution)

| Group | N | OR (vs Low) | 95% CI | p-value |
|---|---|---|---|---|
| Low | 37 | 1.000 (ref) | — | — |
| Medium | 36 | 3.289 | 1.094 – 9.895 | 0.034 |
| High | 37 | 4.893 | 1.651 – 14.500 | 0.004 |

**Trend test:** OR per tertile = 2.125 (95% CI: 1.268 – 3.562), p = 0.0042
**Median split AUC:** 0.640; **Tertile split AUC:** 0.660

---

#### 2. GSE71187 (Chinese CRC, Shanghai, N = 52) — PARTIAL

| Metric | Value | 95% CI | p-value |
|---|---|---|---|
| Odds Ratio (High vs Low) | 2.000 | 0.647 – 6.185 | 0.229 |
| AUC | 0.583 | — | — |

**Logistic regression model:**
```
glm(formula = status ~ risk_group, family = binomial)
Coefficients:
                   Estimate Std. Error z value Pr(>|z|)
(Intercept)       -0.5978     0.3754  -1.593    0.111
risk_groupHigh     0.6931     0.5760   1.203    0.229
```
- Null deviance: 70.852 on 51 df
- Residual deviance: 69.389 on 50 df
- AIC: 73.389
- Events: 22
- Matched signature genes: 139 / 187
- Endpoint: OS (binary status; no time-to-event data available)

**Note:** The odds ratio direction is consistent with the TCGA signature (OR > 1), but the result is not statistically significant. The wide confidence interval reflects the limited sample size (N = 52, 22 events). Classified as PARTIAL because OR > 1 and AUC > 0.5 but p ≥ 0.05.

---

#### 3. GSE64857 (Chinese CRC, Fudan University, N = 75) — PARTIAL

| Metric | Value | 95% CI | p-value |
|---|---|---|---|
| Odds Ratio (High vs Low) | 1.333 | 0.513 – 3.463 | 0.555 |
| AUC | 0.536 | — | — |

**Logistic regression model:**
```
glm(formula = recur_status ~ risk_group, family = binomial)
Coefficients:
                   Estimate Std. Error z value Pr(>|z|)
(Intercept)       -0.7673     0.3356  -2.286   0.0222
risk_groupHigh     0.2877     0.4870   0.591   0.5547
```
- Null deviance: 96.804 on 74 df
- Residual deviance: 96.455 on 73 df
- AIC: 100.45
- Recurred / Non-recurrent: 26 / 49
- Matched signature genes: 147 / 187

**Note:** The odds ratio is directionally consistent (OR > 1) but not statistically significant. The modest AUC (0.536) indicates near-chance discriminative performance in this cohort. The limited gene match rate (147/187) on the Affymetrix platform may partly explain the attenuated signal. Classified as PARTIAL.

---

#### 4. Zenodo #8333650 (Korean CRC, Asan Historical Cohort, N = 176) — PARTIAL

| Metric | Value | 95% CI | p-value |
|---|---|---|---|
| Odds Ratio (High vs Low) | 1.889 | 0.957 – 3.727 | 0.067 |
| AUC (binary group) | 0.578 | — | — |
| AUC (continuous risk score) | 0.569 | — | — |

**Logistic regression model:**
```
glm(formula = recur_binary ~ risk_group, family = binomial)
Coefficients:
                   Estimate Std. Error z value Pr(>|z|)
(Intercept)       -1.3291     0.2651  -5.014 5.32e-07
risk_groupHigh     0.6360     0.3468   1.834   0.0667
```
- Null deviance: 206.26 on 175 df
- Residual deviance: 202.81 on 174 df
- AIC: 206.81
- Recurrence events: 48
- Matched signature genes: 117 / 187

**Note:** This cohort shows a borderline significant association (p = 0.067) with OR = 1.889, the strongest effect among the PARTIAL cohorts. The lower gene match rate (117/187) on a custom Affymetrix array may contribute to signal attenuation. Classified as PARTIAL because p < 0.10 but p ≥ 0.05.

---

### Validation Outcome Summary

- **SUCCESS (1/4 cohorts):** Direction consistent, statistically significant, AUC > 0.5
- **PARTIAL (3/4 cohorts):** Direction consistent but not statistically significant at α = 0.05; all AUCs > 0.5
- **All 4 cohorts** show OR > 1, indicating directionally consistent association between high risk score and adverse outcomes across diverse Asian populations and assay platforms

*Source files: `results/tables/validation_extended/GSE107422_outcome.txt`, `results/tables/validation_extended/GSE71187_outcome.txt`, `results/tables/validation_extended/GSE64857_outcome.txt`, `results/tables/validation_zenodo/Zenodo8333650_recurrence_outcome.txt`*
