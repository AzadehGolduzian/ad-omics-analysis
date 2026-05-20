# Neuroinflammation & Blood-Brain Barrier Analysis

Full Quarto analysis pipeline supporting an in-preparation manuscript on neuroinflammation and BBB disruption in Alzheimer's disease and related dementias (UNM ADRC cohort, Golduzian & Rosenberg et al.).

## Study Context

We examine how CSF and plasma inflammatory biomarkers (cytokines, MMPs, angiogenic factors) relate to blood-brain barrier permeability (albumin index, 3T MRI permeability) across diagnostic groups: AD, LA (Large Artery), SIVD, and Mixed dementia.

## Analysis Pipeline

Run scripts in order — each downstream file loads the `.RData` saved by `00_data_prep.qmd`.

| Step | File | Description |
|------|------|-------------|
| 0 | `00_data_prep.qmd` | Cohort subsetting, variable selection, baseline visit filter, missing data plot |
| 1 | `01_univariate.qmd` | Univariate associations between each biomarker and BBB outcomes |
| 2 | `02_multivariate.qmd` | Multivariate regression; triple dichotomy visualization |
| 3 | `03_csf_vs_plasma.qmd` | CSF vs plasma concordance scatterplots by DX group |
| 4 | `04_composite_score.qmd` | Composite inflammatory score from BBB-disruption predictors |
| 4b | `04b_composite_score_corrfirst.qmd` | Composite score using correlation-first approach |
| 5 | `05_pca.qmd` | PCA-based composite inflammatory score |
| 6 | `06_bbb_cognition.qmd` | BBB permeability vs neuropsychological test scores |
| 7 | `07_multivariate_groups.qmd` | Multivariate models stratified by diagnostic group |
| — | `biol492_bbb_conceptual_model.qmd` | Conceptual two-pathway model paper (Biol 492, PDF output) |

## Key Biomarkers

- **CSF:** Albumin index, IgG index, MMPs (1/2/3/9/10), angiogenic factors (VEGF-A/C/D, TIE2, bFGF), cytokines (IL-1β, IL-6, IL-8, IL-10, IFN-γ, TNF-α), NfL, GFAP, pTau, Aβ42/40
- **Plasma:** Matching MMP, angiogenic, and cytokine panels; pTau217, pTau181, Aβ40/42, NfL, GFAP
- **MRI:** 3T permeability, DTI-PSMD, free water, hippocampal volume, cortical thickness

## Dependencies

```r
# CRAN
install.packages(c("tidyverse", "ggplot2", "patchwork", "corrplot",
                   "broom", "car", "emmeans", "pROC", "labelled",
                   "dunn.test", "ggpubr", "modelsummary"))

# GitHub
remotes::install_github("erikerhardt/erikmisc")
```

## Data Note

Real data is from the UNM Alzheimer's Disease Research Center (ADRC) — restricted access clinical data. All `.qmd` files in this folder use **simulated data** (`00_data_prep.qmd` generates it) with the same variable names, types, and scales as the real dataset.
