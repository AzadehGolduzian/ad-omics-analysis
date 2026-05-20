# Alzheimer's Biomarkers

Biomarker modeling and cut-point analysis for Alzheimer's disease (AD) diagnosis and staging.

## Files

### `pTau217_cutoff_model.Rmd`
**Goal:** Identify the optimal plasma pTau217 threshold that classifies individuals as amyloid-beta positive or negative.

**Methods:**
- Logistic regression (pTau217 → Aβ42/40 positivity)
- ROC curve with AUC and Youden-optimal threshold
- Probability contrast plot with cut-point overlay
- Scatter plot of observed positivity vs pTau217

**Dependencies:** `erikmisc`, `pROC`, `ggplot2`, `patchwork`, `emmeans`

---

### `ad_biomarker_analysis.R`
**Goal:** Identify which CSF, plasma, and MRI biomarkers best predict blood-brain barrier permeability (albumin index, 3T MRI permeability).

**Methods:**
- Random Forest variable importance (`randomForest`, `missForest` for imputation)
- Bidirectional stepwise AIC model selection
- Linear regression with age + sex covariates
- LM diagnostics (residuals, Q-Q, leverage)

**Four analyses run:**
1. Albumin Index ~ CSF biomarkers + MRI
2. 3T Permeability ~ CSF biomarkers + MRI
3. Albumin Index ~ Plasma biomarkers
4. 3T Permeability ~ Plasma biomarkers

**Dependencies:** `tidyverse`, `randomForest`, `missForest`, `broom`, `car`

---

> All scripts use simulated data with the same variable structure as the real UNM ADRC dataset. No patient data included.
