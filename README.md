# Azadeh Golduzian — Research Portfolio

Public code highlights from my research as a Ph.D. candidate in Computational Biology at the University of New Mexico. My work bridges plasma-based biomarker modeling for Alzheimer's disease, neuroinflammation and blood-brain barrier research, and large-scale genomic annotation.

---

## Repository Structure

### `alzheimers_biomarkers/`
Biomarker modeling for Alzheimer's disease and blood-brain barrier permeability.
- `pTau217_cutoff_model.Rmd` — Logistic regression + ROC cut-point analysis for plasma pTau217 predicting Aβ positivity (simulated data)
- `ad_biomarker_analysis.R` — Random Forest and stepwise AIC models predicting albumin index and 3T MRI permeability from CSF, plasma, and MRI biomarkers (simulated data)

### `inflammation_bbb/`
Full analysis pipeline for an in-preparation paper on neuroinflammation and blood-brain barrier disruption in Alzheimer's disease and related dementias (UNM ADRC cohort). Real data is restricted access; analysis files use simulated data for demonstration.

| File | Description |
|------|-------------|
| `00_data_prep.qmd` | Variable selection, cohort subsetting, missing data visualization |
| `01_univariate.qmd` | Univariate biomarker associations with BBB disruption |
| `02_multivariate.qmd` | Multivariate regression with triple dichotomy and refined axes |
| `03_csf_vs_plasma.qmd` | CSF vs plasma biomarker scatterplots by diagnostic group |
| `04_composite_score.qmd` | Composite inflammatory score: BBB disruption factor |
| `04b_composite_score_corrfirst.qmd` | Composite score: correlation-first approach |
| `05_pca.qmd` | PCA-derived composite inflammatory score |
| `06_bbb_cognition.qmd` | BBB permeability and cognitive domain relationships |
| `07_multivariate_groups.qmd` | Multivariate analysis stratified by diagnostic group |
| `biol492_bbb_conceptual_model.qmd` | Conceptual two-pathway model of BBB disruption (course paper) |

### `genomics_annotation/`
Functional annotation pipeline for large-scale variant classification from UK Biobank and All of Us.
- `cisencode_transeqtl.qmd` — Intersects array variants with ENCODE cCREs (cis) and GTEx trans-eQTL windows; outputs CIS, TRANS, and REG variant sets

### `presentations/`
Posters and slide decks from public presentations.
- `Plasma_AD_VD_poster_2025.pdf` — Poster presented at a 2025 symposium on plasma biomarkers in Alzheimer's disease and vascular dementia

---

## Notes

- All scripts use **simulated data** that mirrors the real variable structure. No raw patient data, identifiers, or restricted materials are included.
- The `inflammation_bbb/` analyses require the `erikmisc` package: `remotes::install_github("erikerhardt/erikmisc")`
- Run `00_data_prep.qmd` first to generate the `.RData` files loaded by downstream analysis scripts.

---

📧 agolduzian96@unm.edu
🔗 [LinkedIn](https://www.linkedin.com/in/azadeh-golduzian-48236818b/) | [GitHub](https://github.com/AzadehGolduzian)
