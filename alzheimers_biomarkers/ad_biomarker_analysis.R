# ============================================================
# AD Biomarker Analysis: Albumin Index & BBB Permeability
# ============================================================
# Author  : Azadeh Golduzian | PhD Candidate, UNM
# Lab     : UNM Memory & Aging Center (ADRC)
# Date    : 2025
#
# Research question:
#   Which CSF, plasma, and MRI biomarkers predict
#   blood-brain barrier permeability in aging/dementia?
#
# Outcomes (response variables):
#   - fl_albumin_index_log2   : Albumin Index (CSF/plasma ratio, log2)
#   - mri_perm_3t_tot_log2    : 3T MRI white matter permeability (log2)
#
# Predictor groups:
#   - CSF biomarkers  : cytokines, MMPs, angiogenic factors, AD markers
#   - Plasma biomarkers: cytokines, MMPs, angiogenic factors, AD markers
#   - MRI measures    : DTI, hippocampal volume, cortical thickness
#
# Methods:
#   - Random Forest (randomForest + missForest for missing data)
#   - Stepwise AIC variable selection (bidirectional)
#   - Linear regression with age + sex covariates and interactions
#   - LM diagnostics (residuals, QQ, scale-location, leverage)
#
# NOTE: Real data not included (UNM ADRC clinical data, access restricted).
#       This script uses SIMULATED data with the same variable structure.
#       All variable names, types, and scales match the real dataset.
# ============================================================


# ── 0. Packages ──────────────────────────────────────────────────────────────
library(tidyverse)
library(randomForest)
library(missForest)
library(broom)
library(car)

# ── 1. Simulated Data ─────────────────────────────────────────────────────────
# Simulates the structure of the real ADRC dataset.
# Real data not included — controlled access clinical data.
# n = 120 participants (approximate real sample size)

set.seed(42)
n <- 120

# Helper: generate realistic log2-scale biomarker values
sim_log2 <- function(n, mean = 5, sd = 1.2) rnorm(n, mean, sd)

# Introduce ~15% missing data to mimic real clinical dataset
add_missing <- function(x, prop = 0.15) {
  x[sample(length(x), size = floor(prop * length(x)))] <- NA
  x
}

dat_sheet <- tibble(
  # Demographics
  age = rnorm(n, mean = 72, sd = 8),
  sex = sample(c("Female", "Male"), n, replace = TRUE),

  # ── Outcome variables ──────────────────────────────────────────────────────
  fl_albumin_index_log2   = add_missing(sim_log2(n, 4.2, 1.0)),
  mri_perm_3t_tot_log2    = add_missing(sim_log2(n, 3.8, 0.9)),

  # ── CSF Permeability / AD markers ─────────────────────────────────────────
  fl_igg_index_log2        = add_missing(sim_log2(n)),
  fl_csf_igg_log2          = add_missing(sim_log2(n)),
  fl_abeta_4240_ratio_log2 = add_missing(sim_log2(n, 3.5, 0.8)),
  fl_csf__ALL__ptau_log2   = add_missing(sim_log2(n, 4.0, 1.1)),

  # ── CSF MMPs ──────────────────────────────────────────────────────────────
  fl_msd_csf_mmp_1_log2    = add_missing(sim_log2(n)),
  fl_msd_csf_mmp_2_log2    = add_missing(sim_log2(n)),
  fl_msd_csf_mmp_3_log2    = add_missing(sim_log2(n)),
  fl_msd_csf_mmp_9_log2    = add_missing(sim_log2(n)),
  fl_csf__ALL__mmp_10_log2 = add_missing(sim_log2(n)),

  # ── CSF Angiogenic factors ────────────────────────────────────────────────
  fl_angio_csf_bfgf_log2   = add_missing(sim_log2(n)),
  fl_angio_csf_flt_1_log2  = add_missing(sim_log2(n)),
  fl_angio_csf_tie2_log2   = add_missing(sim_log2(n)),
  fl_angio_csf_vegf_a_log2 = add_missing(sim_log2(n)),
  fl_angio_csf_vegf_c_log2 = add_missing(sim_log2(n)),
  fl_angio_csf_vegf_d_log2 = add_missing(sim_log2(n)),

  # ── CSF Cytokines ─────────────────────────────────────────────────────────
  fl_cyto_csf_ifn_g_log2      = add_missing(sim_log2(n)),
  fl_csf__ALL__il_1b_log2     = add_missing(sim_log2(n)),
  fl_cyto_csf_il_6_log2       = add_missing(sim_log2(n)),
  fl_cyto_csf_il_8_log2       = add_missing(sim_log2(n)),
  fl_cyto_csf_il_10_log2      = add_missing(sim_log2(n)),
  fl_cyto_csf_il_12p70_log2   = add_missing(sim_log2(n)),
  fl_cyto_csf_tnf_a_log2      = add_missing(sim_log2(n)),

  # ── CSF Other ─────────────────────────────────────────────────────────────
  fl_csf__ALL__nfl_log2       = add_missing(sim_log2(n)),
  fl_csf__ALL__plgf_log2      = add_missing(sim_log2(n)),
  fl_csf__ALL__gfap_log2      = add_missing(sim_log2(n)),
  fl_myelin_basic_protein_log2 = add_missing(sim_log2(n)),

  # ── Plasma MMPs ───────────────────────────────────────────────────────────
  fl_plasma__ALL__mmp_1_log2  = add_missing(sim_log2(n)),
  fl_plasma__ALL__mmp_2_log2  = add_missing(sim_log2(n)),
  fl_plasma__ALL__mmp_3_log2  = add_missing(sim_log2(n)),
  fl_plasma__ALL__mmp_9_log2  = add_missing(sim_log2(n)),
  fl_plasma__ALL__mmp_10_log2 = add_missing(sim_log2(n)),

  # ── Plasma Angiogenic factors ─────────────────────────────────────────────
  fl_plasma__ALL__bfgf_log2   = add_missing(sim_log2(n)),
  fl_plasma__ALL__flt_1_log2  = add_missing(sim_log2(n)),
  fl_plasma__ALL__tie2_log2   = add_missing(sim_log2(n)),
  fl_plasma__ALL__vegf_a_log2 = add_missing(sim_log2(n)),
  fl_plasma__ALL__vegf_c_log2 = add_missing(sim_log2(n)),
  fl_plasma__ALL__vegf_d_log2 = add_missing(sim_log2(n)),

  # ── Plasma Cytokines ──────────────────────────────────────────────────────
  fl_cyto_heparin_plasma_ifn_g_log2    = add_missing(sim_log2(n)),
  fl_plasma__ALL__il_1b_log2           = add_missing(sim_log2(n)),
  fl_cyto_heparin_plasma_il_10_log2    = add_missing(sim_log2(n)),
  fl_cyto_heparin_plasma_il_13_log2    = add_missing(sim_log2(n)),

  # ── Plasma AD markers ─────────────────────────────────────────────────────
  fl_plasma__ALL__abeta_40_log2  = add_missing(sim_log2(n, 3.5, 0.8)),
  fl_plasma__ALL__abeta_42_log2  = add_missing(sim_log2(n, 3.2, 0.9)),
  fl_plasma__ALL__ptau_181_log2  = add_missing(sim_log2(n, 4.1, 1.0)),

  # ── Plasma Other ──────────────────────────────────────────────────────────
  fl_serum_igg_log2              = add_missing(sim_log2(n)),
  fl_plasma__ALL__plgf_log2      = add_missing(sim_log2(n)),

  # ── MRI measures ──────────────────────────────────────────────────────────
  mri_dti_psmd             = add_missing(rnorm(n, 0.002, 0.0005)),
  mri_dti_freewater_mean   = add_missing(rnorm(n, 0.18,  0.04)),
  mri_percent_vol_hippo    = add_missing(rnorm(n, 0.25,  0.05)),
  mri_mean_cort_thick_mm   = add_missing(rnorm(n, 2.5,   0.3))
)

cat("Dataset created:", nrow(dat_sheet), "rows,", ncol(dat_sheet), "columns\n")
cat("NOTE: This is SIMULATED data. Real ADRC data not included.\n\n")


# ── 2. Helper functions ───────────────────────────────────────────────────────

# Replacement for lab's custom e_plot_lm_diagnostics()
plot_lm_diagnostics <- function(fit, title = "") {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(fit, main = title)
  par(mfrow = c(1, 1))
}

# Select best complete-case subset from a variable list
# Replacement for lab's custom e_plot_complete_by_variable_subset()
select_complete_subset <- function(dat, var_list, var_resp, min_n = 20) {
  results <- tibble(
    n_vars    = integer(),
    n_complete = integer(),
    vars      = list()
  )
  for (k in seq_along(var_list)) {
    vars_k   <- c(var_resp, var_list[seq_len(k)])
    n_comp_k <- sum(complete.cases(dat[, vars_k]))
    results <- add_row(results,
                       n_vars     = k,
                       n_complete = n_comp_k,
                       vars       = list(vars_k))
  }
  # Return subset with most variables that still has >= min_n complete cases
  best <- results |> filter(n_complete >= min_n) |> slice_max(n_vars, n = 1)
  list(
    var_names  = best$vars[[1]],
    n_complete = best$n_complete
  )
}


# ── 3. Run stepwise AIC model ─────────────────────────────────────────────────
run_aic_model <- function(dat, var_resp, pred_vars) {

  dat_lm <- dat |>
    dplyr::select(all_of(c(var_resp, pred_vars, "age", "sex"))) |>
    tidyr::drop_na()

  cat("  N (complete cases):", nrow(dat_lm), "\n")

  form_main  <- as.formula(paste(var_resp, "~ age + sex +",
                                 paste(pred_vars, collapse = " + ")))
  form_lower <- as.formula(paste(var_resp, "~ 1"))

  fit_step <- step(
    lm(form_main, data = dat_lm),
    scope     = list(lower = form_lower, upper = form_main),
    direction = "both",
    k         = 2,     # AIC
    trace     = 0
  )

  list(fit = fit_step, dat = dat_lm)
}


# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1 — Albumin Index ~ CSF biomarkers + MRI
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== ANALYSIS 1: Albumin Index ~ CSF biomarkers + MRI ===\n")

sig_csf_alb <- c(
  "fl_abeta_4240_ratio_log2",
  "fl_csf__ALL__ptau_log2",
  "fl_angio_csf_bfgf_log2",
  "fl_angio_csf_flt_1_log2",
  "fl_angio_csf_tie2_log2",
  "fl_cyto_csf_ifn_g_log2",
  "fl_csf__ALL__il_1b_log2",
  "fl_cyto_csf_il_6_log2",
  "fl_cyto_csf_il_8_log2",
  "fl_cyto_csf_il_10_log2",
  "fl_msd_csf_mmp_1_log2",
  "fl_msd_csf_mmp_2_log2",
  "fl_msd_csf_mmp_3_log2",
  "fl_csf__ALL__mmp_10_log2",
  "fl_csf__ALL__plgf_log2",
  "fl_csf__ALL__gfap_log2",
  "fl_csf__ALL__nfl_log2",
  "mri_dti_psmd",
  "mri_dti_freewater_mean",
  "mri_percent_vol_hippo",
  "mri_mean_cort_thick_mm"
)

# Select best subset based on completeness
subset_1 <- select_complete_subset(dat_sheet, sig_csf_alb, "fl_albumin_index_log2")
pred_1   <- subset_1$var_names[subset_1$var_names != "fl_albumin_index_log2"]
cat("Variables selected:", length(pred_1), "| Complete cases:", subset_1$n_complete, "\n")

# Random Forest
dat_rf1 <- dat_sheet |>
  select(all_of(c("fl_albumin_index_log2", pred_1))) |>
  tidyr::drop_na()

set.seed(123)
rf1 <- randomForest(fl_albumin_index_log2 ~ ., data = dat_rf1,
                    ntree = 500, importance = TRUE)
cat("\nRandom Forest — Albumin ~ CSF:\n")
print(rf1)
varImpPlot(rf1, main = "Analysis 1: Variable Importance\nAlbumin Index ~ CSF + MRI")

# Stepwise AIC
res1  <- run_aic_model(dat_sheet, "fl_albumin_index_log2", pred_1)
cat("\nStepwise AIC model summary:\n")
print(summary(res1$fit))
plot_lm_diagnostics(res1$fit, "Analysis 1: Albumin Index ~ CSF + MRI")


# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2 — 3T Permeability ~ CSF biomarkers + MRI
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== ANALYSIS 2: 3T Permeability ~ CSF biomarkers + MRI ===\n")

sig_csf_perm <- c(
  "fl_angio_csf_vegf_c_log2",
  "fl_cyto_csf_ifn_g_log2",
  "fl_cyto_csf_il_10_log2",
  "fl_csf__ALL__gfap_log2",
  "fl_csf_igg_log2",
  "fl_csf__ALL__il_1b_log2",
  "fl_msd_csf_mmp_2_log2",
  "fl_csf__ALL__nfl_log2",
  "mri_dti_psmd",
  "mri_dti_freewater_mean"
)

subset_2 <- select_complete_subset(dat_sheet, sig_csf_perm, "mri_perm_3t_tot_log2")
pred_2   <- subset_2$var_names[subset_2$var_names != "mri_perm_3t_tot_log2"]
cat("Variables selected:", length(pred_2), "| Complete cases:", subset_2$n_complete, "\n")

# Random Forest
dat_rf2 <- dat_sheet |>
  select(all_of(c("mri_perm_3t_tot_log2", pred_2))) |>
  tidyr::drop_na()

set.seed(123)
rf2 <- randomForest(mri_perm_3t_tot_log2 ~ ., data = dat_rf2,
                    ntree = 500, importance = TRUE)
cat("\nRandom Forest — Permeability ~ CSF:\n")
print(rf2)
varImpPlot(rf2, main = "Analysis 2: Variable Importance\n3T Permeability ~ CSF + MRI")

# Stepwise AIC
res2 <- run_aic_model(dat_sheet, "mri_perm_3t_tot_log2", pred_2)
cat("\nStepwise AIC model summary:\n")
print(summary(res2$fit))
plot_lm_diagnostics(res2$fit, "Analysis 2: 3T Permeability ~ CSF + MRI")


# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3 — Albumin Index ~ Plasma biomarkers
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== ANALYSIS 3: Albumin Index ~ Plasma biomarkers ===\n")

sig_plasma_alb <- c(
  "fl_plasma__ALL__abeta_40_log2",
  "fl_plasma__ALL__abeta_42_log2",
  "fl_plasma__ALL__vegf_d_log2",
  "fl_plasma__ALL__il_1b_log2",
  "fl_cyto_heparin_plasma_il_10_log2",
  "fl_plasma__ALL__mmp_2_log2",
  "fl_plasma__ALL__mmp_3_log2",
  "fl_plasma__ALL__mmp_10_log2",
  "fl_serum_igg_log2"
)

subset_3 <- select_complete_subset(dat_sheet, sig_plasma_alb, "fl_albumin_index_log2")
pred_3   <- subset_3$var_names[subset_3$var_names != "fl_albumin_index_log2"]
cat("Variables selected:", length(pred_3), "| Complete cases:", subset_3$n_complete, "\n")

# Random Forest with missForest imputation
dat_rf3_raw <- dat_sheet |>
  select(all_of(c("fl_albumin_index_log2", pred_3))) |>
  filter(rowSums(!is.na(.)) > 0) |>
  mutate(across(everything(), ~ as.numeric(as.character(.))))

set.seed(123)
imp3         <- missForest(as.data.frame(dat_rf3_raw), ntree = 100)
dat_rf3_imp  <- imp3$ximp

set.seed(123)
rf3 <- randomForest(fl_albumin_index_log2 ~ ., data = dat_rf3_imp,
                    ntree = 500, importance = TRUE)
cat("\nRandom Forest (imputed) — Albumin ~ Plasma:\n")
print(rf3)
varImpPlot(rf3, main = "Analysis 3: Variable Importance\nAlbumin Index ~ Plasma (imputed)")

# Stepwise AIC
res3 <- run_aic_model(dat_sheet, "fl_albumin_index_log2", pred_3)
cat("\nStepwise AIC model summary:\n")
print(summary(res3$fit))
cat("\nType III ANOVA:\n")
print(car::Anova(res3$fit, type = 3))
cat("\nSorted coefficients:\n")
print(broom::tidy(res3$fit) |> arrange(p.value), n = Inf)
plot_lm_diagnostics(res3$fit, "Analysis 3: Albumin Index ~ Plasma")


# ════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4 — 3T Permeability ~ Plasma biomarkers
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== ANALYSIS 4: 3T Permeability ~ Plasma biomarkers ===\n")

sig_plasma_perm <- c(
  "fl_plasma__ALL__ptau_181_log2",
  "fl_plasma__ALL__abeta_40_log2",
  "fl_plasma__ALL__vegf_a_log2",
  "fl_plasma__ALL__vegf_c_log2",
  "fl_cyto_heparin_plasma_ifn_g_log2",
  "fl_cyto_heparin_plasma_il_13_log2",
  "fl_plasma__ALL__mmp_9_log2",
  "fl_plasma__ALL__plgf_log2"
)

subset_4 <- select_complete_subset(dat_sheet, sig_plasma_perm, "mri_perm_3t_tot_log2")
pred_4   <- subset_4$var_names[subset_4$var_names != "mri_perm_3t_tot_log2"]
cat("Variables selected:", length(pred_4), "| Complete cases:", subset_4$n_complete, "\n")

# Random Forest
dat_rf4 <- dat_sheet |>
  select(all_of(c("mri_perm_3t_tot_log2", pred_4))) |>
  tidyr::drop_na()

set.seed(123)
rf4 <- randomForest(mri_perm_3t_tot_log2 ~ ., data = dat_rf4,
                    ntree = 500, importance = TRUE)
cat("\nRandom Forest — Permeability ~ Plasma:\n")
print(rf4)
varImpPlot(rf4, main = "Analysis 4: Variable Importance\n3T Permeability ~ Plasma")

# Stepwise AIC
res4 <- run_aic_model(dat_sheet, "mri_perm_3t_tot_log2", pred_4)
cat("\nStepwise AIC model summary:\n")
print(summary(res4$fit))
cat("\nType III ANOVA:\n")
print(car::Anova(res4$fit, type = 3))
cat("\nSorted coefficients:\n")
print(broom::tidy(res4$fit) |> arrange(p.value), n = Inf)
plot_lm_diagnostics(res4$fit, "Analysis 4: 3T Permeability ~ Plasma")


# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("Four models fitted:\n")
cat("  1. Albumin Index    ~ CSF biomarkers + MRI\n")
cat("  2. 3T Permeability  ~ CSF biomarkers + MRI\n")
cat("  3. Albumin Index    ~ Plasma biomarkers\n")
cat("  4. 3T Permeability  ~ Plasma biomarkers\n")
cat("\nAll models use AIC-based stepwise selection (bidirectional),\n")
cat("adjusted for age and sex.\n")
cat("\nNote: Results above are from SIMULATED data.\n")
cat("Real results reported in Golduzian et al. (in preparation).\n")
