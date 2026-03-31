#' ==============================================================================
#' Phase 8: Bayesian DPM — PRODUCTION RUN (Full Sample)
#'
#' Runs the Stan DPM model on the full analysis sample (~95K observations).
#' Designed for HPC cluster execution (requires ~64GB RAM, ~24-48 hours).
#'
#' Usage:
#'   Rscript phase8_bayesian_dpm_production.R [N_SUBSAMPLE]
#'   - No argument: full sample
#'   - With argument: subsample of that size (e.g., 50000)
#'
#' Requirements:
#'   - R 4.x with rstan, haven, dplyr
#'   - Stan toolchain (g++, make)
#'   - ~64GB RAM for full sample; ~32GB for N=50K
#'
#' Author: Alina Malkova
#' Date: March 2026
#' ==============================================================================

library(rstan)
library(haven)
library(dplyr)

# Parse command-line argument for subsample size
args <- commandArgs(trailingOnly = TRUE)
N_TARGET <- if (length(args) > 0) as.integer(args[1]) else NULL

# Configure Stan for cluster
n_cores <- min(4, parallel::detectCores())
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)

cat("============================================================\n")
cat("BAYESIAN DPM — PRODUCTION RUN\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Cores: %d\n", n_cores))
cat(sprintf("Target N: %s\n", ifelse(is.null(N_TARGET), "FULL SAMPLE", format(N_TARGET, big.mark = ","))))
cat("============================================================\n\n")

# ==============================================================================
# 1. Paths — adjust DATA_DIR for cluster if needed
# ==============================================================================

# Try cluster path first, then local
cluster_path <- Sys.getenv("DATA_DIR", unset = "")
if (nchar(cluster_path) > 0 && dir.exists(cluster_path)) {
  DATA_DIR <- cluster_path
} else {
  DATA_DIR <- "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data"
}

OUTPUT_DIR <- file.path(DATA_DIR, "output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Data directory: %s\n", DATA_DIR))
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

# ==============================================================================
# 2. Load and Prepare Data
# ==============================================================================

cat("Loading data...\n")
df <- read_dta(file.path(DATA_DIR, "analysis_dataset_with_se.dta"))
cat(sprintf("Raw dataset: %d observations\n", nrow(df)))

df <- df %>%
  filter(age >= 18 & age <= 64) %>%
  filter(employed == 1 | unemployed == 1) %>%
  filter(year >= 2013) %>%
  filter(cbsa > 0 & !is.na(cbsa)) %>%
  filter(!is.na(banking_mode)) %>%
  filter(!is.na(self_employed)) %>%
  filter(!is.na(pct_broadband)) %>%
  filter(!is.na(hsupwgtk)) %>%
  mutate(
    se = as.integer(self_employed == 1),
    branch = as.integer(banking_mode == 3),
    mobile = as.integer(banking_mode == 2),
    emp_status = case_when(
      wage_worker == 1 ~ 1,
      self_employed == 1 ~ 2,
      TRUE ~ 3
    ),
    choice = (banking_mode - 1) * 3 + emp_status,
    age_std = (age - mean(age)) / sd(age),
    female_ind = as.integer(!is.na(female) & female == 1),
    married_ind = as.integer(!is.na(married) & married == 1),
    educ_cat = case_when(
      no_hs == 1 ~ 1,
      hs_diploma == 1 ~ 2,
      some_college == 1 ~ 3,
      college_degree == 1 ~ 4,
      TRUE ~ 2
    )
  )

cat(sprintf("Analysis sample: %d observations\n", nrow(df)))
cat(sprintf("Baseline SE rate: %.2f%%\n", mean(df$se) * 100))

# Subsample if requested
if (!is.null(N_TARGET) && N_TARGET < nrow(df)) {
  set.seed(20260330)
  df <- df %>% sample_n(N_TARGET)
  cat(sprintf("Subsampled to: %d observations\n", nrow(df)))
}

# ==============================================================================
# 3. Prepare Stan Data
# ==============================================================================

X <- model.matrix(~ age_std + female_ind + married_ind + factor(educ_cat) + factor(year),
                  data = df)[, -1]

K_MAX <- 8  # Increased from 6 for production

stan_data <- list(
  N = nrow(df),
  J = 9,
  P = ncol(X),
  K_max = K_MAX,
  y = df$choice,
  X = X,
  branch = df$branch,
  mobile = df$mobile,
  branch_density = as.numeric(scale(df$pct_broadband)),
  weights = df$hsupwgtk / mean(df$hsupwgtk, na.rm = TRUE)
)

cat(sprintf("\nStan data:\n"))
cat(sprintf("  N = %s, J = %d, P = %d, K_max = %d\n",
            format(stan_data$N, big.mark = ","), stan_data$J, stan_data$P, K_MAX))
cat(sprintf("  Memory estimate: ~%.1f GB for gradients\n",
            stan_data$N * K_MAX * 9 * 8 / 1e9 * 3))

# Free memory
rm(df, X); gc(verbose = FALSE)

# ==============================================================================
# 4. Compile and Run Stan
# ==============================================================================

stan_file <- file.path(DATA_DIR, "phase8_bayesian_dpm.stan")
rds_file <- file.path(DATA_DIR, "phase8_bayesian_dpm.rds")

cat("\nCompiling Stan model...\n")
if (file.exists(rds_file)) {
  model <- readRDS(rds_file)
  cat("Loaded pre-compiled model from RDS.\n")
} else {
  model <- stan_model(stan_file)
  saveRDS(model, rds_file)
  cat("Model compiled and saved to RDS.\n")
}

cat("\n============================================================\n")
cat("STARTING MCMC SAMPLING\n")
cat(sprintf("Time: %s\n", Sys.time()))
cat("============================================================\n\n")

# Production MCMC settings
N_CHAINS <- 4
N_ITER <- 3000     # More iterations for convergence with large N
N_WARMUP <- 1500
N_THIN <- 2        # Keep every 2nd draw → 750 post-warmup draws per chain

fit <- sampling(
  model,
  data = stan_data,
  chains = N_CHAINS,
  iter = N_ITER,
  warmup = N_WARMUP,
  thin = N_THIN,
  seed = 20260330,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  )
)

cat(sprintf("\nSampling complete: %s\n", Sys.time()))

# ==============================================================================
# 5. Diagnostics
# ==============================================================================

cat("\n============================================================\n")
cat("MCMC DIAGNOSTICS\n")
cat("============================================================\n")

# Check for divergences
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
cat(sprintf("Divergent transitions: %d\n", n_divergent))

# Check treedepth
n_max_treedepth <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 12)))
cat(sprintf("Max treedepth hits: %d\n", n_max_treedepth))

# Summary of key parameters
cat("\nParameter summary (Rhat, n_eff):\n")
key_params <- c("pi", "gamma_branch", "sigma_gamma", "K_effective",
                "counterfactual_effect", "se_rate_baseline")
param_summary <- summary(fit, pars = key_params)$summary
rhat_range <- range(param_summary[, "Rhat"], na.rm = TRUE)
neff_range <- range(param_summary[, "n_eff"], na.rm = TRUE)
cat(sprintf("  Rhat range: [%.3f, %.3f]\n", rhat_range[1], rhat_range[2]))
cat(sprintf("  n_eff range: [%.0f, %.0f]\n", neff_range[1], neff_range[2]))

if (max(param_summary[, "Rhat"], na.rm = TRUE) > 1.1) {
  cat("  WARNING: Some Rhat > 1.1 — consider more iterations.\n")
}

# ==============================================================================
# 6. Extract Results
# ==============================================================================

cat("\n============================================================\n")
cat("POSTERIOR SUMMARY\n")
cat("============================================================\n")

post <- extract(fit)

# Mixing weights
pi_post <- post$pi
pi_mean <- colMeans(pi_post)

cat("\nMixing Weights (posterior mean):\n")
for (k in 1:K_MAX) {
  cat(sprintf("  Type %d: %.1f%% (95%% CI: [%.1f%%, %.1f%%])%s\n",
              k, pi_mean[k] * 100,
              quantile(pi_post[, k], 0.025) * 100,
              quantile(pi_post[, k], 0.975) * 100,
              ifelse(pi_mean[k] < 0.01, "  [inactive]", "")))
}

# Effective K
K_eff <- post$K_effective
cat(sprintf("\nEffective K: %.2f (95%% CI: [%.0f, %.0f])\n",
            mean(K_eff), quantile(K_eff, 0.025), quantile(K_eff, 0.975)))
cat(sprintf("P(K >= 4): %.3f\n", mean(K_eff >= 4)))
cat(sprintf("P(K >= 3): %.3f\n", mean(K_eff >= 3)))

# Branch effects
gamma_post <- post$gamma_branch
gamma_mean <- colMeans(gamma_post)
cat("\nBranch Effects by Type:\n")
for (k in 1:K_MAX) {
  if (pi_mean[k] > 0.01) {
    cat(sprintf("  Type %d: gamma = %.4f (95%% CI: [%.4f, %.4f])\n",
                k, gamma_mean[k],
                quantile(gamma_post[, k], 0.025),
                quantile(gamma_post[, k], 0.975)))
  }
}

# Counterfactual
cf <- post$counterfactual_effect
cat(sprintf("\nCounterfactual (50%% branch closure):\n"))
cat(sprintf("  Mean: %.2f%%\n", mean(cf)))
cat(sprintf("  Median: %.2f%%\n", median(cf)))
cat(sprintf("  SD: %.2f%%\n", sd(cf)))
cat(sprintf("  95%% CI: [%.2f%%, %.2f%%]\n",
            quantile(cf, 0.025), quantile(cf, 0.975)))
cat(sprintf("  P(effect < 0): %.4f\n", mean(cf < 0)))

# Baseline SE rate
se_base <- post$se_rate_baseline
cat(sprintf("\nBaseline SE rate: %.2f%% (95%% CI: [%.2f%%, %.2f%%])\n",
            mean(se_base) * 100,
            quantile(se_base, 0.025) * 100,
            quantile(se_base, 0.975) * 100))

# ==============================================================================
# 7. Comparison Table
# ==============================================================================

cat("\n============================================================\n")
cat("COMPARISON: Finite Mixture (K=4) vs Bayesian DPM\n")
cat("============================================================\n\n")
cat("                    | MNL K=4 (BIC) | Bayesian DPM (production)\n")
cat("--------------------|---------------|-------------------------\n")
cat(sprintf("Effective K         |       4       |    %.2f\n", mean(K_eff)))
cat(sprintf("P(K >= 4)           |     1.00      |    %.3f\n", mean(K_eff >= 4)))
cat(sprintf("CF Effect           |    -10.8%%     |   %.1f%%\n", mean(cf)))
cat(sprintf("CF 95%% CI           | [-14.1,-7.4]  | [%.1f, %.1f]\n",
            quantile(cf, 0.025), quantile(cf, 0.975)))
cat(sprintf("CI width            |     6.7pp     |   %.1fpp\n",
            quantile(cf, 0.975) - quantile(cf, 0.025)))

# ==============================================================================
# 8. Save Results
# ==============================================================================

# Main results
results <- data.frame(
  parameter = c("N_sample", "K_max", "n_chains", "n_iter", "n_warmup",
                "K_effective_mean", "K_effective_sd",
                "prob_K_ge_3", "prob_K_ge_4", "prob_K_ge_5",
                "cf_effect_mean", "cf_effect_median", "cf_effect_sd",
                "cf_ci_lower", "cf_ci_upper",
                "prob_negative",
                "se_rate_baseline", "se_rate_counterfactual",
                "n_divergent", "max_rhat", "min_neff"),
  value = c(stan_data$N, K_MAX, N_CHAINS, N_ITER, N_WARMUP,
            mean(K_eff), sd(K_eff),
            mean(K_eff >= 3), mean(K_eff >= 4), mean(K_eff >= 5),
            mean(cf), median(cf), sd(cf),
            quantile(cf, 0.025), quantile(cf, 0.975),
            mean(cf < 0),
            mean(se_base), mean(post$se_rate_counterfactual),
            n_divergent, max(param_summary[, "Rhat"], na.rm = TRUE),
            min(param_summary[, "n_eff"], na.rm = TRUE))
)

output_suffix <- ifelse(is.null(N_TARGET), "full", paste0("N", N_TARGET/1000, "k"))

write.csv(results,
          file.path(OUTPUT_DIR, paste0("phase8_dpm_production_", output_suffix, ".csv")),
          row.names = FALSE)

# Save full posterior draws for mixing weights and counterfactual
posterior_draws <- data.frame(
  draw = 1:length(cf),
  K_effective = K_eff,
  cf_effect = cf,
  se_baseline = se_base,
  se_counterfactual = post$se_rate_counterfactual
)
for (k in 1:K_MAX) {
  posterior_draws[[paste0("pi_", k)]] <- pi_post[, k]
  posterior_draws[[paste0("gamma_branch_", k)]] <- gamma_post[, k]
}

write.csv(posterior_draws,
          file.path(OUTPUT_DIR, paste0("phase8_dpm_production_posterior_", output_suffix, ".csv")),
          row.names = FALSE)

# Save Stan fit object for later analysis
saveRDS(fit, file.path(OUTPUT_DIR, paste0("phase8_dpm_production_fit_", output_suffix, ".rds")))

cat(sprintf("\nResults saved to:\n"))
cat(sprintf("  %s/phase8_dpm_production_%s.csv\n", OUTPUT_DIR, output_suffix))
cat(sprintf("  %s/phase8_dpm_production_posterior_%s.csv\n", OUTPUT_DIR, output_suffix))
cat(sprintf("  %s/phase8_dpm_production_fit_%s.rds\n", OUTPUT_DIR, output_suffix))

cat(sprintf("\n============================================================\n"))
cat(sprintf("PRODUCTION RUN COMPLETE: %s\n", Sys.time()))
cat(sprintf("============================================================\n"))
