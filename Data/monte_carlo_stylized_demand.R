#' ==============================================================================
#' Monte Carlo Simulation: Stylized Binary Demand Model (Comment 2b)
#'
#' Demonstrates that K sensitivity generalizes beyond banking/MNL:
#'   - Binary demand (buy/not-buy) with K=3 latent consumer types
#'   - "Partial cancellation" DGP: type-specific treatment effects partly cancel
#'   - Weighted average treatment effect near zero (0.023), masking heterogeneity
#'
#' Model:
#'   Y_i = 1{u_i > 0}, where
#'   u_i = alpha_k + beta_price_k * price_i + beta_quality_k * quality_i
#'         + gamma_k * treatment_i + epsilon_i,  epsilon ~ logistic
#'
#' Estimation: EM with binary logit components for K=1,...,5
#'
#' Output: CSV results + ggplot figure to Data/output/figures_methods/
#'
#' Author: Alina Malkova
#' Date: April 2026
#' ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(20260401)

# Paths
DATA_DIR <- "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data"
OUTPUT_DIR <- file.path(DATA_DIR, "output/figures_methods")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DATA_DIR, "output"), showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("MONTE CARLO: STYLIZED BINARY DEMAND MODEL\n")
cat("Comment 2b -- K sensitivity beyond banking\n")
cat("============================================================\n\n")


# ==============================================================================
# 1. DGP Parameters
# ==============================================================================

dgp_demand <- list(
  name       = "Demand",
  K_true     = 3,
  pi         = c(0.30, 0.40, 0.30),
  # Intercepts: type 1 price-sensitive (low baseline), type 3 quality-seeking
  alpha      = c(-1.0, 0.0, 0.5),
  beta_price = c(-0.40, -0.20, -0.10),
  beta_qual  = c( 0.10,  0.15,  0.30),
  # Treatment effects: partial cancellation
  gamma      = c(-0.15, 0.02, 0.20)
)

weighted_avg <- sum(dgp_demand$pi * dgp_demand$gamma)
cat(sprintf("DGP: K_true = %d, shares = (%s)\n",
            dgp_demand$K_true, paste(dgp_demand$pi, collapse = ", ")))
cat(sprintf("  gamma = (%s)\n", paste(dgp_demand$gamma, collapse = ", ")))
cat(sprintf("  Weighted average gamma = %.4f (near zero)\n", weighted_avg))
cat(sprintf("  Interpretation: treatment effect nearly vanishes in aggregate\n\n"))


# ==============================================================================
# 2. Data Generating Process
# ==============================================================================

simulate_demand_data <- function(N, dgp) {
  K <- dgp$K_true

  # Assign consumer types
  type <- sample(1:K, N, replace = TRUE, prob = dgp$pi)

  # Covariates
  price     <- runif(N, 1, 10)
  quality   <- rnorm(N, 5, 1)
  treatment <- rbinom(N, 1, 0.5)

  # Utility (type-specific parameters)
  u <- dgp$alpha[type] +
       dgp$beta_price[type] * price +
       dgp$beta_qual[type]  * quality +
       dgp$gamma[type]      * treatment

  # Logistic error -> binary logit
  eps <- rlogis(N)
  Y   <- as.integer((u + eps) > 0)

  data.frame(Y = Y, price = price, quality = quality,
             treatment = treatment, type = type)
}

# Quick sanity check
test_data <- simulate_demand_data(50000, dgp_demand)
cat(sprintf("Sanity check (N=50,000): purchase rate = %.3f\n", mean(test_data$Y)))
cat(sprintf("  By type: %s\n",
            paste(round(tapply(test_data$Y, test_data$type, mean), 3), collapse = ", ")))
cat(sprintf("  Treatment rate: %.3f\n\n", mean(test_data$treatment)))
rm(test_data)


# ==============================================================================
# 3. Binary Logit Log-Likelihood (vectorized)
# ==============================================================================

#' Compute log-likelihood contributions for binary logit
#' @param Y     N-vector of 0/1 outcomes
#' @param Xbeta N-vector of linear index (alpha + X*beta)
#' @return N-vector of log P(Y_i | X_i, params)
logit_loglik <- function(Y, Xbeta) {
  # log P(Y=1) = Xbeta - log(1 + exp(Xbeta))
  # log P(Y=0) = -log(1 + exp(Xbeta))
  # Numerically stable: use -log(1+exp(x)) = -softplus(x)
  # softplus(x) = x + log(1+exp(-x)) for x > 0, log(1+exp(x)) for x <= 0
  softplus <- ifelse(Xbeta > 0,
                     Xbeta + log1p(exp(-Xbeta)),
                     log1p(exp(Xbeta)))
  Y * Xbeta - softplus
}


# ==============================================================================
# 4. EM Estimation for Binary Logit Mixture
# ==============================================================================

estimate_demand_mixture <- function(data, K, max_iter = 200, tol = 1e-6) {
  N <- nrow(data)
  Y <- data$Y
  price <- data$price
  quality <- data$quality
  treat <- data$treatment

  # Number of parameters: K intercepts + K beta_price + K beta_quality
  #                       + K gamma + (K-1) mixing proportions
  n_params <- 4 * K + (K - 1)

  # Initialize parameters with perturbation
  params <- list()
  for (k in 1:K) {
    params[[k]] <- list(
      alpha      = rnorm(1, 0, 0.3),
      beta_price = runif(1, -0.5, -0.05),
      beta_qual  = runif(1, 0.05, 0.35),
      gamma      = rnorm(1, 0, 0.1)
    )
  }
  pi_k <- rep(1 / K, K)
  tau  <- matrix(1 / K, N, K)
  ll   <- -Inf

  # Pre-build design vectors (constant across iterations)
  # Xbeta_k = alpha_k + beta_price_k * price + beta_qual_k * quality + gamma_k * treat

  for (iter in 1:max_iter) {

    # ---- E-step: compute posterior type probabilities ----
    log_tau <- matrix(0, N, K)
    for (k in 1:K) {
      Xbeta_k <- params[[k]]$alpha +
                 params[[k]]$beta_price * price +
                 params[[k]]$beta_qual  * quality +
                 params[[k]]$gamma      * treat
      log_tau[, k] <- log(max(pi_k[k], 1e-300)) + logit_loglik(Y, Xbeta_k)
    }

    # Log-sum-exp for numerical stability
    log_max <- apply(log_tau, 1, max)
    log_sum <- log_max + log(rowSums(exp(log_tau - log_max)))
    ll_new  <- sum(log_sum)

    # Posterior probabilities
    tau <- exp(log_tau - log_sum)
    tau <- pmax(tau, 1e-300)
    tau <- tau / rowSums(tau)

    # Convergence check
    if (is.finite(ll) && is.finite(ll_new) &&
        abs(ll_new - ll) / (abs(ll) + 1) < tol) break
    ll <- ll_new

    # ---- M-step: update parameters via weighted IRLS ----
    for (k in 1:K) {
      pi_k[k] <- max(mean(tau[, k]), 1e-6)
      w <- tau[, k]  # posterior weights

      # Run a few Newton-like gradient steps (IRLS)
      for (g in 1:5) {
        Xbeta_k <- params[[k]]$alpha +
                   params[[k]]$beta_price * price +
                   params[[k]]$beta_qual  * quality +
                   params[[k]]$gamma      * treat

        # Predicted probability
        p_hat <- 1 / (1 + exp(-Xbeta_k))
        p_hat <- pmin(pmax(p_hat, 1e-10), 1 - 1e-10)

        # Weighted residual
        resid <- w * (Y - p_hat)

        # Weighted Hessian diagonal approximation
        W_diag <- w * p_hat * (1 - p_hat) + 1e-8

        # Gradient and Hessian-based updates (diagonal approx)
        sum_W <- sum(W_diag)

        # Intercept
        params[[k]]$alpha <- params[[k]]$alpha +
          sum(resid) / sum_W

        # Price
        grad_price <- sum(resid * price)
        H_price    <- sum(W_diag * price^2)
        params[[k]]$beta_price <- params[[k]]$beta_price +
          grad_price / max(H_price, 1e-8)

        # Quality
        grad_qual <- sum(resid * quality)
        H_qual    <- sum(W_diag * quality^2)
        params[[k]]$beta_qual <- params[[k]]$beta_qual +
          grad_qual / max(H_qual, 1e-8)

        # Treatment
        grad_gamma <- sum(resid * treat)
        H_gamma    <- sum(W_diag * treat^2)
        params[[k]]$gamma <- params[[k]]$gamma +
          grad_gamma / max(H_gamma, 1e-8)
      }
    }
  }

  # Normalize mixing proportions
  pi_k <- pi_k / sum(pi_k)

  # ---- BIC ----
  bic_val <- -2 * ll + n_params * log(N)

  # ---- Counterfactual: remove treatment (set treatment = 0 for all) ----
  # Baseline: predicted purchase probability under observed data
  prob_base <- rep(0, N)
  for (k in 1:K) {
    Xbeta_k <- params[[k]]$alpha +
               params[[k]]$beta_price * price +
               params[[k]]$beta_qual  * quality +
               params[[k]]$gamma      * treat
    p_k <- 1 / (1 + exp(-Xbeta_k))
    prob_base <- prob_base + tau[, k] * p_k
  }

  # Counterfactual: set treatment = 0 for everyone
  prob_cf <- rep(0, N)
  for (k in 1:K) {
    Xbeta_cf <- params[[k]]$alpha +
                params[[k]]$beta_price * price +
                params[[k]]$beta_qual  * quality
                # gamma * 0 = 0 (no treatment)
    p_cf <- 1 / (1 + exp(-Xbeta_cf))
    prob_cf <- prob_cf + tau[, k] * p_cf
  }

  purchase_base <- mean(prob_base)
  purchase_cf   <- mean(prob_cf)
  cf_effect_pp  <- (purchase_base - purchase_cf) * 100  # percentage points
  cf_effect_pct <- (purchase_base - purchase_cf) / purchase_cf * 100  # percent change

  # ---- Collect gamma estimates ----
  gamma_est <- sapply(params, function(p) p$gamma)
  pi_est    <- pi_k

  list(K = K, loglik = ll, bic = bic_val, n_params = n_params,
       pi = pi_est, gamma = gamma_est,
       gamma_weighted = sum(pi_est * gamma_est),
       cf_effect_pp = cf_effect_pp, cf_effect_pct = cf_effect_pct,
       purchase_base = purchase_base, purchase_cf = purchase_cf,
       converged = (iter < max_iter), n_iter = min(iter, max_iter),
       params = params)
}


# ==============================================================================
# 5. Three-Pronged Selection (adapted for demand model)
# ==============================================================================

three_prong_select_demand <- function(results) {
  K_values  <- sapply(results, function(r) r$K)
  bic_values <- sapply(results, function(r) r$bic)
  cf_values  <- sapply(results, function(r) r$cf_effect_pp)

  # Prong 1: BIC
  k_bic <- K_values[which.min(bic_values)]

  # Prong 2: CF stability (first K where |CF(K) - CF(K-1)| < threshold)
  cf_changes <- abs(diff(cf_values))
  k_stable   <- max(K_values)
  if (length(cf_changes) >= 2) {
    for (idx in 2:length(cf_changes)) {
      if (cf_changes[idx] < 0.5) {  # 0.5 pp threshold
        k_stable <- K_values[idx + 1]
        break
      }
    }
  }

  # Prong 3: OSCE (number of economically significant types)
  r_max  <- results[[length(results)]]
  k_osce <- max(1, sum(abs(r_max$gamma) > 0.01 & r_max$pi > 0.03))

  # Consensus: modal vote
  votes <- c(k_bic, k_stable, k_osce)
  k_consensus <- as.integer(names(sort(table(votes), decreasing = TRUE))[1])

  list(k_bic = k_bic, k_stable = k_stable, k_osce = k_osce,
       k_consensus = k_consensus)
}


# ==============================================================================
# 6. Main Simulation Loop
# ==============================================================================

run_demand_simulation <- function(dgp, N, R, K_max = 5) {
  cat(sprintf("\n=== Demand MC: N = %d, R = %d, K_max = %d ===\n", N, R, K_max))
  all_results <- list()
  idx <- 0
  t0 <- proc.time()

  for (r in 1:R) {
    if (r %% 10 == 0 || r == 1) {
      elapsed <- (proc.time() - t0)[3]
      eta <- if (r > 1) elapsed / (r - 1) * (R - r + 1) else NA
      cat(sprintf("  Rep %d/%d (%.0fs elapsed, ETA %.0fs)\n",
                  r, R, elapsed, ifelse(is.na(eta), 0, eta)))
    }

    data <- simulate_demand_data(N, dgp)

    results_k <- list()
    for (K in 1:K_max) {
      tryCatch({
        results_k[[K]] <- estimate_demand_mixture(data, K)
      }, error = function(e) {
        results_k[[K]] <<- list(
          K = K, loglik = NA, bic = Inf,
          cf_effect_pp = NA, cf_effect_pct = NA,
          gamma = rep(NA, K), pi = rep(1/K, K),
          gamma_weighted = NA, purchase_base = NA, purchase_cf = NA,
          converged = FALSE, n_iter = NA, n_params = NA
        )
      })
    }

    # Three-pronged selection
    valid <- results_k[!sapply(results_k, function(r)
      is.null(r) || is.na(r$loglik))]
    sel <- if (length(valid) >= 2) {
      three_prong_select_demand(valid)
    } else {
      list(k_bic = NA, k_stable = NA, k_osce = NA, k_consensus = NA)
    }

    # Collect results for each K
    for (K in 1:K_max) {
      rk <- results_k[[K]]
      if (!is.null(rk) && !is.na(rk$loglik)) {
        idx <- idx + 1
        all_results[[idx]] <- data.frame(
          N = N, rep = r, K_est = K,
          loglik = rk$loglik, bic = rk$bic, n_params = rk$n_params,
          cf_effect_pp  = rk$cf_effect_pp,
          cf_effect_pct = rk$cf_effect_pct,
          purchase_base = rk$purchase_base,
          purchase_cf   = rk$purchase_cf,
          gamma_weighted = rk$gamma_weighted,
          converged = rk$converged,
          k_selected_bic    = sel$k_bic,
          k_selected_stable = sel$k_stable,
          k_selected_osce   = sel$k_osce,
          k_selected_consensus = sel$k_consensus,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  elapsed_total <- (proc.time() - t0)[3]
  cat(sprintf("  Complete: %.0f sec (%.1f min)\n", elapsed_total, elapsed_total / 60))
  bind_rows(all_results)
}


# ==============================================================================
# 7. Run Simulations
# ==============================================================================

N_sim <- 10000
R_sim <- 100
K_max <- 5

cat("\n============================================================\n")
cat(sprintf("RUNNING: N = %d, R = %d, K = 1,...,%d\n", N_sim, R_sim, K_max))
cat("============================================================\n")

t_start <- proc.time()
results_all <- run_demand_simulation(dgp_demand, N = N_sim, R = R_sim, K_max = K_max)
t_total <- (proc.time() - t_start)[3]
cat(sprintf("\nTotal simulation time: %.0f sec (%.1f min)\n", t_total, t_total / 60))


# ==============================================================================
# 8. Summary Statistics
# ==============================================================================

cat("\n============================================================\n")
cat("SUMMARY STATISTICS\n")
cat("============================================================\n")

# K selection frequency
k_freq <- results_all %>%
  filter(K_est == k_selected_consensus) %>%
  distinct(rep, .keep_all = TRUE) %>%
  group_by(k_selected_consensus) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(freq = count / sum(count))

cat("\nK Selection Frequency (Consensus, True K=3):\n")
print(as.data.frame(k_freq))

k_recovery_rate <- sum(k_freq$count[k_freq$k_selected_consensus == 3]) /
                   sum(k_freq$count) * 100
cat(sprintf("\nK=3 recovery rate: %.1f%%\n", k_recovery_rate))

# BIC-based selection
k_bic_freq <- results_all %>%
  filter(K_est == k_selected_bic) %>%
  distinct(rep, .keep_all = TRUE) %>%
  group_by(k_selected_bic) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(freq = count / sum(count))

cat("\nK Selection by BIC:\n")
print(as.data.frame(k_bic_freq))

# CF statistics by K
cf_stats <- results_all %>%
  group_by(K_est) %>%
  summarise(
    cf_pp_mean = mean(cf_effect_pp, na.rm = TRUE),
    cf_pp_sd   = sd(cf_effect_pp, na.rm = TRUE),
    cf_pp_q25  = quantile(cf_effect_pp, 0.25, na.rm = TRUE),
    cf_pp_q75  = quantile(cf_effect_pp, 0.75, na.rm = TRUE),
    cf_pct_mean = mean(cf_effect_pct, na.rm = TRUE),
    gamma_wt_mean = mean(gamma_weighted, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nCounterfactual by K (percentage points):\n")
print(as.data.frame(cf_stats))

# CF range across K for each replication
cf_range <- results_all %>%
  group_by(rep) %>%
  summarise(
    cf_range_pp = max(cf_effect_pp, na.rm = TRUE) - min(cf_effect_pp, na.rm = TRUE),
    cf_min = min(cf_effect_pp, na.rm = TRUE),
    cf_max = max(cf_effect_pp, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nCounterfactual Range Across K (pp):\n")
cat(sprintf("  Mean range:   %.2f pp\n", mean(cf_range$cf_range_pp, na.rm = TRUE)))
cat(sprintf("  Median range: %.2f pp\n", median(cf_range$cf_range_pp, na.rm = TRUE)))
cat(sprintf("  Min range:    %.2f pp\n", min(cf_range$cf_range_pp, na.rm = TRUE)))
cat(sprintf("  Max range:    %.2f pp\n", max(cf_range$cf_range_pp, na.rm = TRUE)))

# Gamma recovery
gamma_stats <- results_all %>%
  filter(K_est == 3) %>%
  summarise(
    gamma_wt_mean = mean(gamma_weighted, na.rm = TRUE),
    gamma_wt_sd   = sd(gamma_weighted, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("\nGamma recovery (K=3): weighted avg = %.4f (sd = %.4f), true = %.4f\n",
            gamma_stats$gamma_wt_mean, gamma_stats$gamma_wt_sd, weighted_avg))


# ==============================================================================
# 9. Build Summary Table
# ==============================================================================

summary_df <- data.frame(
  metric = c("N", "R", "K_true", "K_recovery_rate_pct",
             "weighted_gamma_true",
             paste0("cf_pp_K", 1:K_max),
             "cf_range_mean_pp", "cf_range_median_pp"),
  value = c(N_sim, R_sim, dgp_demand$K_true, k_recovery_rate,
            weighted_avg,
            cf_stats$cf_pp_mean,
            mean(cf_range$cf_range_pp, na.rm = TRUE),
            median(cf_range$cf_range_pp, na.rm = TRUE)),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 10. Save Results
# ==============================================================================

write.csv(results_all,
          file.path(DATA_DIR, "output/monte_carlo_demand_results.csv"),
          row.names = FALSE)

write.csv(summary_df,
          file.path(DATA_DIR, "output/monte_carlo_demand_summary.csv"),
          row.names = FALSE)

cat("\nCSV results saved:\n")
cat(sprintf("  %s\n", file.path(DATA_DIR, "output/monte_carlo_demand_results.csv")))
cat(sprintf("  %s\n", file.path(DATA_DIR, "output/monte_carlo_demand_summary.csv")))


# ==============================================================================
# 11. Generate Figure: fig_mc7_demand_cf
# ==============================================================================

cat("\nGENERATING FIGURE\n")

tryCatch({
  # Panel A: CF sensitivity by K (main point)
  cf_by_k <- results_all %>%
    group_by(K_est) %>%
    summarise(
      cf_mean = mean(cf_effect_pp, na.rm = TRUE),
      cf_q10  = quantile(cf_effect_pp, 0.10, na.rm = TRUE),
      cf_q25  = quantile(cf_effect_pp, 0.25, na.rm = TRUE),
      cf_q75  = quantile(cf_effect_pp, 0.75, na.rm = TRUE),
      cf_q90  = quantile(cf_effect_pp, 0.90, na.rm = TRUE),
      .groups = "drop"
    )

  fig_mc7 <- ggplot(cf_by_k, aes(x = K_est)) +
    # 80% interval (10th-90th percentile)
    geom_ribbon(aes(ymin = cf_q10, ymax = cf_q90),
                fill = "#4575b4", alpha = 0.15) +
    # IQR (25th-75th percentile)
    geom_ribbon(aes(ymin = cf_q25, ymax = cf_q75),
                fill = "#4575b4", alpha = 0.3) +
    # Mean line
    geom_line(aes(y = cf_mean), color = "#2166ac", size = 1.3) +
    geom_point(aes(y = cf_mean), color = "#2166ac", size = 3.5) +
    # Reference line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.5) +
    # True K indicator
    geom_vline(xintercept = 3, linetype = "dotted", color = "#d73027",
               size = 0.7, alpha = 0.7) +
    annotate("text", x = 3.15, y = max(cf_by_k$cf_q90) * 0.95,
             label = "True K=3", color = "#d73027", hjust = 0, size = 3.5) +
    scale_x_continuous(breaks = 1:K_max) +
    labs(
      title = "Counterfactual Sensitivity to K: Binary Demand Model",
      subtitle = paste0("Treatment removal effect (pp), N=",
                        format(N_sim, big.mark = ","), ", R=", R_sim,
                        " | Shading: IQR (dark), 80% interval (light)"),
      x = "Number of Consumer Types (K)",
      y = "Treatment Effect on Purchase Rate (pp)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(OUTPUT_DIR, "fig_mc7_demand_cf.pdf"),
         fig_mc7, width = 8, height = 5.5)
  ggsave(file.path(OUTPUT_DIR, "fig_mc7_demand_cf.png"),
         fig_mc7, width = 8, height = 5.5, dpi = 300)

  cat("  Saved: fig_mc7_demand_cf.pdf / .png\n")

}, error = function(e) {
  cat(sprintf("  fig_mc7 error: %s\n", e$message))
})


# ==============================================================================
# 12. Final Summary
# ==============================================================================

cat("\n============================================================\n")
cat("STYLIZED DEMAND MC COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Total time: %.1f min\n", t_total / 60))
cat(sprintf("K=3 recovery rate: %.1f%%\n", k_recovery_rate))
cat(sprintf("CF range across K: %.2f pp (median)\n",
            median(cf_range$cf_range_pp, na.rm = TRUE)))
cat(sprintf("Key finding: even in a simple binary demand model,\n"))
cat(sprintf("  K choice drives substantial CF variation (%.2f pp range),\n",
            mean(cf_range$cf_range_pp, na.rm = TRUE)))
cat(sprintf("  confirming K sensitivity generalizes beyond banking/MNL.\n"))
cat("============================================================\n")
