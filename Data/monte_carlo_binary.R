#' ==============================================================================
#' Monte Carlo Simulation: Binary Logit Mixture (Comment 3c)
#'
#' Binary outcome version of the main MC simulation.
#' Instead of a 9-alternative MNL, this uses a binary logit mixture
#' where the outcome is self-employed (1) vs not (0).
#'
#' Two DGP designs (same gamma vectors as main MC):
#'   Design 1 (Cancellation): Type-specific branch effects cancel in aggregate
#'   Design 2 (Same-Sign): All type branch effects are positive
#'
#' Vectorized EM with gradient-descent M-step. Runtime: ~1-2 hours.
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
cat("MONTE CARLO SIMULATION: BINARY LOGIT MIXTURE (Comment 3c)\n")
cat("Vectorized EM implementation\n")
cat("============================================================\n\n")


# ==============================================================================
# 1. Constants and DGP Parameters
# ==============================================================================

# Common parameters across both designs
K_TRUE <- 4
PI_TRUE <- c(0.13, 0.20, 0.35, 0.32)
BETA_0_TRUE <- c(-1.5, -1.2, -1.0, -0.8)   # type-specific intercepts
BETA_1_TRUE <- c(0.3, 0.2, 0.1, -0.1)       # type-specific slopes on x1
BETA_2_TRUE <- c(0.4, 0.3, 0.2, 0.1)        # type-specific slopes on x2

dgp_cancellation <- list(
  name = "Cancellation", K_true = K_TRUE,
  pi = PI_TRUE,
  gamma = c(0.002, -0.030, 0.026, 0.138),
  beta_0 = BETA_0_TRUE,
  beta_1 = BETA_1_TRUE,
  beta_2 = BETA_2_TRUE
)

dgp_samesign <- list(
  name = "Same-Sign", K_true = K_TRUE,
  pi = PI_TRUE,
  gamma = c(0.05, 0.08, 0.12, 0.15),
  beta_0 = BETA_0_TRUE,
  beta_1 = BETA_1_TRUE,
  beta_2 = BETA_2_TRUE
)

for (dgp in list(dgp_cancellation, dgp_samesign)) {
  cat(sprintf("Design: %s\n", dgp$name))
  cat(sprintf("  gamma = (%s), weighted avg = %.4f\n",
              paste(dgp$gamma, collapse = ", "), sum(dgp$pi * dgp$gamma)))
  cat(sprintf("  beta_0 = (%s)\n", paste(dgp$beta_0, collapse = ", ")))
}
cat("\n")


# ==============================================================================
# 2. Vectorized DGP: Binary Logit Mixture
# ==============================================================================

#' Simulate binary SE outcome from a K-component logit mixture
#'
#' P(SE=1 | type=k, x) = logistic(beta_0^k + beta_1^k * x1 + beta_2^k * x2
#'                                  + gamma^k * branch)
simulate_binary_data <- function(N, dgp, n_cbsa = 50) {
  K <- dgp$K_true

  # Draw latent types
  type <- sample(1:K, N, replace = TRUE, prob = dgp$pi)

  # Covariates
  cbsa <- sample(1:n_cbsa, N, replace = TRUE)
  x1 <- rnorm(N, 0, 1)
  x2 <- rbinom(N, 1, 0.4)
  branch <- rbinom(N, 1, 0.5)  # branch indicator

  # Compute latent utility for SE (vectorized over types)
  u <- dgp$beta_0[type] +
       dgp$beta_1[type] * x1 +
       dgp$beta_2[type] * x2 +
       dgp$gamma[type] * branch

  # Binary outcome via logistic
  prob_se <- 1 / (1 + exp(-u))
  y <- rbinom(N, 1, prob_se)

  data.frame(y = y, x1 = x1, x2 = x2, branch = branch, cbsa = cbsa)
}


# ==============================================================================
# 3. Logistic Helper
# ==============================================================================

logistic <- function(z) {
  # Numerically stable logistic function
  pos <- z >= 0
  out <- numeric(length(z))
  out[pos] <- 1 / (1 + exp(-z[pos]))
  ez <- exp(z[!pos])
  out[!pos] <- ez / (1 + ez)
  out
}


# ==============================================================================
# 4. Binary EM Estimation
# ==============================================================================

#' Compute P(SE=1 | type k, data) for all N individuals
#' Returns N-vector of probabilities
compute_prob_k <- function(x1, x2, branch, pk) {
  u <- pk$beta_0 + pk$beta_1 * x1 + pk$beta_2 * x2 + pk$gamma * branch
  logistic(u)
}

#' Log-likelihood contribution for individual i under type k
#' l_ik = y_i * log(p_ik) + (1 - y_i) * log(1 - p_ik)
log_lik_k <- function(y, prob_k) {
  prob_k <- pmax(pmin(prob_k, 1 - 1e-15), 1e-15)
  y * log(prob_k) + (1 - y) * log(1 - prob_k)
}

#' Estimate K-component binary logit mixture via EM
#'
#' E-step: posterior type responsibilities
#' M-step: gradient descent on type-specific logistic regression parameters
estimate_binary_mixture <- function(data, K, max_iter = 150, tol = 1e-6) {
  N <- nrow(data)
  y <- data$y
  x1 <- data$x1
  x2 <- data$x2
  branch <- data$branch

  # Initialize parameters with small random perturbations
  params <- list()
  pi_k <- rep(1 / K, K)
  for (k in 1:K) {
    params[[k]] <- list(
      beta_0 = rnorm(1, -1.0, 0.3),
      beta_1 = rnorm(1, 0, 0.1),
      beta_2 = rnorm(1, 0, 0.1),
      gamma  = rnorm(1, 0, 0.02)
    )
  }

  # Responsibility matrix N x K
  tau <- matrix(1 / K, N, K)
  ll <- -Inf

  for (iter in 1:max_iter) {
    # ---- E-step ----
    # Compute log P(y_i | type=k, x_i) for each k, then weight by pi_k
    log_tau <- matrix(0, N, K)
    for (k in 1:K) {
      prob_k <- compute_prob_k(x1, x2, branch, params[[k]])
      log_tau[, k] <- log(pi_k[k]) + log_lik_k(y, prob_k)
    }

    # Log-sum-exp for numerical stability
    log_max <- apply(log_tau, 1, max)
    log_sum <- log_max + log(rowSums(exp(log_tau - log_max)))
    ll_new <- sum(log_sum)

    # Posterior responsibilities
    tau <- exp(log_tau - log_sum)

    # Check convergence
    if (is.finite(ll) && is.finite(ll_new) &&
        abs(ll_new - ll) / (abs(ll) + 1) < tol) break
    ll <- ll_new

    # ---- M-step ----
    for (k in 1:K) {
      # Update mixing proportions
      pi_k[k] <- max(mean(tau[, k]), 1e-6)

      w <- tau[, k]  # N-vector of responsibilities
      lr <- 0.05 / (1 + iter * 0.03)  # decaying learning rate

      # Multiple gradient steps per M-step for faster convergence
      for (g in 1:3) {
        prob_k <- compute_prob_k(x1, x2, branch, params[[k]])

        # Weighted residuals: r_i = w_i * (y_i - p_ik)
        resid <- w * (y - prob_k)

        # Gradient updates for logistic regression parameters
        params[[k]]$beta_0 <- params[[k]]$beta_0 + lr * sum(resid) / N
        params[[k]]$beta_1 <- params[[k]]$beta_1 + lr * sum(resid * x1) / N
        params[[k]]$beta_2 <- params[[k]]$beta_2 + lr * sum(resid * x2) / N
        params[[k]]$gamma  <- params[[k]]$gamma  + lr * sum(resid * branch) / N
      }
    }

    # Renormalize mixing proportions
    pi_k <- pi_k / sum(pi_k)
  }

  # Final log-likelihood
  log_tau_final <- matrix(0, N, K)
  for (k in 1:K) {
    prob_k <- compute_prob_k(x1, x2, branch, params[[k]])
    log_tau_final[, k] <- log(pi_k[k]) + log_lik_k(y, prob_k)
  }
  log_max <- apply(log_tau_final, 1, max)
  ll_final <- sum(log_max + log(rowSums(exp(log_tau_final - log_max))))

  # BIC computation
  # Parameters: K intercepts + K slopes on x1 + K slopes on x2 + K gammas + (K-1) mixing
  np <- K * 4 + (K - 1)
  n_cbsa <- length(unique(data$cbsa))
  bic_std <- -2 * ll_final + np * log(N)
  c_T <- 1 + log(N / n_cbsa) / log(N)
  bic_panel <- -2 * ll_final + np * c_T * log(n_cbsa)

  # ---- Counterfactual: 50% reduction in branch indicator ----
  # Baseline SE rate
  se_base <- mean(y)

  # Counterfactual: reduce branch by 50% (multiply branch by 0.5 in utility)
  prob_se_cf <- rep(0, N)
  for (k in 1:K) {
    pk_cf <- params[[k]]
    # 50% closure = halve the branch effect (equivalent to removing half of branches)
    u_cf <- pk_cf$beta_0 + pk_cf$beta_1 * x1 + pk_cf$beta_2 * x2 +
            pk_cf$gamma * branch * 0.5
    p_cf_k <- logistic(u_cf)
    prob_se_cf <- prob_se_cf + tau[, k] * p_cf_k
  }
  se_cf <- mean(prob_se_cf)
  cf_effect <- (se_cf - se_base) / se_base * 100

  list(
    K = K,
    loglik = ll_final,
    bic_standard = bic_std,
    bic_panel = bic_panel,
    n_params = np,
    pi = pi_k,
    gamma = sapply(params, function(p) p$gamma),
    beta_0 = sapply(params, function(p) p$beta_0),
    cf_effect = cf_effect,
    se_rate_baseline = se_base,
    se_rate_cf = se_cf,
    converged = iter < max_iter
  )
}


# ==============================================================================
# 5. Three-Pronged Selection (same logic as main MC)
# ==============================================================================

three_prong_select <- function(results) {
  K_values <- sapply(results, function(r) r$K)
  bic_values <- sapply(results, function(r) r$bic_panel)
  k_bic <- K_values[which.min(bic_values)]

  cf_values <- sapply(results, function(r) r$cf_effect)
  cf_changes <- abs(diff(cf_values))
  k_stable <- max(K_values)
  if (length(cf_changes) >= 2) {
    for (idx in 2:length(cf_changes)) {
      if (cf_changes[idx] < 2.0) { k_stable <- K_values[idx + 1]; break }
    }
  }

  r_max <- results[[length(results)]]
  k_osce <- max(1, sum(abs(r_max$gamma) > 0.005 & r_max$pi > 0.02))

  votes <- c(k_bic, k_stable, k_osce)
  k_consensus <- as.integer(names(sort(table(votes), decreasing = TRUE))[1])
  list(k_bic = k_bic, k_stable = k_stable, k_osce = k_osce,
       k_consensus = k_consensus)
}


# ==============================================================================
# 6. Main Simulation Loop
# ==============================================================================

run_binary_simulation <- function(dgp, N, R, K_max = 5) {
  cat(sprintf("\n=== Design: %s (N=%d, R=%d) ===\n", dgp$name, N, R))
  all_results <- list()
  idx <- 0
  t0 <- proc.time()

  for (r in 1:R) {
    if (r %% 10 == 0 || r == 1) {
      elapsed <- (proc.time() - t0)[3]
      eta <- if (r > 1) elapsed / (r - 1) * (R - r + 1) else NA
      cat(sprintf("  Rep %d/%d (%.0fs elapsed, ETA %.0fs)\n", r, R, elapsed,
                  ifelse(is.na(eta), 0, eta)))
    }

    data <- simulate_binary_data(N, dgp)
    results_k <- list()

    for (K in 1:K_max) {
      tryCatch({
        results_k[[K]] <- estimate_binary_mixture(data, K)
      }, error = function(e) {
        results_k[[K]] <<- list(
          K = K, loglik = NA, bic_standard = Inf, bic_panel = Inf,
          cf_effect = NA, gamma = rep(NA, K), pi = rep(1/K, K),
          beta_0 = rep(NA, K), converged = FALSE,
          se_rate_baseline = NA, se_rate_cf = NA
        )
      })
    }

    # Three-pronged selection
    valid <- results_k[!sapply(results_k, function(r)
      is.null(r) || is.na(r$loglik))]
    sel <- if (length(valid) >= 2) three_prong_select(valid) else
      list(k_bic = NA, k_stable = NA, k_osce = NA, k_consensus = NA)

    for (K in 1:K_max) {
      rk <- results_k[[K]]
      if (!is.null(rk) && !is.na(rk$loglik)) {
        idx <- idx + 1
        all_results[[idx]] <- data.frame(
          design = dgp$name, N = N, rep = r, K_est = K,
          loglik = rk$loglik,
          bic_standard = rk$bic_standard,
          bic_panel = rk$bic_panel,
          cf_effect = rk$cf_effect,
          se_rate_baseline = rk$se_rate_baseline,
          se_rate_cf = rk$se_rate_cf,
          gamma_weighted = sum(rk$pi * rk$gamma),
          k_selected_bic = sel$k_bic,
          k_selected_stable = sel$k_stable,
          k_selected_osce = sel$k_osce,
          k_selected_consensus = sel$k_consensus,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  elapsed_total <- (proc.time() - t0)[3]
  cat(sprintf("  Design '%s' complete: %.0f sec (%.1f min)\n",
              dgp$name, elapsed_total, elapsed_total / 60))

  bind_rows(all_results)
}


# ==============================================================================
# 7. Run Simulations: R=100 at N=10,000 for both designs
# ==============================================================================

N_SIM <- 10000
R_SIM <- 100
K_MAX <- 5

cat("\n============================================================\n")
cat("RUNNING BINARY LOGIT MIXTURE SIMULATIONS\n")
cat(sprintf("  N = %d, R = %d, K_max = %d\n", N_SIM, R_SIM, K_MAX))
cat("============================================================\n")

t_start <- proc.time()

results_cancel <- run_binary_simulation(dgp_cancellation, N_SIM, R_SIM, K_MAX)
results_same   <- run_binary_simulation(dgp_samesign, N_SIM, R_SIM, K_MAX)
results_all    <- bind_rows(results_cancel, results_same)

t_total <- (proc.time() - t_start)[3]
cat(sprintf("\nTotal simulation time: %.0f sec (%.1f hours)\n",
            t_total, t_total / 3600))


# ==============================================================================
# 8. Summary Statistics
# ==============================================================================

cat("\n============================================================\n")
cat("SUMMARY STATISTICS\n")
cat("============================================================\n")

# K selection frequency (consensus)
k_freq <- results_all %>%
  filter(K_est == k_selected_consensus) %>%
  distinct(design, N, rep, .keep_all = TRUE) %>%
  group_by(design, N, k_selected_consensus) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(design, N) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

cat("\nK Selection Frequency (Consensus):\n")
print(as.data.frame(k_freq))

# CF statistics by design and K
cf_stats <- results_all %>%
  group_by(design, N, K_est) %>%
  summarise(
    cf_mean = mean(cf_effect, na.rm = TRUE),
    cf_sd   = sd(cf_effect, na.rm = TRUE),
    cf_q25  = quantile(cf_effect, 0.25, na.rm = TRUE),
    cf_q75  = quantile(cf_effect, 0.75, na.rm = TRUE),
    se_rate_mean = mean(se_rate_baseline, na.rm = TRUE),
    gamma_wt_mean = mean(gamma_weighted, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nCounterfactual Statistics by K:\n")
print(as.data.frame(cf_stats))

# CF range across K within each replication
cf_range <- results_all %>%
  group_by(design, N, rep) %>%
  summarise(
    cf_range = max(cf_effect, na.rm = TRUE) - min(cf_effect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(design, N) %>%
  summarise(
    mean_range   = mean(cf_range, na.rm = TRUE),
    median_range = median(cf_range, na.rm = TRUE),
    sd_range     = sd(cf_range, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nCounterfactual Range Across K:\n")
print(as.data.frame(cf_range))

# BIC accuracy
bic_acc <- results_all %>%
  filter(K_est == k_selected_bic) %>%
  distinct(design, N, rep, .keep_all = TRUE) %>%
  mutate(correct = k_selected_bic == 4) %>%
  group_by(design, N) %>%
  summarise(accuracy = mean(correct, na.rm = TRUE), .groups = "drop")

cat("\nBIC Accuracy (K=4 recovery):\n")
print(as.data.frame(bic_acc))


# ==============================================================================
# 9. Save Results
# ==============================================================================

write.csv(results_all,
          file.path(DATA_DIR, "output/monte_carlo_binary_results.csv"),
          row.names = FALSE)

# Create summary table combining cf_stats and cf_range
summary_combined <- cf_stats %>%
  left_join(
    cf_range %>% select(design, N, mean_range, median_range),
    by = c("design", "N")
  )

write.csv(summary_combined,
          file.path(DATA_DIR, "output/monte_carlo_binary_summary.csv"),
          row.names = FALSE)

cat("\nCSV results saved:\n")
cat("  output/monte_carlo_binary_results.csv (all replications)\n")
cat("  output/monte_carlo_binary_summary.csv (summary stats by design x K)\n")


# ==============================================================================
# 10. Figure: fig_mc6_binary_cf
# ==============================================================================

cat("\nGENERATING FIGURE\n")

tryCatch({
  # Counterfactual effect by K, for both designs
  cf_plot_data <- results_all %>%
    group_by(design, K_est) %>%
    summarise(
      cf_mean = mean(cf_effect, na.rm = TRUE),
      cf_q25  = quantile(cf_effect, 0.25, na.rm = TRUE),
      cf_q75  = quantile(cf_effect, 0.75, na.rm = TRUE),
      cf_q10  = quantile(cf_effect, 0.10, na.rm = TRUE),
      cf_q90  = quantile(cf_effect, 0.90, na.rm = TRUE),
      .groups = "drop"
    )

  fig_mc6 <- ggplot(cf_plot_data, aes(x = K_est, y = cf_mean, color = design)) +
    # 10-90 percentile ribbon
    geom_ribbon(aes(ymin = cf_q10, ymax = cf_q90, fill = design),
                alpha = 0.10, color = NA) +
    # IQR ribbon
    geom_ribbon(aes(ymin = cf_q25, ymax = cf_q75, fill = design),
                alpha = 0.25, color = NA) +
    # Mean line and points
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    # Reference line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    # Vertical line at true K
    geom_vline(xintercept = 4, linetype = "dotted", color = "gray40") +
    annotate("text", x = 4.15, y = max(cf_plot_data$cf_q90) * 0.95,
             label = "True K=4", hjust = 0, size = 3.5, color = "gray30") +
    # Styling
    scale_color_manual(
      values = c("Cancellation" = "#d73027", "Same-Sign" = "#4575b4"),
      name = "DGP Design"
    ) +
    scale_fill_manual(
      values = c("Cancellation" = "#d73027", "Same-Sign" = "#4575b4"),
      name = "DGP Design"
    ) +
    scale_x_continuous(breaks = 1:K_MAX) +
    labs(
      title = "Binary Logit Mixture: Counterfactual Sensitivity to K",
      subtitle = sprintf("50%% branch reduction, N=%s, R=%d replications",
                         format(N_SIM, big.mark = ","), R_SIM),
      x = "Number of Mixture Components (K)",
      y = "Counterfactual Effect on SE Rate (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(OUTPUT_DIR, "fig_mc6_binary_cf.pdf"),
         fig_mc6, width = 8, height = 5.5)
  ggsave(file.path(OUTPUT_DIR, "fig_mc6_binary_cf.png"),
         fig_mc6, width = 8, height = 5.5, dpi = 300)

  cat("  Saved: fig_mc6_binary_cf.pdf / .png\n")

}, error = function(e) cat(sprintf("  fig_mc6 error: %s\n", e$message)))


# ==============================================================================
# Done
# ==============================================================================

cat("\n============================================================\n")
cat("BINARY LOGIT MIXTURE MONTE CARLO COMPLETE\n")
cat(sprintf("Total time: %.1f minutes (%.2f hours)\n",
            t_total / 60, t_total / 3600))
cat("============================================================\n")
