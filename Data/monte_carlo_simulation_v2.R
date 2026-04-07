#' ==============================================================================
#' Monte Carlo Simulation v2: When K Matters — R&R Version
#'
#' Revision with increased replications (doubled R values) and Monte Carlo
#' standard error reporting, prepared for the Revise & Resubmit.
#'
#' Changes from v1:
#'   - R values doubled: c(200, 100, 20) instead of c(100, 50, 10)
#'   - MC standard errors (mc_se) added to cf_stats output
#'   - All output files use _v2 suffix
#'   - All figure files use _v2 suffix
#'
#' Two DGP designs to validate the three-pronged selection framework:
#'   Design 1 (Cancellation): Type-specific effects cancel in the aggregate
#'   Design 2 (Same-Sign): All type effects are positive (no cancellation)
#'
#' Vectorized EM implementation. Runtime: ~6-14 hours depending on hardware.
#'
#' Output: CSV results + ggplot figures to Data/output/figures_methods/
#'
#' Author: Alina Malkova
#' Date: March 2026 (v2: April 2026, R&R)
#' ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(20260320)

# Paths
DATA_DIR <- "/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/Data"
OUTPUT_DIR <- file.path(DATA_DIR, "output/figures_methods")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DATA_DIR, "output"), showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("MONTE CARLO SIMULATION v2: WHEN K MATTERS (R&R VERSION)\n")
cat("Vectorized EM implementation — doubled replications\n")
cat("============================================================\n\n")


# ==============================================================================
# 1. Constants and DGP Parameters
# ==============================================================================

J <- 9
alt_b <- ((1:J) - 1) %/% 3 + 1   # Banking mode for each alternative
alt_d <- ((1:J) - 1) %% 3 + 1    # Employment for each alternative
is_branch_se <- (alt_d == 2) & (alt_b == 3)  # Alt 8: Branch x SE
is_se <- (alt_d == 2)                         # Alts 2, 5, 8: SE
se_branch_col <- which(is_branch_se)          # Column index for branch-SE

# Pre-compute base utility offsets for each alternative
alt_b_offset <- (alt_b - 2) * 0.3
alt_d_offset <- (alt_d - 2) * (-0.5)
alt_base_offset <- alt_b_offset + alt_d_offset  # J-vector

dgp_cancellation <- list(
  name = "Cancellation", K_true = 4,
  pi = c(0.13, 0.20, 0.35, 0.32),
  gamma = c(0.002, -0.030, 0.026, 0.138),
  beta_intercept = c(-1.5, -1.2, -1.0, -0.8),
  sigma_x = 1.0
)

dgp_samesign <- list(
  name = "Same-Sign", K_true = 4,
  pi = c(0.13, 0.20, 0.35, 0.32),
  gamma = c(0.05, 0.08, 0.12, 0.15),
  beta_intercept = c(-1.5, -1.2, -1.0, -0.8),
  sigma_x = 1.0
)

for (dgp in list(dgp_cancellation, dgp_samesign)) {
  cat(sprintf("Design: %s\n", dgp$name))
  cat(sprintf("  gamma = (%s), weighted avg = %.4f\n",
              paste(dgp$gamma, collapse=", "), sum(dgp$pi * dgp$gamma)))
}


# ==============================================================================
# 2. Vectorized DGP
# ==============================================================================

simulate_data <- function(N, dgp, n_cbsa = 50) {
  K <- dgp$K_true
  type <- sample(1:K, N, replace = TRUE, prob = dgp$pi)
  cbsa <- sample(1:n_cbsa, N, replace = TRUE)
  x1 <- rnorm(N, 0, dgp$sigma_x)
  x2 <- rbinom(N, 1, 0.4)

  # Build N x J utility matrix using vectorized operations
  # u[i,j] = beta_intercept[type_i] + alt_base_offset[j] + 0.1*x1_i + 0.2*x2_i
  #         + gamma[type_i] * is_branch_se[j] + eps[i,j]
  xb <- 0.1 * x1 + 0.2 * x2  # N-vector
  base_by_type <- dgp$beta_intercept[type]  # N-vector

  # N x J matrix: each row is base_by_type[i] + alt_base_offset + xb[i]
  utils <- outer(base_by_type + xb, alt_base_offset, "+")

  # Add gamma for branch-SE alternative
  utils[, se_branch_col] <- utils[, se_branch_col] + dgp$gamma[type]

  # Gumbel noise
  eps <- -log(-log(matrix(runif(N * J), N, J)))
  choice <- max.col(utils + eps, ties.method = "first")

  data.frame(choice = choice, x1 = x1, x2 = x2, cbsa = cbsa)
}


# ==============================================================================
# 3. Vectorized Softmax
# ==============================================================================

softmax_mat <- function(U) {
  # Numerically stable row-wise softmax using vectorized pmax
  rm <- do.call(pmax, as.data.frame(U))  # row max, much faster than apply
  e <- exp(U - rm)
  rs <- rowSums(e)
  rs[rs == 0] <- 1e-300
  e / rs
}


# ==============================================================================
# 4. Vectorized EM Estimation
# ==============================================================================

#' Compute N x J utility matrix for all individuals under type k
#' Using broadcasting: outer product + column adjustment
compute_utils_k <- function(x1, x2, N, pk) {
  # alpha is J-vector; beta1, beta2, gamma are scalars
  # u[i,j] = alpha[j] + beta1*x1[i] + beta2*x2[i] + gamma*is_branch_se[j]
  xb <- pk$beta1 * x1 + pk$beta2 * x2  # N-vector
  u <- outer(xb, pk$alpha, "+")  # N x J: xb[i] + alpha[j]
  u[, se_branch_col] <- u[, se_branch_col] + pk$gamma
  u
}

estimate_mixture <- function(data, K, max_iter = 80, tol = 1e-5) {
  N <- nrow(data)
  x1 <- data$x1; x2 <- data$x2; ch <- data$choice

  # One-hot
  Y <- matrix(0, N, J); Y[cbind(1:N, ch)] <- 1

  # Initialize
  params <- list(); pi_k <- rep(1/K, K)
  for (k in 1:K) {
    params[[k]] <- list(
      alpha = rnorm(J, 0, 0.2), beta1 = 0, beta2 = 0,
      gamma = rnorm(1, 0, 0.02)
    )
  }
  tau <- matrix(1/K, N, K)
  ll <- -Inf

  for (iter in 1:max_iter) {
    # ---- E-step ----
    for (k in 1:K) {
      Pk <- softmax_mat(compute_utils_k(x1, x2, N, params[[k]]))
      tau[, k] <- pi_k[k] * pmax(Pk[cbind(1:N, ch)], 1e-300)
    }
    rs <- pmax(rowSums(tau), 1e-300)
    ll_new <- sum(log(rs))
    tau <- tau / rs

    if (is.finite(ll) && is.finite(ll_new) &&
        abs(ll_new - ll) / (abs(ll) + 1) < tol) break
    ll <- ll_new

    # ---- M-step ----
    for (k in 1:K) {
      pi_k[k] <- max(mean(tau[, k]), 1e-6)
      w <- tau[, k]
      lr <- 0.02 / (1 + iter * 0.02)

      # 2 gradient descent steps per M-step
      for (g in 1:2) {
        Pk <- softmax_mat(compute_utils_k(x1, x2, N, params[[k]]))
        wr <- (Y - Pk) * w  # N x J weighted residuals

        params[[k]]$alpha <- params[[k]]$alpha + lr * colSums(wr) / N
        params[[k]]$beta1 <- params[[k]]$beta1 + lr * sum(rowSums(wr) * x1) / N
        params[[k]]$beta2 <- params[[k]]$beta2 + lr * sum(rowSums(wr) * x2) / N
        params[[k]]$gamma <- params[[k]]$gamma + lr * sum(wr[, se_branch_col]) / N
      }
    }
  }

  # BIC
  np <- K * (J + 3) + (K - 1)
  n_cbsa <- length(unique(data$cbsa))
  bic_std <- -2 * ll + np * log(N)
  c_T <- 1 + log(N / n_cbsa) / log(N)
  bic_panel <- -2 * ll + np * c_T * log(n_cbsa)

  # Counterfactual: 50% branch closure (reduce branch-SE utility)
  se_base <- mean(ch %in% which(is_se))
  prob_se_cf <- rep(0, N)
  for (k in 1:K) {
    u_cf <- compute_utils_k(x1, x2, N, params[[k]])
    u_cf[, se_branch_col] <- u_cf[, se_branch_col] - 0.5 * params[[k]]$gamma
    P_cf <- softmax_mat(u_cf)
    prob_se_cf <- prob_se_cf + tau[, k] * rowSums(P_cf[, is_se, drop = FALSE])
  }
  se_cf <- mean(prob_se_cf)
  cf_effect <- (se_cf - se_base) / se_base * 100

  list(K = K, loglik = ll, bic_standard = bic_std, bic_panel = bic_panel,
       n_params = np, pi = pi_k,
       gamma = sapply(params, function(p) p$gamma),
       cf_effect = cf_effect, se_rate_baseline = se_base,
       converged = iter < max_iter)
}


# ==============================================================================
# 5. Three-Pronged Selection
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

run_simulation <- function(dgp, N_values, R_values, K_max = 5) {
  cat(sprintf("\n=== Design: %s ===\n", dgp$name))
  all_results <- list(); idx <- 0

  for (n_i in seq_along(N_values)) {
    N <- N_values[n_i]; R <- R_values[n_i]
    cat(sprintf("\n  N = %d, R = %d:\n", N, R))
    t0 <- proc.time()

    for (r in 1:R) {
      if (r %% 10 == 0 || r == 1) {
        elapsed <- (proc.time() - t0)[3]
        eta <- if (r > 1) elapsed / (r - 1) * (R - r + 1) else NA
        cat(sprintf("    Rep %d/%d (%.0fs elapsed, ETA %.0fs)\n", r, R, elapsed,
                    ifelse(is.na(eta), 0, eta)))
      }

      data <- simulate_data(N, dgp)
      results_k <- list()
      for (K in 1:K_max) {
        tryCatch({
          results_k[[K]] <- estimate_mixture(data, K)
        }, error = function(e) {
          results_k[[K]] <<- list(K = K, loglik = NA, bic_standard = Inf,
            bic_panel = Inf, cf_effect = NA, gamma = rep(NA, K),
            pi = rep(1/K, K), converged = FALSE)
        })
      }

      valid <- results_k[!sapply(results_k, function(r)
        is.null(r) || is.na(r$loglik))]
      sel <- if (length(valid) >= 2) three_prong_select(valid) else
        list(k_bic=NA, k_stable=NA, k_osce=NA, k_consensus=NA)

      for (K in 1:K_max) {
        rk <- results_k[[K]]
        if (!is.null(rk) && !is.na(rk$loglik)) {
          idx <- idx + 1
          all_results[[idx]] <- data.frame(
            design = dgp$name, N = N, rep = r, K_est = K,
            loglik = rk$loglik, bic_standard = rk$bic_standard,
            bic_panel = rk$bic_panel, cf_effect = rk$cf_effect,
            gamma_weighted = sum(rk$pi * rk$gamma),
            k_selected_bic = sel$k_bic, k_selected_stable = sel$k_stable,
            k_selected_osce = sel$k_osce, k_selected_consensus = sel$k_consensus,
            stringsAsFactors = FALSE)
        }
      }
    }
    cat(sprintf("  N=%d complete: %.0f sec\n", N, (proc.time()-t0)[3]))
  }
  bind_rows(all_results)
}


# ==============================================================================
# 7. Run Simulations (v2: doubled replications)
# ==============================================================================

N_values <- c(2000, 10000, 50000)
R_values <- c(200, 100, 20)   # v2: doubled reps (was 100, 50, 10)

cat("\n============================================================\n")
cat("RUNNING SIMULATIONS (v2: doubled replications)\n")
cat(sprintf("  N values: %s\n", paste(N_values, collapse=", ")))
cat(sprintf("  R values: %s\n", paste(R_values, collapse=", ")))
cat("============================================================\n")

t_start <- proc.time()
results_cancel <- run_simulation(dgp_cancellation, N_values, R_values)
results_same <- run_simulation(dgp_samesign, N_values, R_values)
results_all <- bind_rows(results_cancel, results_same)
t_total <- (proc.time() - t_start)[3]
cat(sprintf("\nTotal simulation time: %.0f sec (%.1f hours)\n", t_total, t_total/3600))


# ==============================================================================
# 8. Summary Statistics (v2: with MC standard errors)
# ==============================================================================

cat("\n============================================================\n")
cat("SUMMARY STATISTICS (v2)\n")
cat("============================================================\n")

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

# v2: cf_stats now includes MC standard error (mc_se = sd / sqrt(n))
cf_stats <- results_all %>%
  group_by(design, N, K_est) %>%
  summarise(cf_mean = mean(cf_effect, na.rm=TRUE),
            cf_sd = sd(cf_effect, na.rm=TRUE),
            mc_se = sd(cf_effect, na.rm=TRUE) / sqrt(n()),
            n_reps = n(),
            .groups="drop")

cat("\nCounterfactual Statistics by K (with MC SE):\n")
print(as.data.frame(cf_stats))

cf_range <- results_all %>%
  group_by(design, N, rep) %>%
  summarise(cf_range = max(cf_effect, na.rm=TRUE) - min(cf_effect, na.rm=TRUE),
            .groups="drop") %>%
  group_by(design, N) %>%
  summarise(mean_range = mean(cf_range, na.rm=TRUE),
            median_range = median(cf_range, na.rm=TRUE),
            mc_se_range = sd(cf_range, na.rm=TRUE) / sqrt(n()),
            .groups="drop")

cat("\nCounterfactual Range Across K (with MC SE):\n")
print(as.data.frame(cf_range))


# ==============================================================================
# 9. Save Results (v2 filenames)
# ==============================================================================

write.csv(results_all, file.path(DATA_DIR, "output/monte_carlo_results_v2.csv"), row.names=FALSE)
write.csv(k_freq, file.path(DATA_DIR, "output/monte_carlo_k_selection_v2.csv"), row.names=FALSE)
write.csv(cf_stats, file.path(DATA_DIR, "output/monte_carlo_cf_stats_v2.csv"), row.names=FALSE)
write.csv(cf_range, file.path(DATA_DIR, "output/monte_carlo_cf_range_v2.csv"), row.names=FALSE)
cat("\nCSV results saved to Data/output/ (v2 files)\n")


# ==============================================================================
# 10. Generate Figures (v2 filenames)
# ==============================================================================

cat("\nGENERATING FIGURES (v2)\n")

# MC1: K selection frequency
tryCatch({
  fig_mc1 <- k_freq %>%
    filter(design == "Cancellation") %>%
    ggplot(aes(x = factor(k_selected_consensus), y = freq, fill = factor(N))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("#91bfdb", "#4575b4", "#2166ac"),
                      labels = c("N=2,000", "N=10,000", "N=50,000"), name = "Sample Size") +
    labs(title = "K Selection Frequency: Cancellation DGP (True K=4)",
         x = "Selected K", y = "Frequency") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTPUT_DIR, "fig_mc1_k_selection_freq_v2.pdf"), fig_mc1, width=8, height=5)
  ggsave(file.path(OUTPUT_DIR, "fig_mc1_k_selection_freq_v2.png"), fig_mc1, width=8, height=5, dpi=300)
  cat("  Saved: fig_mc1_v2\n")
}, error = function(e) cat(sprintf("  fig_mc1 error: %s\n", e$message)))

# MC2: CF by design
tryCatch({
  cf_by_k <- results_all %>%
    group_by(design, N, K_est) %>%
    summarise(cf_mean = mean(cf_effect, na.rm=TRUE),
              cf_q25 = quantile(cf_effect, 0.25, na.rm=TRUE),
              cf_q75 = quantile(cf_effect, 0.75, na.rm=TRUE), .groups="drop")
  fig_mc2 <- cf_by_k %>% filter(N == 10000) %>%
    ggplot(aes(x = K_est, y = cf_mean, color = design)) +
    geom_line(size = 1.2) + geom_point(size = 3) +
    geom_ribbon(aes(ymin=cf_q25, ymax=cf_q75, fill=design), alpha=0.2, color=NA) +
    scale_color_manual(values = c("#d73027", "#4575b4"), name = "DGP") +
    scale_fill_manual(values = c("#d73027", "#4575b4"), name = "DGP") +
    labs(title = "Counterfactual Sensitivity by DGP (N=10,000)",
         x = "K", y = "Mean CF Effect (%)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTPUT_DIR, "fig_mc2_cf_by_design_v2.pdf"), fig_mc2, width=8, height=5)
  ggsave(file.path(OUTPUT_DIR, "fig_mc2_cf_by_design_v2.png"), fig_mc2, width=8, height=5, dpi=300)
  cat("  Saved: fig_mc2_v2\n")
}, error = function(e) cat(sprintf("  fig_mc2 error: %s\n", e$message)))

# MC3: CF range distribution
tryCatch({
  cf_range_dist <- results_all %>%
    group_by(design, N, rep) %>%
    summarise(cf_range = max(cf_effect, na.rm=TRUE) - min(cf_effect, na.rm=TRUE), .groups="drop")
  fig_mc3 <- cf_range_dist %>% filter(N == 10000) %>%
    ggplot(aes(x = cf_range, fill = design)) +
    geom_histogram(bins = 25, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("#d73027", "#4575b4"), name = "DGP") +
    labs(title = "Distribution of CF Range Across K",
         x = "Range (pp)", y = "Count") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTPUT_DIR, "fig_mc3_cf_range_dist_v2.pdf"), fig_mc3, width=8, height=5)
  ggsave(file.path(OUTPUT_DIR, "fig_mc3_cf_range_dist_v2.png"), fig_mc3, width=8, height=5, dpi=300)
  cat("  Saved: fig_mc3_v2\n")
}, error = function(e) cat(sprintf("  fig_mc3 error: %s\n", e$message)))

# MC4: BIC accuracy
tryCatch({
  bic_acc <- results_all %>%
    filter(K_est == k_selected_bic) %>%
    distinct(design, N, rep, .keep_all=TRUE) %>%
    mutate(correct = k_selected_bic == 4) %>%
    group_by(design, N) %>%
    summarise(accuracy = mean(correct, na.rm=TRUE), .groups="drop")
  fig_mc4 <- bic_acc %>%
    ggplot(aes(x = factor(N, labels=c("2,000","10,000","50,000")),
               y = accuracy, fill = design)) +
    geom_bar(stat="identity", position="dodge", width=0.6) +
    scale_fill_manual(values = c("#d73027", "#4575b4"), name = "DGP") +
    geom_hline(yintercept=1, linetype="dashed", color="gray50") +
    scale_y_continuous(limits=c(0,1), labels=function(x) paste0(x*100,"%")) +
    labs(title = "BIC Selection Accuracy (True K=4)", x = "Sample Size",
         y = "Proportion Selecting K=4") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTPUT_DIR, "fig_mc4_bic_accuracy_v2.pdf"), fig_mc4, width=8, height=5)
  ggsave(file.path(OUTPUT_DIR, "fig_mc4_bic_accuracy_v2.png"), fig_mc4, width=8, height=5, dpi=300)
  cat("  Saved: fig_mc4_v2\n")
}, error = function(e) cat(sprintf("  fig_mc4 error: %s\n", e$message)))


cat("\n============================================================\n")
cat("MONTE CARLO SIMULATION v2 COMPLETE\n")
cat(sprintf("Total time: %.1f hours\n", t_total / 3600))
cat("============================================================\n")
