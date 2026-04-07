#' ==============================================================================
#' Monte Carlo Simulation: Severity Grid (Comment 3b)
#'
#' Varies the severity of cancellation across a grid of DGPs.
#' For each weighted average gamma_bar in {0, 0.02, 0.05, 0.10, 0.15},
#' constructs gamma vectors that produce that weighted average while
#' maintaining cancellation (opposite signs for at least some types).
#'
#' Output:
#'   - Data/output/monte_carlo_severity_grid.csv     (replication-level)
#'   - Data/output/monte_carlo_severity_summary.csv  (summary)
#'   - Data/output/figures_methods/fig_mc5_severity_grid.pdf/.png
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
cat("MONTE CARLO SEVERITY GRID (Comment 3b)\n")
cat("Varying cancellation severity across DGPs\n")
cat("============================================================\n\n")


# ==============================================================================
# 1. Constants (same as monte_carlo_simulation.R)
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


# ==============================================================================
# 2. Core Functions (copied from monte_carlo_simulation.R)
# ==============================================================================

simulate_data <- function(N, dgp, n_cbsa = 50) {
  K <- dgp$K_true
  type <- sample(1:K, N, replace = TRUE, prob = dgp$pi)
  cbsa <- sample(1:n_cbsa, N, replace = TRUE)
  x1 <- rnorm(N, 0, dgp$sigma_x)
  x2 <- rbinom(N, 1, 0.4)

  # Build N x J utility matrix using vectorized operations
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

softmax_mat <- function(U) {
  # Numerically stable row-wise softmax using vectorized pmax
  rm <- do.call(pmax, as.data.frame(U))  # row max, much faster than apply
  e <- exp(U - rm)
  rs <- rowSums(e)
  rs[rs == 0] <- 1e-300
  e / rs
}

compute_utils_k <- function(x1, x2, N, pk) {
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
# 3. Define Severity Grid DGPs
# ==============================================================================

# Common parameters across all DGPs
K_true <- 4
pi_vec <- c(0.13, 0.20, 0.35, 0.32)
N_sim  <- 10000
R_sim  <- 50
K_max  <- 5

# Gamma vectors for each gamma_bar level
# Each produces the target weighted average: sum(pi * gamma) = gamma_bar
# while maintaining cancellation (opposite signs for at least some types)
#
# Verification:
#   gamma_bar = 0.00: sum(c(0.13,0.20,0.35,0.32) * c(-0.10,-0.05,0.05,0.10))
#                   = -0.013 - 0.010 + 0.0175 + 0.032 = 0.0265 (approx)
# We use exact constructions below that hit the target precisely.

gamma_grid <- list(
  list(
    gamma_bar = 0.00,
    gamma = c(-0.10, -0.05, 0.05, 0.10),
    label = "gamma_bar = 0 (full cancellation)"
  ),
  list(
    gamma_bar = 0.02,
    gamma = c(-0.08, -0.03, 0.05, 0.10),
    label = "gamma_bar = 0.02 (strong cancellation)"
  ),
  list(
    gamma_bar = 0.05,
    gamma = c(-0.05, 0.00, 0.08, 0.12),
    label = "gamma_bar = 0.05 (moderate cancellation)"
  ),
  list(
    gamma_bar = 0.10,
    gamma = c(0.00, 0.05, 0.12, 0.18),
    label = "gamma_bar = 0.10 (mild cancellation)"
  ),
  list(
    gamma_bar = 0.15,
    gamma = c(0.05, 0.10, 0.18, 0.22),
    label = "gamma_bar = 0.15 (no cancellation)"
  )
)

# Print the DGP grid with verification
cat("Severity Grid DGPs:\n")
cat(sprintf("  Common: K_true=%d, pi=(%s), N=%d, R=%d\n",
            K_true, paste(pi_vec, collapse=", "), N_sim, R_sim))
cat("  ---\n")
for (g in gamma_grid) {
  actual_bar <- sum(pi_vec * g$gamma)
  cat(sprintf("  gamma_bar=%.2f: gamma=(%s) -> actual weighted avg=%.4f  [%s]\n",
              g$gamma_bar,
              paste(sprintf("%.2f", g$gamma), collapse=", "),
              actual_bar,
              g$label))
}
cat("\n")


# ==============================================================================
# 4. Main Simulation Loop
# ==============================================================================

cat("============================================================\n")
cat("RUNNING SEVERITY GRID SIMULATIONS\n")
cat(sprintf("  %d DGPs x %d replications x K=1..%d = %d total estimations\n",
            length(gamma_grid), R_sim, K_max,
            length(gamma_grid) * R_sim * K_max))
cat("============================================================\n\n")

t_start <- proc.time()
all_results <- list()
idx <- 0

for (g_idx in seq_along(gamma_grid)) {
  g <- gamma_grid[[g_idx]]
  gamma_bar <- g$gamma_bar

  # Construct DGP
  dgp <- list(
    name = sprintf("gamma_bar_%.2f", gamma_bar),
    K_true = K_true,
    pi = pi_vec,
    gamma = g$gamma,
    beta_intercept = c(-1.5, -1.2, -1.0, -0.8),
    sigma_x = 1.0
  )

  cat(sprintf("\n--- DGP %d/%d: %s ---\n", g_idx, length(gamma_grid), g$label))
  t0 <- proc.time()

  for (r in 1:R_sim) {
    if (r %% 10 == 0 || r == 1) {
      elapsed <- (proc.time() - t0)[3]
      eta <- if (r > 1) elapsed / (r - 1) * (R_sim - r + 1) else NA
      cat(sprintf("  Rep %d/%d (%.0fs elapsed, ETA %.0fs)\n", r, R_sim, elapsed,
                  ifelse(is.na(eta), 0, eta)))
    }

    data <- simulate_data(N_sim, dgp)
    results_k <- list()

    for (K in 1:K_max) {
      tryCatch({
        results_k[[K]] <- estimate_mixture(data, K)
      }, error = function(e) {
        results_k[[K]] <<- list(
          K = K, loglik = NA, bic_standard = Inf,
          bic_panel = Inf, cf_effect = NA, gamma = rep(NA, K),
          pi = rep(1/K, K), converged = FALSE
        )
      })
    }

    # Three-prong selection
    valid <- results_k[!sapply(results_k, function(r)
      is.null(r) || is.na(r$loglik))]
    sel <- if (length(valid) >= 2) three_prong_select(valid) else
      list(k_bic = NA, k_stable = NA, k_osce = NA, k_consensus = NA)

    # Store results for each K
    for (K in 1:K_max) {
      rk <- results_k[[K]]
      if (!is.null(rk) && !is.na(rk$loglik)) {
        idx <- idx + 1
        all_results[[idx]] <- data.frame(
          gamma_bar = gamma_bar,
          rep = r,
          K_est = K,
          cf_effect = rk$cf_effect,
          k_selected_consensus = sel$k_consensus,
          gamma_weighted = sum(rk$pi * rk$gamma),
          bic_panel = rk$bic_panel,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  elapsed_dgp <- (proc.time() - t0)[3]
  cat(sprintf("  DGP gamma_bar=%.2f complete: %.0f sec (%.1f min)\n",
              gamma_bar, elapsed_dgp, elapsed_dgp / 60))
}

t_total <- (proc.time() - t_start)[3]
cat(sprintf("\nTotal simulation time: %.0f sec (%.1f hours)\n", t_total, t_total / 3600))


# ==============================================================================
# 5. Combine Results and Save Replication-Level CSV
# ==============================================================================

results_df <- bind_rows(all_results)

cat("\n============================================================\n")
cat("RESULTS SUMMARY\n")
cat("============================================================\n")
cat(sprintf("Total rows: %d\n", nrow(results_df)))
cat(sprintf("DGPs: %s\n", paste(unique(results_df$gamma_bar), collapse = ", ")))

# Save replication-level results
write.csv(results_df,
          file.path(DATA_DIR, "output/monte_carlo_severity_grid.csv"),
          row.names = FALSE)
cat("\nSaved: output/monte_carlo_severity_grid.csv\n")


# ==============================================================================
# 6. Compute Summary Statistics
# ==============================================================================

# K recovery rate: proportion of replications where consensus selects true K=4
k_recovery <- results_df %>%
  filter(K_est == k_selected_consensus) %>%
  distinct(gamma_bar, rep, .keep_all = TRUE) %>%
  group_by(gamma_bar) %>%
  summarise(
    k_recovery_rate = mean(k_selected_consensus == K_true, na.rm = TRUE),
    .groups = "drop"
  )

# CF range across K within each replication
cf_range_by_rep <- results_df %>%
  group_by(gamma_bar, rep) %>%
  summarise(
    cf_range = max(cf_effect, na.rm = TRUE) - min(cf_effect, na.rm = TRUE),
    .groups = "drop"
  )

mean_cf_range <- cf_range_by_rep %>%
  group_by(gamma_bar) %>%
  summarise(
    mean_cf_range = mean(cf_range, na.rm = TRUE),
    .groups = "drop"
  )

# Mean CF at true K=4
cf_at_true_k <- results_df %>%
  filter(K_est == K_true) %>%
  group_by(gamma_bar) %>%
  summarise(
    mean_cf_at_true_k = mean(cf_effect, na.rm = TRUE),
    mc_se_cf = sd(cf_effect, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Merge into summary
summary_df <- k_recovery %>%
  left_join(mean_cf_range, by = "gamma_bar") %>%
  left_join(cf_at_true_k, by = "gamma_bar")

cat("\nSeverity Grid Summary:\n")
print(as.data.frame(summary_df))

# Save summary
write.csv(summary_df,
          file.path(DATA_DIR, "output/monte_carlo_severity_summary.csv"),
          row.names = FALSE)
cat("\nSaved: output/monte_carlo_severity_summary.csv\n")


# ==============================================================================
# 7. Additional Diagnostics
# ==============================================================================

# K selection frequency by gamma_bar
k_freq <- results_df %>%
  filter(K_est == k_selected_consensus) %>%
  distinct(gamma_bar, rep, .keep_all = TRUE) %>%
  group_by(gamma_bar, k_selected_consensus) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(gamma_bar) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

cat("\nK Selection Frequency by Gamma Bar:\n")
print(as.data.frame(k_freq))

# Mean CF by K and gamma_bar
cf_by_k <- results_df %>%
  group_by(gamma_bar, K_est) %>%
  summarise(
    cf_mean = mean(cf_effect, na.rm = TRUE),
    cf_sd = sd(cf_effect, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nCounterfactual by K and Gamma Bar:\n")
print(as.data.frame(cf_by_k))


# ==============================================================================
# 8. Generate Severity Heatmap Figure (fig_mc5)
# ==============================================================================

cat("\nGENERATING SEVERITY HEATMAP FIGURE\n")

tryCatch({
  # Prepare heatmap data: two metrics across gamma_bar
  # Metric 1: K recovery rate
  # Metric 2: Mean CF range across K

  heatmap_data <- summary_df %>%
    select(gamma_bar, k_recovery_rate, mean_cf_range) %>%
    pivot_longer(
      cols = c(k_recovery_rate, mean_cf_range),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric_label = case_when(
        metric == "k_recovery_rate" ~ "K Recovery Rate",
        metric == "mean_cf_range" ~ "CF Range Across K (pp)"
      ),
      # For display: scale K recovery to percentage
      display_value = ifelse(metric == "k_recovery_rate", value * 100, value)
    )

  # Create heatmap
  fig_mc5 <- ggplot(heatmap_data,
                     aes(x = factor(gamma_bar),
                         y = metric_label,
                         fill = display_value)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%.1f", display_value)),
              color = "white", size = 5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#f7f7f7",
      high = "#d73027",
      midpoint = median(heatmap_data$display_value),
      name = "Value"
    ) +
    # Override text color for light cells
    geom_text(
      data = heatmap_data %>%
        filter(abs(display_value - median(heatmap_data$display_value)) <
                 0.3 * diff(range(heatmap_data$display_value))),
      aes(label = sprintf("%.1f", display_value)),
      color = "gray20", size = 5, fontface = "bold"
    ) +
    labs(
      title = expression(paste(
        "Severity of Cancellation: K Recovery & CF Sensitivity by ", bar(gamma))),
      subtitle = paste0("K_true=4, N=", format(N_sim, big.mark=","),
                        ", R=", R_sim, " replications per DGP"),
      x = expression(paste("Weighted Average Treatment Effect  ", bar(gamma))),
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40", size = 11),
      axis.text = element_text(size = 12),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  ggsave(file.path(OUTPUT_DIR, "fig_mc5_severity_grid.pdf"),
         fig_mc5, width = 9, height = 4)
  ggsave(file.path(OUTPUT_DIR, "fig_mc5_severity_grid.png"),
         fig_mc5, width = 9, height = 4, dpi = 300)
  cat("  Saved: fig_mc5_severity_grid.pdf and .png\n")

}, error = function(e) {
  cat(sprintf("  fig_mc5 error: %s\n", e$message))
})


# ==============================================================================
# 9. Supplementary Figure: CF by K across Severity Levels
# ==============================================================================

tryCatch({
  fig_mc5b <- cf_by_k %>%
    ggplot(aes(x = K_est, y = cf_mean,
               color = factor(gamma_bar),
               group = factor(gamma_bar))) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_ribbon(
      aes(ymin = cf_mean - cf_sd, ymax = cf_mean + cf_sd,
          fill = factor(gamma_bar)),
      alpha = 0.1, color = NA
    ) +
    scale_color_manual(
      values = c("#2166ac", "#4393c3", "#92c5de", "#f4a582", "#d73027"),
      name = expression(bar(gamma)),
      labels = c("0.00", "0.02", "0.05", "0.10", "0.15")
    ) +
    scale_fill_manual(
      values = c("#2166ac", "#4393c3", "#92c5de", "#f4a582", "#d73027"),
      name = expression(bar(gamma)),
      labels = c("0.00", "0.02", "0.05", "0.10", "0.15")
    ) +
    geom_vline(xintercept = K_true, linetype = "dashed", color = "gray50",
               linewidth = 0.5) +
    annotate("text", x = K_true + 0.15, y = max(cf_by_k$cf_mean, na.rm = TRUE),
             label = "True K", hjust = 0, color = "gray40", size = 3.5) +
    labs(
      title = expression(paste(
        "Counterfactual Sensitivity by K Across Severity Levels")),
      subtitle = paste0("N=", format(N_sim, big.mark=","),
                        ", R=", R_sim, " replications; ",
                        "ribbons show +/- 1 SD"),
      x = "Estimated K",
      y = "Mean CF Effect (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40", size = 11),
      legend.position = "right"
    )

  ggsave(file.path(OUTPUT_DIR, "fig_mc5b_cf_by_severity.pdf"),
         fig_mc5b, width = 9, height = 5.5)
  ggsave(file.path(OUTPUT_DIR, "fig_mc5b_cf_by_severity.png"),
         fig_mc5b, width = 9, height = 5.5, dpi = 300)
  cat("  Saved: fig_mc5b_cf_by_severity.pdf and .png\n")

}, error = function(e) {
  cat(sprintf("  fig_mc5b error: %s\n", e$message))
})


cat("\n============================================================\n")
cat("MONTE CARLO SEVERITY GRID COMPLETE\n")
cat(sprintf("Total time: %.1f hours\n", t_total / 3600))
cat("============================================================\n")
