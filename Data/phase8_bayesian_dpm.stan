/*
Phase 8: Bayesian Nonparametric Mixture Model (Dirichlet Process)

Following:
- Frühwirth-Schnatter (2006)
- Malsiner-Walli, Frühwirth-Schnatter, Grün (2016, Statistics and Computing)

Key idea: Instead of selecting K via BIC, place a Dirichlet prior on mixing weights.
Redundant components get shrunk to near-zero posterior weight.

Output: Posterior distribution over K and counterfactual effects.
Report: "posterior probability of K >= 4 is 0.XX" rather than "BIC selects K=4"
*/

data {
  int<lower=1> N;                    // Number of observations
  int<lower=1> J;                    // Number of alternatives (9)
  int<lower=1> P;                    // Number of covariates
  int<lower=1> K_max;                // Maximum number of types (e.g., 10)

  array[N] int<lower=1,upper=J> y;   // Observed choices
  matrix[N, P] X;                    // Covariates
  vector[N] branch;                  // Branch banking indicator
  vector[N] mobile;                  // Mobile banking indicator
  vector[N] branch_density;          // CBSA branch density
  vector[N] weights;                 // Survey weights
}

parameters {
  // Mixing weights with symmetric Dirichlet prior
  simplex[K_max] pi;

  // Type-specific branch effects
  vector[K_max] gamma_branch;

  // Common parameters across types
  matrix[P, J-1] beta;               // Covariate effects (base = alternative J)
  real gamma_mobile;                 // Mobile effect (common)
  real delta_density;                // Density effect (common)

  // Hyperparameters
  real<lower=0> sigma_gamma;         // SD of branch effects across types
}

transformed parameters {
  // Compute choice probabilities for each observation and type
  array[N, K_max] vector[J] log_prob;

  for (i in 1:N) {
    for (k in 1:K_max) {
      vector[J] V;

      // Utility for each alternative
      // Base alternative (J) has V = 0
      V[J] = 0;

      for (j in 1:(J-1)) {
        V[j] = X[i] * beta[, j];

        // Add banking mode effects for SE alternatives (j = 2, 5, 8)
        // Assuming: 1-3 = Unbanked, 4-6 = Mobile, 7-9 = Branch
        // SE alternatives: 2, 5, 8
        if (j == 2 || j == 5 || j == 8) {
          // Branch effect (type-specific)
          if (j == 8) {  // Branch x SE
            V[j] += gamma_branch[k] * branch_density[i];
          }
          // Mobile effect
          if (j == 5) {  // Mobile x SE
            V[j] += gamma_mobile;
          }
        }
      }

      // Log softmax for numerical stability
      log_prob[i, k] = log_softmax(V);
    }
  }
}

model {
  // Priors

  // Symmetric Dirichlet prior on mixing weights (sparse)
  // Small alpha encourages few active components
  pi ~ dirichlet(rep_vector(0.1, K_max));

  // Hierarchical prior on type-specific branch effects
  sigma_gamma ~ cauchy(0, 1);
  gamma_branch ~ normal(0, sigma_gamma);

  // Weakly informative priors on common parameters
  to_vector(beta) ~ normal(0, 2);
  gamma_mobile ~ normal(0, 1);
  delta_density ~ normal(0, 1);

  // Likelihood (mixture)
  for (i in 1:N) {
    vector[K_max] log_lik_k;

    for (k in 1:K_max) {
      log_lik_k[k] = log(pi[k]) + log_prob[i, k, y[i]];
    }

    // Weighted log-likelihood
    target += weights[i] * log_sum_exp(log_lik_k);
  }
}

generated quantities {
  // Posterior predictive checks
  array[N] int y_rep;

  // Effective number of components
  real K_effective;

  // Type assignments (posterior mode)
  array[N] int<lower=1,upper=K_max> type_assignment;

  // Counterfactual: SE rate under 50% branch closure
  real se_rate_baseline;
  real se_rate_counterfactual;
  real counterfactual_effect;

  // Compute effective K (number of types with pi > 0.01)
  {
    int count = 0;
    for (k in 1:K_max) {
      if (pi[k] > 0.01) count += 1;
    }
    K_effective = count;
  }

  // Type assignments and predictions
  for (i in 1:N) {
    vector[K_max] log_lik_k;

    for (k in 1:K_max) {
      log_lik_k[k] = log(pi[k]) + log_prob[i, k, y[i]];
    }

    // Assign to most likely type
    type_assignment[i] = 1;
    for (k in 2:K_max) {
      if (log_lik_k[k] > log_lik_k[type_assignment[i]]) {
        type_assignment[i] = k;
      }
    }

    // Posterior predictive
    {
      int k_draw = categorical_rng(pi);
      y_rep[i] = categorical_rng(softmax(log_prob[i, k_draw]));
    }
  }

  // Counterfactual calculation
  {
    real se_base = 0;
    real se_cf = 0;
    real total_weight = 0;

    for (i in 1:N) {
      // Baseline SE probability (weighted average over types)
      real p_se_base = 0;
      real p_se_cf = 0;

      for (k in 1:K_max) {
        // SE alternatives are 2, 5, 8
        real p_se_k = exp(log_prob[i, k, 2]) + exp(log_prob[i, k, 5]) + exp(log_prob[i, k, 8]);
        p_se_base += pi[k] * p_se_k;

        // Counterfactual: reduce branch density by 50% for branch users
        // Simplified: reduce gamma_branch effect
        real p_se_cf_k = p_se_k;
        if (branch[i] > 0.5) {
          // Approximate counterfactual effect
          p_se_cf_k = p_se_k * exp(-0.5 * gamma_branch[k] * branch_density[i]);
        }
        p_se_cf += pi[k] * p_se_cf_k;
      }

      se_base += weights[i] * p_se_base;
      se_cf += weights[i] * p_se_cf;
      total_weight += weights[i];
    }

    se_rate_baseline = se_base / total_weight;
    se_rate_counterfactual = se_cf / total_weight;
    counterfactual_effect = (se_rate_counterfactual - se_rate_baseline) / se_rate_baseline * 100;
  }
}
