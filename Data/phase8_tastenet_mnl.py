"""
Phase 8: TasteNet-MNL - Neural Network-Enhanced Discrete Choice Model

Following Han et al. (2022, Transportation Research Part B)

The standard MNL has linear utility:
    V_ij = X_i * beta_j + C_j * gamma_C + Z_j * delta_j

TasteNet replaces fixed gamma_C with a neural network that learns gamma_C(X_i):
    V_ij = X_i * beta_j + C_j * NN_theta(X_i) + Z_j * delta_j

Key innovation:
- Keep MNL choice probability structure (counterfactuals well-defined)
- NN discovers which individuals are most sensitive to branch access
- Continuous analog to finite mixture

If NN learns a function that clusters into ~4 groups, confirms finite mixture.
If NN finds smooth continuous function, suggests mixed logit is better.
"""

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# ============================================================================
# 1. TasteNet Module
# ============================================================================

class TasteNet(nn.Module):
    """
    Neural network that outputs individual-specific taste parameters.

    Input: Individual demographics X_i
    Output: gamma_C(X_i) - individual's sensitivity to branch access
    """

    def __init__(self, input_dim, hidden_dims=[32, 16], output_dim=1):
        super(TasteNet, self).__init__()

        layers = []
        prev_dim = input_dim

        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(prev_dim, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.BatchNorm1d(hidden_dim))
            layers.append(nn.Dropout(0.2))
            prev_dim = hidden_dim

        layers.append(nn.Linear(prev_dim, output_dim))

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)


# ============================================================================
# 2. TasteNet-MNL Model
# ============================================================================

class TasteNetMNL(nn.Module):
    """
    Multinomial Logit with TasteNet for heterogeneous taste parameters.

    P(y_i = j) = exp(V_ij) / sum_k exp(V_ik)

    where V_ij = X_i @ beta_j + C_j * TasteNet(X_i) + Z_j @ delta_j
    """

    def __init__(self, n_alternatives, x_dim, z_dim, tastenet_hidden=[32, 16]):
        super(TasteNetMNL, self).__init__()

        self.n_alternatives = n_alternatives

        # Alternative-specific parameters (beta)
        self.beta = nn.Parameter(torch.randn(x_dim, n_alternatives - 1) * 0.01)

        # Infrastructure parameters (delta)
        self.delta = nn.Parameter(torch.randn(z_dim, n_alternatives - 1) * 0.01)

        # TasteNet for heterogeneous credit access sensitivity
        self.tastenet = TasteNet(x_dim, hidden_dims=tastenet_hidden, output_dim=1)

        # Credit access indicator for each alternative (fixed)
        # Assumes: alternatives 0-2 = Unbanked x {Wage, SE, NotWork}
        #          alternatives 3-5 = Mobile x {Wage, SE, NotWork}
        #          alternatives 6-8 = Branch x {Wage, SE, NotWork}
        # SE alternatives: 1, 4, 7
        # Branch alternatives: 6, 7, 8

    def forward(self, X, Z, credit_access, se_indicator):
        """
        X: Individual demographics [batch, x_dim]
        Z: CBSA infrastructure [batch, z_dim]
        credit_access: Credit access level by alternative [batch, n_alt]
        se_indicator: 1 if alternative is self-employment [batch, n_alt]

        Returns: Log-probabilities [batch, n_alternatives]
        """
        batch_size = X.shape[0]

        # Base utility from demographics
        V_base = X @ self.beta  # [batch, n_alt-1]
        V_base = torch.cat([torch.zeros(batch_size, 1), V_base], dim=1)  # Add base alt

        # Infrastructure utility
        V_infra = Z @ self.delta  # [batch, n_alt-1]
        V_infra = torch.cat([torch.zeros(batch_size, 1), V_infra], dim=1)

        # TasteNet: individual-specific branch sensitivity
        gamma_i = self.tastenet(X)  # [batch, 1]

        # Credit access effect: only for SE alternatives
        V_credit = gamma_i * credit_access * se_indicator  # [batch, n_alt]

        # Total utility
        V = V_base + V_infra + V_credit

        # Log-softmax for numerical stability
        log_probs = torch.log_softmax(V, dim=1)

        return log_probs, gamma_i

    def predict_proba(self, X, Z, credit_access, se_indicator):
        """Return choice probabilities"""
        with torch.no_grad():
            log_probs, _ = self.forward(X, Z, credit_access, se_indicator)
            return torch.exp(log_probs)


# ============================================================================
# 3. Training Loop
# ============================================================================

def train_tastenet_mnl(model, train_loader, val_loader, epochs=100, lr=0.001):
    """
    Train TasteNet-MNL with early stopping.
    """
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=10)
    criterion = nn.NLLLoss()

    best_val_loss = float('inf')
    patience_counter = 0

    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0
        for X, Z, credit, se_ind, y in train_loader:
            optimizer.zero_grad()
            log_probs, _ = model(X, Z, credit, se_ind)
            loss = criterion(log_probs, y)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()

        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for X, Z, credit, se_ind, y in val_loader:
                log_probs, _ = model(X, Z, credit, se_ind)
                loss = criterion(log_probs, y)
                val_loss += loss.item()

        train_loss /= len(train_loader)
        val_loss /= len(val_loader)

        scheduler.step(val_loss)

        if epoch % 10 == 0:
            print(f"Epoch {epoch}: Train Loss = {train_loss:.4f}, Val Loss = {val_loss:.4f}")

        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            torch.save(model.state_dict(), 'best_tastenet_mnl.pt')
        else:
            patience_counter += 1
            if patience_counter >= 20:
                print(f"Early stopping at epoch {epoch}")
                break

    # Load best model
    model.load_state_dict(torch.load('best_tastenet_mnl.pt'))
    return model


# ============================================================================
# 4. Analyze Learned Taste Heterogeneity
# ============================================================================

def analyze_taste_heterogeneity(model, X, feature_names):
    """
    Analyze the learned gamma_C(X) function from TasteNet.
    """
    model.eval()
    with torch.no_grad():
        gamma_i = model.tastenet(X).numpy().flatten()

    print("\n" + "="*60)
    print("LEARNED TASTE HETEROGENEITY")
    print("="*60)

    # Distribution of learned coefficients
    print(f"\nDistribution of gamma_C(X_i):")
    print(f"  Mean:   {np.mean(gamma_i):.4f}")
    print(f"  Std:    {np.std(gamma_i):.4f}")
    print(f"  Min:    {np.min(gamma_i):.4f}")
    print(f"  Max:    {np.max(gamma_i):.4f}")

    # Percentiles
    percentiles = [5, 10, 25, 50, 75, 90, 95]
    print(f"\nPercentiles:")
    for p in percentiles:
        val = np.percentile(gamma_i, p)
        print(f"  {p}th: {val:.4f}")

    # Check for clustering (comparison to finite mixture)
    from sklearn.cluster import KMeans

    print(f"\nK-Means Clustering (checking for discrete types):")
    for k in [2, 3, 4, 5]:
        kmeans = KMeans(n_clusters=k, random_state=42)
        labels = kmeans.fit_predict(gamma_i.reshape(-1, 1))

        # Compute cluster centers and shares
        centers = kmeans.cluster_centers_.flatten()
        shares = [np.mean(labels == i) for i in range(k)]

        print(f"\n  K={k}:")
        for i in range(k):
            print(f"    Type {i+1}: center={centers[i]:.4f}, share={shares[i]*100:.1f}%")

    return gamma_i


# ============================================================================
# 5. Counterfactual Analysis
# ============================================================================

def counterfactual_analysis(model, X, Z, credit_access, se_indicator,
                            closure_rate=0.5):
    """
    Compute counterfactual SE rate under branch closure.

    For branch users, reduce credit_access by closure_rate.
    """
    model.eval()

    with torch.no_grad():
        # Baseline
        probs_base = model.predict_proba(X, Z, credit_access, se_indicator)
        se_alts = [1, 4, 7]  # SE alternatives
        se_rate_base = probs_base[:, se_alts].sum(dim=1).mean().item()

        # Counterfactual: reduce branch credit access
        credit_cf = credit_access.clone()
        branch_alts = [6, 7, 8]
        credit_cf[:, branch_alts] *= (1 - closure_rate)

        probs_cf = model.predict_proba(X, Z, credit_cf, se_indicator)
        se_rate_cf = probs_cf[:, se_alts].sum(dim=1).mean().item()

        effect = (se_rate_cf - se_rate_base) / se_rate_base * 100

    print(f"\nCounterfactual Analysis ({closure_rate*100:.0f}% branch closure):")
    print(f"  Baseline SE rate: {se_rate_base*100:.2f}%")
    print(f"  Counterfactual SE rate: {se_rate_cf*100:.2f}%")
    print(f"  Effect: {effect:+.1f}%")

    return se_rate_base, se_rate_cf, effect


# ============================================================================
# 6. Main Script
# ============================================================================

if __name__ == "__main__":
    print("="*60)
    print("TasteNet-MNL: Neural Network-Enhanced Discrete Choice")
    print("="*60)

    # Note: This script requires actual data loading
    # Below is a demonstration with synthetic data

    print("\nNote: This script demonstrates the TasteNet-MNL architecture.")
    print("For actual estimation, load the FDIC/CPS data and prepare tensors.")
    print("\nKey outputs:")
    print("1. Learned gamma_C(X_i) function from neural network")
    print("2. Distribution of individual-level branch sensitivities")
    print("3. Comparison to K=4 finite mixture (clustering analysis)")
    print("4. Counterfactual SE rates under branch closure")

    # Synthetic demonstration
    np.random.seed(42)
    n_obs = 10000
    x_dim = 5
    z_dim = 2
    n_alternatives = 9

    X = torch.randn(n_obs, x_dim)
    Z = torch.randn(n_obs, z_dim)

    # Credit access: higher for branch alternatives
    credit_access = torch.zeros(n_obs, n_alternatives)
    credit_access[:, 6:9] = 1.0  # Branch alternatives
    credit_access[:, 3:6] = 0.3  # Mobile alternatives (partial)

    # SE indicator
    se_indicator = torch.zeros(n_obs, n_alternatives)
    se_indicator[:, [1, 4, 7]] = 1.0  # SE alternatives

    # Initialize model
    model = TasteNetMNL(n_alternatives, x_dim, z_dim, tastenet_hidden=[32, 16])

    print(f"\nModel architecture:")
    print(model)

    print("\nTo run full estimation:")
    print("1. Load analysis_dataset_with_se.dta")
    print("2. Create choice variable (1-9)")
    print("3. Prepare X, Z, credit_access, se_indicator tensors")
    print("4. Call train_tastenet_mnl()")
    print("5. Analyze heterogeneity with analyze_taste_heterogeneity()")
    print("6. Compute counterfactuals with counterfactual_analysis()")
