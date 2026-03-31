# When K Matters: Model Selection in Finite Mixtures

## Paper
- **Title**: When K Matters: Model Selection in Finite Mixtures and the Magnitude of Policy Counterfactuals
- **Author**: Alina Malkova (Florida Institute of Technology)
- **Target**: JBES / JoE
- **Status**: Draft (March 2026), seeking perspective/comments. Numbers audited Mar 30.
- **Companion paper**: Malkova (2026) — reduced-form mobile banking & self-employment paper (in `Paper/` subfolder)

## Core Argument
The number of mixture components K — typically a footnote — can be the **dominant source of uncertainty** in counterfactual policy analysis. In the application (bank branch access → self-employment), counterfactuals range from −0.6% (K=1) to −11% (K=4), an **18× difference** dwarfing all other specification uncertainty by an order of magnitude.

## Three-Pronged K Selection Framework
1. **Panel BIC** (Hao & Kasahara 2025) — statistical fit → selects K=4
2. **Counterfactual stability** (Bonhomme et al. 2022) — policy relevance → selects K=4/5
3. **OSCE / Penalized MLE** (Budanova 2025) — type significance → identifies 4 active types

## Key Results
- Reduced-form: null effect (0.34pp, p=0.36)
- Structural K=4: −11.0% (95% CI: [−15.6%, −6.4%])
- Bayesian DPM: −8.5% (95% CI: [−14.2%, −3.1%]), P(effect<0)=0.99
- Bounds under model uncertainty: [−11%, −0.6%]
- K selection generates 19× range; next largest source is 2.4× (sampling uncertainty)

## Data
- **Primary**: FDIC National Survey of Unbanked/Underbanked Households (CPS supplement, 2013–2023)
- **Geographic**: FDIC Summary of Deposits (branch-level), ACS CBSA controls
- **Sample**: N=94,886 working-age adults (18–64) in labor force
- **Joint choice**: 9 alternatives (3 banking modes × 3 employment statuses)
- **Analysis dataset**: `Data/analysis_dataset_with_se.dta`

## Project Structure
```
Mobile banking USA/
├── Methods Paper Seeking Perspective.tex/pdf   # Main methods paper
├── Seeking Perspective Methods.pdf             # Earlier draft (Feb 2026)
├── Paper/mobile_banking_self_employment.tex    # Companion reduced-form paper
├── CLAUDE.md                                   # This file
└── Data/
    ├── analysis_dataset_with_se.dta            # Main analysis dataset
    ├── construct_analysis_dataset.do/.R         # Data construction
    ├── compute_ccps.do                         # CCP computation
    │
    ├── # REDUCED FORM (Phase 0)
    ├── phase0_reduced_form.do                  # OLS benchmark
    ├── phase0_continuation.do                  # Continuation analysis
    │
    ├── # STRUCTURAL ESTIMATION (Phases 1-3)
    ├── phase1_structural_mlogit.do             # Base MNL estimation
    ├── phase2_type_selection.do                # K=1–4 BIC comparison
    ├── phase2_panel_bic.do                     # Hao-Kasahara panel BIC
    ├── phase2_unobs_het.do                     # Unobserved heterogeneity
    ├── phase2_ccp_estimation.do                # CCP estimation
    ├── phase2_dynamic_ccp.do                   # Dynamic CCP
    ├── phase3_counterfactual_final.do          # K=4 counterfactual analysis
    ├── phase3_counterfactuals.do               # Counterfactual variants
    ├── phase3_structural_cf.do                 # Structural counterfactual
    │
    ├── # ML ROBUSTNESS (Phase 5)
    ├── phase5_lasso_reduced_form.do            # Post-double-selection LASSO
    ├── phase5_causal_forest.do                 # Causal forest
    ├── phase5_ml_ccps.do                       # ML-based CCPs
    ├── phase5_structural_bounded.do            # Bounded structural
    │
    ├── # EXTENSIONS (Phases 6-7)
    ├── phase6_cohort_pseudopanel.do            # Pseudo-panel
    ├── phase6_static_structural.do             # Static structural
    ├── phase7_mixed_logit.do                   # Mixed logit (continuous het)
    ├── phase7_mte.do                           # Marginal treatment effects
    ├── phase7_dml.do                           # Double ML
    ├── phase7_cps_panel_linkage.do             # CPS panel linkage
    │
    ├── # BAYESIAN & SENSITIVITY (Phase 8)
    ├── phase8_bayesian_dpm.R                   # Bayesian DPM (Stan)
    ├── phase8_bayesian_dpm_sensitivity.R       # DPM sensitivity
    ├── phase8_conformal_inference.do           # Conformal inference
    ├── phase8_sensitivity_analysis.do          # Specification sensitivity
    ├── phase8_distributional_synth.do          # Distributional synthetic
    ├── phase8_tastenet_mnl.py                  # Neural net MNL
    │
    ├── # MONTE CARLO
    ├── monte_carlo_simulation.R                # Two-design MC (cancellation vs same-sign)
    │
    ├── # FIGURES
    ├── create_methods_paper_figures.py          # All 14 paper figures
    ├── create_paper1_figures.py                 # Companion paper figures
    │
    └── output/
        ├── figures_methods/                    # 14 figures (PDF + PNG)
        ├── monte_carlo_*.csv                   # MC results
        ├── phase2_*.csv                        # Model selection results
        ├── phase3_*.csv                        # Counterfactual results
        └── phase8_dpm_*.csv                    # Bayesian results
```

## Completion Status

### DONE (logs closed successfully)
- [x] phase0_reduced_form.do — OLS benchmark, tables 1-3, IV (CBSA FE fails, state FE in continuation)
- [x] phase0_continuation.do — IV (state FE), heterogeneity, CBSA-level, CCPs
- [x] phase1_completion.do — base MNL estimation
- [x] phase1_structural_mlogit.do — 9-alternative MNL, margins, conditional logit
- [x] phase2_type_selection.do — K=1–4 BIC comparison (K=4 selected)
- [x] phase2_panel_bic.do — Hao-Kasahara panel BIC (K=4 confirmed)
- [x] phase2_unobs_het.do — unobserved heterogeneity estimation
- [x] phase2_ccp_estimation.do — CCP estimation with structural parameters
- [x] phase2_dynamic_ccp.do — dynamic CCP estimation
- [x] phase3_simple.do — simple counterfactual
- [x] phase3_consistent_cf.do — consistent counterfactual
- [x] phase3_structural_cf.do — structural counterfactual
- [x] phase3_counterfactuals.do — counterfactual variants (K=1–4)
- [x] phase3_counterfactual_final.do — K=4 counterfactual with LPM and delta method CIs
- [x] phase5_lasso_reduced_form.do — post-double-selection LASSO
- [x] phase5_causal_forest.do — causal forest heterogeneity
- [x] phase5_ml_ccps.do — ML-based CCPs
- [x] phase5_structural_bounded.do — bounded structural estimation
- [x] phase6_cohort_pseudopanel.do — pseudo-panel analysis
- [x] phase6_static_structural.do — static structural estimation
- [x] phase7_mixed_logit.do — mixed logit (continuous heterogeneity)
- [x] phase7_mte.do — marginal treatment effects
- [x] phase7_dml.do — double ML + AIPW
- [x] phase7_cps_panel_linkage.do — CPS panel linkage
- [x] phase8_conformal_inference.do — conformal inference
- [x] phase8_sensitivity_analysis.do — Oster sensitivity analysis
- [x] phase8_distributional_synth.do — distributional synthetic control
- [x] phase8_bayesian_dpm_sensitivity.R — DPM sensitivity (20 subsamples)
- [x] monte_carlo_simulation.R — full MC (cancellation + same-sign, N=2K/10K/50K)
- [x] create_methods_paper_figures.py — all 14 figures generated
- [x] Methods Paper Seeking Perspective.tex — compiled (Mar 26)
- [x] extract_demographics.do — extracted PESEX, PEMARITL, PRNMCHLD from raw CPS into analysis dataset

### RUN (with reduced sample for verification)
- [x] phase8_bayesian_dpm.R — Stan DPM model verified with N=500 (K_eff=2.6, pipeline works). Production run (N=10K, 4 chains) needs cluster (~5+ hrs).

### RE-RUN with real demographics (Mar 26)
8 do-files re-run after adding female, married, has_children from raw CPS:
- [x] phase5_lasso_reduced_form.do
- [x] phase5_causal_forest.do
- [x] phase5_ml_ccps.do
- [x] phase7_mte.do
- [x] phase7_mixed_logit.do
- [x] phase7_dml.do
- [x] phase8_sensitivity_analysis.do
- [x] phase8_conformal_inference.do

## Key Verified Results (from output CSVs — audited Mar 30)
- Panel BIC selects K=4 (both standard and HK)
- Standard BIC: K=1(-78728), K=2(-79073), K=3(-80700), K=4(-82110)
- Panel BIC: K=1(-78810), K=2(-79169), K=3(-80811), K=4(-82235)
- OSCE: 5 active types
- Baseline SE rate: 10.56%
- Type shares: 12.1%, 24.1%, 33.0%, 30.8%
- Type betas: −0.037, −0.020, 0.010, 0.093
- 50% closure counterfactual (K=4): −10.8%, 95% CI [−14.1%, −7.4%]
- BLM counterfactual stability: K=1→−0.57%, K=2→−2.49%, K=3→−7.32%, K=4→−9.75%
- MC cancellation DGP: true K=4 recovered 3% (N=2K), 8% (N=10K) of replications
- MC CF range: median ~15–17pp across K
- Bayesian DPM: K_eff≈3.8, P(K≥4)≈0.79, CF≈−8.5%, P(neg)≈0.99
- Demographics: Female SE rate 7.7% vs Male 11.8%; Married SE rate 11.9% vs Unmarried 8.0%
- Gender heterogeneity: mobile×female=0.0253*** (SE=0.0072), joint F=6.21 (p=0.002)

## Dataset Variables
Variables in `analysis_dataset_with_se.dta` (after demographics merge):
- **Demographics**: age, pesex/female, pemaritl/married, prnmchld/has_children, praceeth3, peducgrp, hhincome (categorical 1-5)
- **Employment**: self_employed, wage_worker, employed, unemployed, emp_status
- **Banking**: banking_mode (1=unbanked, 2=mobile, 3=branch), mobile_user, branch_user, banked, unbanked
- **Geography**: cbsa, statefips, metro, pct_broadband
- **Survey**: year (2013-2023), weight, hsupwgtk
- **Household**: hhtype, hhtypev2, hhtypev3

## Development Rules
1. DO NOT overwrite existing do-files — create new versions or separate files
2. Preserve all log files and output CSVs
3. Use Stata 17+ for .do files, R 4.x for .R files, Python 3.10+ for .py files
4. All paths use the OneDrive base: `/Users/amalkova/Library/CloudStorage/OneDrive-FloridaInstituteofTechnology/_Research/Mobile_Money_Banking/Mobile banking USA/`
5. Figures go to `Data/output/figures_methods/` in both PDF and PNG
6. Survey weights: `hsupwgtk`; clustering: CBSA level
