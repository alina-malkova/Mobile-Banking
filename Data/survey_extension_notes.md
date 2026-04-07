# Survey Extension: Cancellation Classification (Comments 1 & 7)

## Current Survey Scope
- 40 papers, 2018-2024, 7 top journals (Econometrica, AER, QJE, REStud, JPE, QE, JBES)
- Keywords: "finite mixture," "latent class," "unobserved types," "number of components"
- Finding: only 6/40 (15%) report CF sensitivity to K

## Referee Request
- **Expand keywords**: add "random coefficients," "discrete heterogeneity," "grouped fixed effects"
- **Classify each paper** by cancellation plausibility
- Create appendix table with classification

## Cancellation Classification Framework

**Cancellation plausible** = the treatment/policy variable plausibly has opposite-signed effects across latent types. This can happen when:

1. **Selection vs. treatment effects**: Some types benefit from treatment, others are harmed (e.g., branch access helps credit-constrained entrepreneurs but enables over-leveraging for others)
2. **Heterogeneous comparative advantage**: Policy shifts composition — beneficial for some types' comparative advantage, harmful for others
3. **Offsetting behavioral responses**: Subsidy increases one type's participation but crowds out another type
4. **Sign-switching coefficients**: Type-specific parameters have literal opposite signs

**Cancellation unlikely** = all types plausibly respond in the same direction (possibly different magnitudes)

---

## Candidate Papers for Expanded Survey

### A. ALREADY IN BIBLIOGRAPHY (classify these first)

| # | Paper | Journal | Year | Field | Cancellation? |
|---|-------|---------|------|-------|---------------|
| 1 | Arcidiacono & Miller — CCP with unobserved het | Econometrica | 2011 | Labor | Possible |
| 2 | Bonhomme & Manresa — Grouped fixed effects | Econometrica | 2015 | Methods/Labor | Possible |
| 3 | Bonhomme, Lamadon, Manresa — Distributional framework | Econometrica | 2019 | Labor | Possible |
| 4 | Bonhomme, Lamadon, Manresa — Discretizing unobs het | Econometrica | 2022 | Methods | Possible |
| 5 | Fox, Kim, Ryan, Bajari — Random coefficients | QE | 2011 | IO/Demand | Likely |
| 6 | Keane & Wolpin — Career decisions | JPE | 1997 | Labor | Possible |
| 7 | Heckman & Singer — Duration data | Econometrica | 1984 | Labor | Possible |
| 8 | Christensen & Connault — CF sensitivity & robustness | Econometrica | 2023 | Methods | Likely |
| 9 | Budanova — OSCE for finite mixtures | J Econometrics | 2025 | Methods | N/A |
| 10 | Kasahara & Shimotsu — Identification of FMM | Econometrica | 2009 | Methods | N/A |

### B. NEW PAPERS TO ADD (from web search, 2018-2025)

#### B1. Applied papers with finite mixtures and counterfactuals

| # | Paper | Journal | Year | Field | Why relevant |
|---|-------|---------|------|-------|-------------|
| 11 | **Arcidiacono, Aucejo, Maurel, Ransom** — College attrition dynamics | JPE | 2025 | Education | Dynamic discrete choice with types; dropout CF |
| 12 | **Kalouptsidi, Kitamura, Scott** — ID of CFs in dynamic DC | QE | 2021 | IO | Identification of CFs with permanent types |
| 13 | **Saltiel** — Dynamic treatment effects of job training | WP | 2025 | Labor | Latent types + training counterfactuals |
| 14 | **Aguirregabiria, Iaria, Sokullu** — Demand with endogenous entry | WP | 2023 | IO | Finite mixture demand + counterfactuals |

#### B2. Random coefficients papers (new keyword)

| # | Paper | Journal | Year | Field | Why relevant |
|---|-------|---------|------|-------|-------------|
| 15 | **BLP-style demand** — various applied IO papers | Various | 2018-24 | IO | Random coefficients in demand; price effects could differ in sign across consumer types |
| 16 | **QE 2025** — ID of random coefficient latent utility | QE | 2025 | Methods | Random coefficients identification + welfare |

#### B3. Grouped fixed effects papers (new keyword)

| # | Paper | Journal | Year | Field | Why relevant |
|---|-------|---------|------|-------|-------------|
| 17 | **Various GFE applications** — labor, macro, finance | Various | 2019-25 | Various | Following Bonhomme-Manresa; group-specific treatment effects |
| 18 | **Time-varying group unobs het** — finance | JBES | 2025 | Finance | GFE in finance with time-varying types |

#### B4. Counterfactual sensitivity papers

| # | Paper | Journal | Year | Field | Why relevant |
|---|-------|---------|------|-------|-------------|
| 19 | **Christensen & Connault** — CF sensitivity | Econometrica | 2023 | Methods | Bounds on CFs under distributional perturbation |
| 20 | **Chen** — Robust structural estimation under misspecified latent dynamics | WP | 2025 | Methods | Sensitivity of welfare to latent state dynamics |

---

## Specific Search Queries for Author to Run on Google Scholar/JSTOR

The following searches should be run to find the remaining ~20-30 papers needed to reach 40-60 total:

### Query 1: "finite mixture" + counterfactual
```
("finite mixture" OR "latent class") AND (counterfactual OR "treatment effect" OR "policy simulation")
site:aeaweb.org OR site:onlinelibrary.wiley.com/journal/14680262
```
**Expected yield**: 10-15 applied papers

### Query 2: "random coefficients" + applied
```
"random coefficients" AND (counterfactual OR welfare OR "policy evaluation")
journal:(AER OR "American Economic Review" OR Econometrica OR "Review of Economic Studies" OR "Journal of Political Economy")
2018-2024
```
**Expected yield**: 5-10 demand/IO papers using BLP-type random coefficients

### Query 3: "discrete heterogeneity" OR "grouped fixed effects" + applied
```
("discrete heterogeneity" OR "grouped fixed effects" OR "latent groups") AND (counterfactual OR "treatment effect")
2018-2024
```
**Expected yield**: 5-10 papers following Bonhomme-Manresa

### Query 4: "unobserved types" + structural
```
"unobserved types" AND structural AND (counterfactual OR "policy simulation")
2018-2024
```
**Expected yield**: 5-8 dynamic discrete choice papers

---

## Classification Protocol

For each paper, determine:

1. **K_used**: What K did they use? (e.g., K=3, K=2-5)
2. **K_selection_method**: How was K chosen? (BIC, AIC, LRT, cross-validation, ad hoc, not discussed)
3. **reports_cf_by_K**: Does the paper show how the counterfactual changes with K? (Yes/Partial/No)
4. **cancellation_plausible**: Could the treatment effects plausibly have opposite signs across types? (Likely/Possible/Unlikely)
5. **cancellation_reasoning**: Brief justification (1-2 sentences)

## Key Insight for Response Letter

The expanded survey should support the argument that:
- Cancellation is not exotic — it's plausible in many applied settings (labor sorting, demand heterogeneity, education returns)
- Most papers don't check for it (only 15% report CF by K)
- The "when K matters" diagnostic (checking for sign changes across types) is cheap and should be standard practice
