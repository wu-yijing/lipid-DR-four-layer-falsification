# Revision notes, 2 November 2026 (exposure-association covariance in the multivariable MR)

Reviewer 1, point 6 asked for three things: the jointly clumped union of instruments, the
conditional strength and heterogeneity diagnostics, and **the exposure-association covariance
matrix** — either estimated, or replaced by a justified assumption. The v6 revision supplied
the first two and, for the third, supplied an argument rather than the matrix. That is the
gap this round closes.

`1_mr/revision_20261102_mvmr_covariance.R`

## Why the matrix cannot simply be handed to an existing estimator

Both packaged implementations accept a covariance argument and then do not use it for
inference, which is worth recording because a reply that says "the covariance was supplied"
would otherwise be untrue in substance:

| call | what `gencov` / `correl.x` actually does |
|---|---|
| `MVMR::ivw_mvmr(r_input, gencov)` | nothing. The body is `lm(betaYG ~ -1 + exposures, weights = 1/seY^2)`; `gencov` is never read. |
| `MendelianRandomization::mr_mvivw(correl.x = ...)` | enters `condFstat()` only. The estimator's variance is computed from `correl`/`rho`, which is the **between-variant** LD matrix, not the exposure covariance. |

`MVMR::strength_mvmr()` does use it, but only as a list of one 3 × 3 matrix per variant;
passing a plain matrix trips `if (!is.list(gencov) && gencov == 0)`, and passing the
per-variant matrix in the wrong scale produces conditional F statistics below 1. Both
failure modes are easy to hit and neither is loud.

## The matrix and the model

The three exposures are measured in the same GLGC participants, so their association
estimates at a given variant are correlated. Across the 133 jointly clumped variants the
correlation matrix of the per-variant exposure effects is

|  | TG | LDL-C | HDL-C |
|---|---|---|---|
| TG | 1 | −0.0707 | −0.2915 |
| LDL-C | −0.0707 | 1 | −0.2289 |
| HDL-C | −0.2915 | −0.2289 | 1 |

and it enters the model as Σ_j = R ∘ s_j s_j′, the Hadamard product of R with the outer
product of the three standard errors at variant j.

The point estimate is weighted least squares and does not involve Σ, so it cannot move. The
variance does:

    Var(θ̂) = s² A⁻¹ + A⁻¹ [ Σ_j w_j² (θ̂′ Σ_j θ̂) Σ_j ] A⁻¹ ,   A = X′ W X

with w_j = 1/se_Yj² and s² the weighted residual variance. Setting Σ = 0 kills the second
term and reduces the expression to s² A⁻¹, which is what `ivw_mvmr()` reports — the script
checks this and reproduces 0.055587 / 0.032479 / 0.042957 to six decimals before adding the
covariance term.

## Result

| specification | conditional F (TG / LDL-C / HDL-C) | TG OR (95% CI) |
|---|---|---|
| exposures treated as independent (Σ = 0) | 42.0 / 82.1 / 60.8 | 0.927 (0.831–1.034) |
| covariance supplied (R as estimated) | 54.0 / 109.4 / 84.5 | 0.927 (0.831–1.034) |

LDL-C 1.006 (0.944–1.072) and HDL-C 1.041 (0.957–1.132) are likewise unchanged. The two
standard errors differ by 5.3 × 10⁻⁶ relative: with 133 instruments the exposure-uncertainty
term is negligible beside the residual term.

Conditional F **rises**, by 29% to 39%. That is the expected direction when the exposures are
negatively correlated — partialling the others out is more informative, not less — and it
means the attenuation of the triglyceride estimate does not rest on an unexamined
independence assumption. Sensitivity over R × 0.5, R × 2 and all off-diagonals set to ±0.60
gives conditional F between 27.0 and 108.0 for TG and leaves every estimate and standard
error identical.

## Files

| file | content |
|---|---|
| `1_mr/revision_20261102_mvmr_covariance.R` | the analysis; takes the joint-set table as its first argument |
| `output_tables/mvmr_exposure_correlation.csv` | R, the matrix above |
| `output_tables/mvmr_covariance_adjusted.csv` | conditional estimates and standard errors under both specifications |
| `output_tables/mvmr_covariance_sensitivity.csv` | the sensitivity grid, conditional F and standard errors per cell |

The joint-set table itself (`1_mr/union_joint_set.csv`, 133 variants × 8 columns of public
GLGC 2013 and FinnGen R13 summary statistics) is shipped, through an explicit `!`
exception in `.gitignore`, so that the script runs from a clean clone without
re-downloading the source GWAS. `*.csv` remains excluded everywhere else; the table
contains no individual-level records or identifiers.

## What changed in the manuscript

Methods 2.3 now reports the matrix, the covariance, both sets of conditional F statistics and
the invariance of the estimates; Table 2's caption carries both sets of conditional F;
Discussion 4.5 registers the covariance rather than an independence assumption; Supplementary
Table S49 tabulates the matrix, both specifications and the sensitivity grid; the response to
point 6 and the Q2 checklist reply are updated to match.
