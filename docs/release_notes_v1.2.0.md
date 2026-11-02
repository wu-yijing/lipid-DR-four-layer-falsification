# v1.2.0 — exposure-association covariance in the multivariable MR

Released 2 November 2026. Supersedes [v1.1.0](https://doi.org/10.5281/zenodo.23857272)
(1 November 2026).

This release adds one analysis and the data it needs. No estimate in the manuscript moves:
the conditional odds ratios for triglyceride, LDL-C and HDL-C are unchanged to four decimal
places, and their standard errors change by five parts per million. What changes is that the
model no longer assumes the three exposures are independent when it reports their conditional
strength.

## Added — analysis

| analysis | script |
|---|---|
| exposure-association covariance in the multivariable MR, with conditional F statistics under both specifications and a sensitivity grid over the assumed correlation | `1_mr/revision_20261102_mvmr_covariance.R` |

Reviewer 1 asked for the covariance matrix of the exposure associations, not for an argument
that it can be ignored. Neither packaged estimator can carry it:

- `MVMR::ivw_mvmr(r_input, gencov)` accepts `gencov` and never reads it — the body is a
  weighted linear model, and the reported standard errors come from that model;
- `MendelianRandomization::mr_mvivw(correl.x = ...)` uses `correl.x` only inside
  `condFstat()`; its `correl`/`rho` argument is the between-variant LD matrix, a different
  object entirely.

`MVMR::strength_mvmr()` does consume it, but only as a list holding one 3 × 3 matrix per
variant and only at the scale of standard errors squared; passing the correlation-times-SE
matrix instead yields conditional F statistics below 1, without an error.

The covariance-adjusted variance is therefore computed directly by M-estimation and validated
against the published standard errors by setting the covariance to zero:

    Var(theta) = s^2 A^-1 + A^-1 [ sum_j w_j^2 (theta' Sigma_j theta) Sigma_j ] A^-1

with Σ_j = R ∘ s_j s_j′ and R the correlation matrix of the three exposure estimates across the
133 jointly clumped variants (TG–LDL-C −0.071, TG–HDL-C −0.292, LDL-C–HDL-C −0.229). At
Σ = 0 the second term vanishes and the expression reduces to what `ivw_mvmr()` reports, which
the script checks to six decimals before adding it.

| specification | conditional F (TG / LDL-C / HDL-C) |
|---|---|
| exposures treated as independent | 42.0 / 82.1 / 60.8 |
| covariance supplied | 54.0 / 109.4 / 84.5 |

Conditional strength rises by 29% to 39%, the expected direction for negatively correlated
exposures, so the attenuation of the triglyceride estimate no longer rests on an unexamined
independence assumption.

## Added — data

| file | content |
|---|---|
| `output_tables/mvmr_exposure_correlation.csv` | R, the matrix above |
| `output_tables/mvmr_covariance_adjusted.csv` | conditional estimates and standard errors under both specifications |
| `output_tables/mvmr_covariance_sensitivity.csv` | the sensitivity grid, conditional F and standard errors per cell |
| `1_mr/union_joint_set.csv` | the 133-variant joint set the script reads: rsid, three exposure effect/SE pairs, the outcome effect/SE and a palindromic flag; public GLGC 2013 and FinnGen R13 summary statistics, no individual-level records |

The joint-set table is shipped through an explicit `!` exception in `.gitignore` so that the
script runs from a clean clone without re-downloading the source GWAS; `*.csv` remains
excluded everywhere else.

## Documentation

`REVISION_NOTES_20261102.md` records the derivation, the two failure modes of the packaged
estimators, and the self-check that ties the covariance-adjusted variance back to the
published one.

No manuscript text, no response to reviewers, no reviewer correspondence and no
individual-level clinical records are included.

## What to read first

1. `REVISION_NOTES_20261102.md` — why the matrix is needed and what it changes.
2. `README.md` — data roles and the four-layer framework.
3. `output_tables/PROVENANCE.md` — which script produces which table.

## Citation

Concept DOI (all versions, always latest): **10.5281/zenodo.23331501**.
Version DOI for this release: **10.5281/zenodo.23877907**
([doi.org/10.5281/zenodo.23877907](https://doi.org/10.5281/zenodo.23877907)).

The archived snapshot was cut before this DOI existed, so `CITATION.cff` inside the zip lists the
v1.1.0 DOI as superseded but cannot carry the v1.2.0 DOI itself; the repository commit that follows
the deposit adds it. This is inherent to a version DOI that is only reserved at publication.
