# Revision notes, 1 November 2026 (second pass on the v6 package)

Two additions were made after the v6 package was assembled; neither changes a point
estimate. Both answer Reviewer 1 point 4, where the manuscript claimed a corrected
survey design but tabulated only record-level intervals for discrimination.

## Design-based intervals for the three layer-4 AUCs

`4_external/revision_20261101_psu_auc_bootstrap.R`

- frame: ungradable photographs as a missing outcome (Methods 2.2), n = 756, 267 events,
  identical fold seed and model specification as `revision_20261030_figure4_redraw.R`;
  the point estimates reproduce `layer4_three_row_AUC_primary.csv` exactly
  (0.7856 / 0.5397 / 0.4674 unweighted; 0.8039 weighted for the full model);
- weighted AUC: PSU subbootstrap replicate weights (`as.svrepdesign`, B = 1000);
- unweighted AUC: PSU cluster bootstrap, PSUs resampled within strata (B = 1000);
- output: `layer4_auc_design_ci.csv`, `layer4_delta_auc_design_ci.csv`, written to
  `LDF_OUT` (default: a session temp directory).

| model | AUC (unweighted) | PSU cluster 95% CI | AUC (weighted) | design 95% CI |
|---|---|---|---|---|
| full reference | 0.7856 | 0.7562–0.8157 | 0.8039 | 0.7599–0.8385 |
| local lipid-only | 0.5397 | 0.5137–0.5682 | 0.5529 | 0.5046–0.6046 |
| transported Chinese | 0.4674 | 0.4352–0.4990 | 0.4506 | 0.4031–0.4916 |

delta-AUC, paired PSU resamples: full − lipid 0.2459 (0.2066–0.2876);
full − transported 0.3182 (0.2747–0.3609); lipid − transported 0.0723 (0.0213–0.1238).

The design-based intervals are narrower than the archived record-level ones
(0.751–0.819, 0.498–0.580, 0.424–0.512) by 10–20%. The two resample different units and
the archived run is not in this repository, so both are reported in the manuscript rather
than one being substituted for the other: the record-level interval stays in the text and
the design-based interval is given beside it. Anyone quoting a single interval should use
the design-based one, which is the one the script below reproduces.

## Decision-first ordering in Discussion 4.1

The opening of 4.1 now leads with the decision-level result (net benefit at or below
"refer everyone"; 97% and 100% referred at a threshold of 0.20; 22–30% of cases detected
at a 25% burden against 51%) and presents the 0.70 margin afterwards as the formal
criterion for that same decision, so that the argument no longer depends on the margin.
The decision numbers it quotes are those of Supplementary Table S48 and Section 3.4.
