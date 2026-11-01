# v1.1.0 — revision analyses, a decision analysis, and an audit record

Released 1 November 2026. Supersedes [v1.0.0](https://doi.org/10.5281/zenodo.23331502)
(12 October 2026).

This release adds the code and aggregate outputs for the analyses carried out during peer
review, plus the documentation that records what changed and how it was checked. The study,
the data sources and the conclusions are the same; several headline estimates were recomputed
under a corrected outcome convention, and every one of them is traceable in
`REVISION_NOTES_20261030.md`, `…_20261031.md` and `…_20261101.md`.

## Added — analyses

| analysis | script |
|---|---|
| multivariable MR on a jointly clumped instrument set, with conditional strength and heterogeneity diagnostics | `1_mr/revision_20261030_mvmr_joint_set.R` |
| NHANES survey design: strata, PSUs, four-year and fasting-subsample weights, design-based intervals | `4_external/revision_20261030_nhanes_design.R` |
| decision analysis: net benefit, operating points, referral burden | `4_external/revision_20261101_decision_analysis.py` |
| design-based intervals for the three benchmark models (PSU subbootstrap and cluster bootstrap) | `4_external/revision_20261101_psu_auc_bootstrap.R` |
| temporal validation under the primary outcome definition | `4_external/revision_20261030_temporal_primary.R` |
| bounds across the three conventions for ungradable photographs | `4_external/phase20_outcome_ungradable_bounds.R` |
| paired AUC intervals | `4_external/phase21_paired_auc_ci.R` |
| cis-arm LD-clumping cross-check (OpenGWAS service against local PLINK 2) | `2_cis_mr/phase19_pparg_ld_clump_server.R`, `2_cis_mr/phase19b_six_target_ld_clump_crosscheck.R` |

## Added — aggregate tables

`dca_net_benefit_curve.csv`, `dca_operating_points.csv`, `dca_fixed_operating.csv`,
`NHANES_outcome_conventions.csv`, `NHANES_temporal_validation.csv`,
`layer4_three_row_AUC_primary.csv`, `layer2_linkage_window.csv`,
`tableS1_primary_frame.csv`, `tableS42_primary_frame.csv`,
`MVMR_joint_set_diagnostics.csv`, `cis_six_target_ld_clump_crosscheck.csv`,
`cis_pparg_ld_clump_crosscheck.csv`.

All are de-identified aggregates. `output_tables/PROVENANCE.md` records, per file, the script
that writes it and the analysis variant it belongs to.

## Added — documentation

- `docs/audit_20261101.md` — what was compared, how, the findings and how each was dispositioned,
  and what remains open. No manuscript text and no correspondence.
- `REVISION_NOTES_20261030.md`, `…_20261031.md`, `…_20261101.md` — the revision rounds, per change.

## Changed

- Figure 1 redrawn so that each layer carries its role in the argument, and so that the layer-2
  verdict describes coefficient instability without implying that stability is a criterion
  (`figures/figure1_roles.py`).
- Figure 3 and Figure 4 regenerated on the primary outcome frame.
- The licence is split and recorded: MIT for code, CC BY 4.0 for the aggregate tables
  (`output_tables/LICENSE`).
- `CITATION.cff` now points at the concept DOI, so the file stays correct across versions;
  version DOIs are listed as identifiers.

## Not in this release

No manuscript text, no response to reviewers, no individual-level clinical records, no
reviewer correspondence. `.gitignore` enforces the first of these and the rest were never
committed.

## What to read first

1. `README.md` — data roles and the four-layer framework.
2. `docs/audit_20261101.md` — what the revision changed and how it was checked.
3. `output_tables/PROVENANCE.md` — which script produces which table.

## Citation

Concept DOI (all versions, always latest): **10.5281/zenodo.23331501**.
Version DOI for this release: **[10.5281/zenodo.23857272](https://doi.org/10.5281/zenodo.23857272)**.
