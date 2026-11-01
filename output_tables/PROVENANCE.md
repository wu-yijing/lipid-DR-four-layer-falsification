# Aggregate result tables: provenance register

Every table under `output_tables/` is listed below with the script that produces it, the
analysis variant it belongs to, and which part of the manuscript it supports. This
register exists because the `layer*` / `MR_*` file names are **not** the names the
scripts write: the scripts emit phase-named files into `phase*_out/` in the analysis
workspace, and this repository holds a documented selection of those outputs under
stable names. `tools/artefact_map.json` is the machine-readable form and
`tools/make_provenance.py` regenerates this file from it.

**Run `tools/audit_repo_consistency.py` before quoting any number from this directory.**
The audit re-derives every quantity the manuscript reports and fails if a committed
table stops supporting it.

## License

The tables in this directory are aggregate results and are released under the
Creative Commons Attribution 4.0 International licence (CC BY 4.0); the terms are
in the `LICENSE` file beside this register. The analysis code is released under the
MIT licence, in the `LICENSE` file at the repository root. Manuscript text and
individual-level data are not covered by either licence.

## Provenance classes

| class | n | meaning |
|---|---|---|
| `direct` | 51 | a committed script writes this exact name into `output_tables/` |
| `collected` | 47 | byte-identical to a named output of a named script (rename step made explicit) |
| `hand-assembled` | 11 | built by hand from script outputs; values re-derived by the audit |
| `local-only` | 8 | record-level file held locally only; listed in `.gitignore`, never committed |

## Superseded or conflicting artefacts (do not quote)

### `layer4_increment_NRI_IDI.csv`

full-model out-of-fold AUC 0.7354 and dAUC 0.2494, from the earlier phase-4 fold split; the current values are 0.7331 and 0.2470 in layer4_three_row_AUC.csv. The NRI and IDI it reports (0.410 / 0.152) are therefore computed on a different dAUC basis than the dAUC quoted alongside them. Either recompute NRI/IDI on the current split or state the provenance when citing 0.410 / 0.152.

### `layer2_outcome_rebuild_crosstab.csv`

SUPERSEDED and actively misleading if read as counts. Its cells are the binary comparison 'source text states DR / does not' at an earlier pipeline stage, with analysed_as_non_case n = 2,207 rather than the current 1,691. The labels analysed_as_case = 271 + 4 and analysed_as_non_case = 1,752 + 455 were read in an earlier draft of the manuscript as '4 cases had an explicit diagnosis and 455 diagnosed records were analysed as non-cases'. That reading is WRONG. The current source of truth is layer2_outcome_crosstab_full.csv: 474 of the 475 records carrying the miscoded value are true DR cases, and none of the 275 records analysed as cases carries a gradeable DR diagnosis (50 are grade 0, 225 are not gradeable).

## Rows that are known to be defective (partial)

These files are cited by the manuscript and are correct except for the rows named
below. `tools/audit_repo_consistency.py` asserts that each of these flags is present,
so the defect cannot be forgotten or silently repaired out of the register.

### `layer4_three_row_AUC.csv` — row `local_lipid_only_std`

The `local_lipid_only_std` row reports n = 888 and AUC_unweighted = 0.4861. Both are artefacts of an alignment defect in the `oof()` helper of phase4e_local_v2.R: `model.matrix()` drops rows with missing predictors, so for the extended lipid model (terms tg_log and ldl_mmol, missing in 473 of the 888 reference records) the fold-wise replacement vector was shorter than `p[fold == k]` and R recycled it with a warning rather than failing. The quantity genuinely exists on 415 records and equals 0.4988 (95% CI 0.4368-0.5575). The other three rows are unaffected and still reproduce exactly. Corrected values: layer4_three_row_AUC_corrected.csv, computed by 4_external/phase10_local_comparator_increment.R, which asserts the reproduction of the three unaffected rows and asserts that the corrected std row is not 888.

**Manuscript impact:** Section 3.4's sentence about the extended variable set must drop 'AUC 0.486' and either quote 0.4988 on 415 records or remove the sentence. The conclusion - the lipid-only model performs no better - is unchanged.

### `layer4_dedup_three_row_AUC.csv` — row `local_lipid_only_std`

The same alignment defect, in the repeat-screening sensitivity run (phase4e_local_v2_dedup.R). Not corrected here because no manuscript statement depends on this row; re-running the patched script regenerates it correctly.

**Manuscript impact:** none (no manuscript statement cites this row)

## Hand-assembled summaries

These have no producing script. They are the highest-risk class in a data-integrity
screen, so each is stated explicitly.

| file | assembled from |
|---|---|
| `layer2_cohort_repeat_screening.csv` | transcribed from the console output of 3_prediction/build_dedup_devset.py |
| `layer2_dedup_cohort_and_outcome.csv` | assembled by hand from the layer2_dedup_* outputs |
| `layer2_dedup_discrimination_and_missingness.csv` | assembled by hand: primary and repeat-screening values placed side by side |
| `layer2_dedup_screening_selection_rule.csv` | assembled by hand from 3_prediction/build_dedup_devset.py and the dedup outcome comparison |
| `layer2_outcome_rebuild_crosstab.csv` | assembled by hand at an earlier pipeline stage (see SUPERSEDED) |
| `layer2_outcome_rebuild_exclusions.csv` | assembled by hand from severity_rebuild_audit.csv |
| `layer2_outcome_rebuild_eye_categories.csv` | assembled by hand from severity_rebuild_audit.csv |
| `layer2_outcome_rebuild_patient_grades.csv` | assembled by hand from severity_rebuild_audit.csv (3_prediction/rebuild_outcome.py) |
| `layer2_pmsampsize_grid_fine.csv` | assembled by hand from the pmsampsize routine in 3_prediction/phase2e_tail2.R |
| `layer2_text_mode_5level_matrix_r2r3_adj.csv` | assembled by hand from 3_prediction/phase9c_text_mode.py (adjudicated matrix) |
| `layer34_dedup_transport_and_benchmark.csv` | assembled by hand from layer3_* and layer4_* outputs; every value in it is independently re-derived by tools/audit_repo_consistency.py, which is why it is safe to cite |

## Written directly by a committed script

| file | variant | producing script / source | supports |
|---|---|---|---|
| `MR_exmore_sensitivity.csv` | primary | `1_mr\phase8_exmore_sensitivity.py` | Section 3.1 / Supplementary Table S12 |
| `MVMR_core_TG_LDL_HDL.csv` | primary | `1_mr/revision_20261030_mvmr_joint_set.R` | Table 2 MVMR block |
| `MVMR_joint_set_diagnostics.csv` | primary | `1_mr/revision_20261030_mvmr_joint_set.R` |  |
| `NHANES_temporal_validation.csv` | primary | `4_external/revision_20261030_temporal_primary.R` |  |
| `cis_drug_target_MR_PPARG.csv` | primary | `2_cis_mr\phase18_pparg_cis_mr.R` |  |
| `cis_instrument_strength_PPARG.csv` | primary | `2_cis_mr\phase18_pparg_cis_mr.R` |  |
| `cis_pparg_correlated_ivw.csv` | primary | `2_cis_mr\phase18_pparg_cis_mr.R` |  |
| `cis_pparg_ld_clump_crosscheck.csv` | primary | `2_cis_mr\phase19_pparg_ld_clump_server.R` |  |
| `cis_pparg_per_instrument_wald.csv` | primary | `2_cis_mr\phase18_pparg_cis_mr.R` |  |
| `cis_six_target_ld_clump_crosscheck.csv` | primary | `2_cis_mr\phase19b_six_target_ld_clump_crosscheck.R` |  |
| `dca_fixed_operating.csv` | primary | `4_external/revision_20261101_decision_analysis.py` |  |
| `dca_net_benefit_curve.csv` | primary | `4_external/revision_20261101_decision_analysis.py` |  |
| `dca_operating_points.csv` | primary | `4_external/revision_20261101_decision_analysis.py` |  |
| `layer1_MDE_by_exposure.csv` | primary | `1_mr\phase12_power_mde.py` |  |
| `layer1_cis_target_MDE.csv` | primary | `2_cis_mr\phase13_cis_MDE.py` |  |
| `layer1_cis_target_MDE_PPARG.csv` | primary | `2_cis_mr\phase18_pparg_cis_mr.R` |  |
| `layer1_external_appendicitis_negctrl.csv` | primary | `4_external\phase15_external_outcome_mr.R` |  |
| `layer1_external_height_negctrl.csv` | primary | `4_external\phase15_external_outcome_mr.R` |  |
| `layer1_external_ukb_dr_replication.csv` | primary | `4_external\phase15_external_outcome_mr.R` |  |
| `layer1_instrument_detail_F_R2.csv` | primary | `1_mr\phase12_power_mde.py` |  |
| `layer1_instrument_strength_by_exposure.csv` | primary | `1_mr\phase12_power_mde.py` |  |
| `layer1_negctrl_estimates_appendicitis_negctrl.csv` | primary | `4_external\phase17_negctrl_diagnostics.R` |  |
| `layer1_negctrl_heterogeneity_appendicitis_negctrl.csv` | primary | `4_external\phase17_negctrl_diagnostics.R` |  |
| `layer1_negctrl_ivw_fe_re_appendicitis_negctrl.csv` | primary | `4_external\phase17_negctrl_diagnostics.R` |  |
| `layer1_negctrl_pleiotropy_appendicitis_negctrl.csv` | primary | `4_external\phase17_negctrl_diagnostics.R` |  |
| `layer1_robust_estimators.csv` | primary | `1_mr\phase14_robust_estimators.R` |  |
| `layer1_winners_curse_threshold_sensitivity.csv` | primary | `1_mr\phase12_power_mde.py` |  |
| `layer2_design_condition_numbers.csv` | primary | `3_prediction/phase10_reporting_gaps.py` | Section 3.2 - basis of the condition number |
| `layer2_independent_reextraction_summary.csv` | primary | `3_prediction\phase9_independent_reextraction.py` | Supplementary Table S14a |
| `layer2_linkage_window.csv` | primary | `3_prediction/phase22_linkage_window.R` |  |
| `layer2_missingness_patterns.csv` | primary | `3_prediction\phase8_review_revision_analyses.py` | Supplementary Table S10 |
| `layer2_outcome_audit_and_bounds.csv` | primary | `3_prediction\phase8b_outcome_error_bounds.py` | Supplementary Table S13 |
| `layer2_outcome_crosstab_full.csv` | primary | `3_prediction\phase8_review_revision_analyses.py` | Methods 2.4 / Supplementary Table S11 |
| `layer2_parsimonious_markers_classified.csv` | primary | `3_prediction/phase10_reporting_gaps.py` | Supplementary Table S9 - explicit marker classes |
| `layer2_parsimonious_models.csv` | primary | `3_prediction\phase8_review_revision_analyses.py` | Supplementary Table S9 |
| `layer2_pmsampsize_verification.csv` | primary | `3_prediction/phase10_reporting_gaps.py` | limitation 6 - verification of the quoted requirements |
| `layer2_readset_subsets.csv` | primary | `3_prediction\phase9d_aprime.py` |  |
| `layer2_readset_subsets_summary.json` | primary | `3_prediction\phase9d_aprime.py` |  |
| `layer2_reporting_gap_inventory.csv` | primary | `3_prediction/phase10_reporting_gaps.py` | audit trail - statistics still without an interval |
| `layer2_reporting_gaps.csv` | primary | `3_prediction/phase10_reporting_gaps.py` | Section 3.2 / Table 3 - intervals for AUC and IDI |
| `layer2_second_grader_composition.csv` | primary | `3_prediction/phase10_reporting_gaps.py` | limitation 12 - composition of the review sample |
| `layer2_second_grader_sampling_summary.json` | primary | `3_prediction/phase9b_second_grader_sampling.py` | limitation 12 re-abstraction plan |
| `layer2_text_mode_agreement_r2r3_adj.json` | primary | `3_prediction/phase9c_text_mode.py` | Supplementary Table S14a |
| `layer3_paired_AUC_CI.csv` | primary | `4_external\phase21_paired_auc_ci.R` | Table 4 confidence intervals (paired bootstrap; supersedes the unpaired auc_ci() of phase4_external_v2.R, see REVISION_NOTES_20261011.md) |
| `layer4_lipid_increment_negative_control.csv` | primary | `4_external/phase11_lipid_increment_over_full_model.R` | Section 3.4 / Supplementary Table S21 - same-number random-variable control |
| `layer4_lipid_increment_over_full_model.csv` | primary | `4_external/phase11_lipid_increment_over_full_model.R` | Section 3.4 / Supplementary Table S21 - lipids added to the full reference model |
| `layer4_local_comparator_increment.csv` | primary | `4_external/phase10_local_comparator_increment.R` | Section 3.4 - attributable increment with intervals |
| `layer4_three_row_AUC_corrected.csv` | primary | `4_external/phase10_local_comparator_increment.R` | Section 3.4 - corrected three-row comparison |
| `phase20_outcome_ungradable_bounds.csv` | primary | `4_external\phase20_outcome_ungradable_bounds.R` |  |
| `tableS1_primary_frame.csv` | primary | `4_external/revision_20261031_tableS1_primary_frame.R` |  |
| `tableS42_primary_frame.csv` | primary | `4_external/revision_20261031_tableS42_primary_frame.R` |  |

## Collected and renamed from a script output

| file | variant | producing script / source | supports |
|---|---|---|---|
| `MR_Egger_intercepts.csv` | primary | `<SRC>\phase1_out\mr_egger_intercept.csv` | Table 2 note / Supplementary Table S3 |
| `MR_PRESSO_10k_summary.csv` | primary | `<SRC>\phase1_out\presso10k_summary.csv` | Table 2 MR-PRESSO row |
| `MR_heterogeneity_CochranQ.csv` | primary | `<SRC>\phase1_out\mr_heterogeneity.csv` | Table 2 note / Supplementary Table S3 |
| `MR_two_sample_all_methods.csv` | primary | `<SRC>\phase1_out\mr_main.csv` | Table 2 (all methods) |
| `NHANES_outcome_conventions.csv` | primary | `<SRC>\phase4_out_v2d\phase20_outcome_ungradable_bounds.csv` |  |
| `cis_drug_target_MR_six_targets.csv` | primary | `<SRC>\phase1_out\druggable_mr_results.csv` | Figure 2B / Supplementary Table S4 |
| `cis_instrument_strength_by_target.csv` | primary | `<SRC>\phase2_out\phase2f_cis_strength.csv` | Methods 2.3 instrument strength |
| `cis_marginal_F_TG_LDL_HDL.csv` | primary | `<SRC>\phase2_out\phase2f_marginal_f.csv` | Methods 2.3 conditional/marginal F |
| `emerging_markers_NHANES.csv` | primary | `<SRC>\phase4_out\phase4c_nhanes_rc_tyg.csv`, `<SRC>\phase4_out_v2\phase4c_nhanes_rc_tyg.csv` | Section 3.5 / Supplementary Table S1 |
| `emerging_markers_chinese_cohort.csv` | primary | `<SRC>\phase4_out_v2\phase4c_chinese_rc.csv` | Section 3.5 / Supplementary Table S1 |
| `instrument_mean_F_by_exposure.csv` | primary | `<SRC>\phase1_out\iv_mean_fstat.csv` | Methods 2.3 mean F by exposure |
| `layer2_MI_pooled_coefficients_3seeds.csv` | primary | `<SRC>\phase2_out\phase_v2_mi_coef.csv` | Supplementary Table S7 |
| `layer2_complete_case_coefficients_VIF.csv` | primary | `<SRC>\phase2_out\phase_v2_cc_model.csv` | Table 3 coefficients + VIF |
| `layer2_complete_case_performance.csv` | primary | `<SRC>\phase2_out_v2\phase2_performance.csv` | Table 3 performance block |
| `layer2_dedup_complete_case_coefficients_VIF.csv` | sensitivity | `<SRC>\phase2_out_v2d\phase_v2d_cc_model.csv` | limitation 12 |
| `layer2_dedup_complete_case_performance.csv` | sensitivity | `<SRC>\phase2_out_v2d\phase2_performance.csv` | limitation 12 |
| `layer2_dedup_delta_auc_bootstrap.csv` | sensitivity | `<SRC>\phase2_out_v2d\phase_v2d_deltaauc.csv` | limitation 12 (+0.005) |
| `layer2_dedup_missingness_rates.csv` | sensitivity | `<SRC>\phase2_out_v2d\phase_v2d_missing_rates.csv` | limitation 12 |
| `layer2_dedup_pmsampsize_grid.csv` | sensitivity | `<SRC>\phase2_out_v2d\phase_v2d_pmsampsize.csv` | limitation 12 |
| `layer2_dedup_selection_model.csv` | sensitivity | `<SRC>\phase2_out_v2d\phase_v2d_selection_model.csv` | limitation 12 |
| `layer2_dedup_three_row_AUC_decomposition.csv` | sensitivity | `<SRC>\phase2_out_v2d\phase_v2d_decomposition.csv` | limitation 12 |
| `layer2_delta_auc_bootstrap.csv` | primary | `<SRC>\phase2_out\phase_v2_deltaauc.csv` | Section 3.2 bootstrap difference |
| `layer2_fraction_missing_information.csv` | primary | `<SRC>\phase7_out_v2\imp_fmi.csv` | Supplementary Figure S1 / Table S6 |
| `layer2_imputation_observed_vs_imputed.csv` | primary | `<SRC>\phase7_out_v2\imp_diag_stats.csv` | Supplementary Figure S1 / Table S6 |
| `layer2_missingness_mechanism.csv` | primary | `<SRC>\phase7_out_v2\missing_assoc_or.csv` | Methods 2.4 missingness mechanism |
| `layer2_missingness_rates.csv` | primary | `<SRC>\phase2_out\phase_v2_missing_rates.csv`, `<SRC>\phase7_out_v2\missing_rates.csv` | Methods 2.4 missingness rates |
| `layer2_pmsampsize_grid.csv` | primary | `<SRC>\phase2_out\phase_v2_pmsampsize.csv` | Methods 2.4 / limitation 6 sample size |
| `layer2_selection_model.csv` | primary | `<SRC>\phase2_out\phase_v2_selection_model.csv` | Section 3.2 completeness model (OR 1.44) |
| `layer2_severity_coding_table.csv` | primary | `<SRC>\phase2_out\phase2f_severity_table.csv` | Methods 2.4 severity coding |
| `layer2_three_row_AUC_decomposition.csv` | primary | `<SRC>\phase2_out\phase_v2_decomposition.csv` | Section 3.2 three-row decomposition |
| `layer3_MI_pooled_export_sensitivity.csv` | primary | `<SRC>\phase4_out_v2\phase4d_mi_export.csv` | Section 3.3 MI-pooled export (0.517 / 0.518) |
| `layer3_chinese_coefficients_by_cycle.csv` | primary | `<SRC>\phase4_out\v2_phase4b_chinese_by_cycle.csv` | Section 3.3 temporal subsets |
| `layer3_dedup_external_validation_AUC.csv` | sensitivity | `<SRC>\phase4_out_v2d\phase4_external_auc.csv` | limitation 12 (core 0.490 / std 0.485) |
| `layer3_dedup_weighted_sensitivity_MEC.csv` | sensitivity | `<SRC>\phase7_out_v2d\v2d_phase7h_weighted_ext.csv` | limitation 12 (slopes -1.65 / -3.19) |
| `layer3_external_validation_AUC.csv` | primary | `<SRC>\phase4_out_v2\phase4_external_auc.csv` | Table 4 core / std / LASSO rows |
| `layer3_recomputed_AUC_CI_B2000.csv` | primary | `<SRC>\phase4_out_v2\phase4e_auc_ci.csv` | Table 4 confidence intervals |
| `layer3_weighted_sensitivity_MEC.csv` | primary | `<SRC>\phase7_out\v2_phase7h_weighted_ext.csv` | Table 4 MEC-weighted columns |
| `layer4_dedup_three_row_AUC.csv` | sensitivity | `<SRC>\phase4_out_v2d\phase4e_three_row.csv` | limitation 12 (0.218 preserved) |
| `layer4_dedup_weighted_calibration.csv` | sensitivity | `<SRC>\phase7_out_v2d\v2d_phase7h_calib_weighted.csv` | limitation 12 |
| `layer4_dedup_weighted_full_vs_lipid.csv` | sensitivity | `<SRC>\phase7_out_v2d\v2d_phase7h_full_vs_lipid_weighted.csv` | limitation 12 |
| `layer4_full_reference_model.csv` | primary | `<SRC>\phase4_out\phase4_reference_model.csv`, `<SRC>\phase4_out_v2\phase4_reference_model.csv`, `<SRC>\phase4_out_v2d\phase4_reference_model.csv` | Table 5 |
| `layer4_increment_NRI_IDI.csv` | legacy (superseded) | `<SRC>\phase4_out_v2\phase4_increment.csv` |  |
| `layer4_temporal_validation.csv` | primary | `<SRC>\phase4_out\phase4b_temporal.csv`, `<SRC>\phase4_out\v2_phase4b_temporal.csv` | Section 3.4 temporal validation |
| `layer4_three_row_AUC.csv` | primary | `<SRC>\phase4_out_v2\phase4e_three_row.csv` | Section 3.4 / Figure 4 three-row comparison |
| `layer4_three_row_AUC_primary.csv` | primary | `<SRC>\..\_r1_review_20261029\pipeB\out_B\layer4_three_row_AUC_corrected.csv` |  |
| `layer4_weighted_calibration.csv` | primary | `<SRC>\phase7_out\v2_phase7h_calib_weighted.csv` | Section 3.4 weighted calibration |
| `layer4_weighted_full_vs_lipid.csv` | primary | `<SRC>\phase7_out\v2_phase7h_full_vs_lipid_weighted.csv` | Section 3.4 weighted re-check + Brier |

## Hand-assembled summaries

| file | variant | producing script / source | supports |
|---|---|---|---|
| `layer2_cohort_repeat_screening.csv` | primary | _(none)_ | repeat-screening accounting |
| `layer2_dedup_cohort_and_outcome.csv` | sensitivity | _(none)_ | limitation 12 sensitivity cohort |
| `layer2_dedup_discrimination_and_missingness.csv` | sensitivity | _(none)_ | limitation 12 summary |
| `layer2_dedup_screening_selection_rule.csv` | sensitivity | _(none)_ | limitation 12 selection rule |
| `layer2_outcome_rebuild_crosstab.csv` | legacy (superseded) | _(none)_ |  |
| `layer2_outcome_rebuild_exclusions.csv` | primary | _(none)_ | Methods 2.4 / limitation 12 exclusions |
| `layer2_outcome_rebuild_eye_categories.csv` | primary | _(none)_ | Methods 2.4 eye-level categories |
| `layer2_outcome_rebuild_patient_grades.csv` | primary | _(none)_ | Methods 2.4 grade distribution |
| `layer2_pmsampsize_grid_fine.csv` | primary | _(none)_ | limitation 6 grid (C = 0.56-0.70) |
| `layer2_text_mode_5level_matrix_r2r3_adj.csv` | primary | _(none)_ | Supplementary Table S14a (adjudicated) |
| `layer34_dedup_transport_and_benchmark.csv` | primary | _(none)_ | Table 1 / README transport + benchmark summary |

## Held locally only (never committed)

These carry record-level identifiers or source-report text and are excluded by
`.gitignore`. They are listed here so that the gap between the local workspace and the
repository is explicit rather than implied.

| file | producing script |
|---|---|
| `layer2_independent_reextraction.csv` | `3_prediction/phase9_independent_reextraction.py` |
| `layer2_independent_reextraction_discordances.csv` | `3_prediction/phase9_independent_reextraction.py` |
| `layer2_independent_reextraction_discordances_preadjudication.csv` | `3_prediction/phase9_independent_reextraction.py` |
| `layer2_second_grader_sampling.csv` | `3_prediction/phase9b_second_grader_sampling.py` |
| `layer2_text_mode_5level_matrix__selftest_returnedxlsx.csv` | `3_prediction/phase9c_text_mode.py` |
| `layer2_text_mode_agreement__selftest_returnedxlsx.json` | `3_prediction/phase9c_text_mode.py` |
| `layer2_text_mode_discordances.csv` | `3_prediction/phase9c_text_mode.py` |
| `layer2_text_mode_practice_key.csv` | `3_prediction/phase9c_text_mode.py` |

## How to cite a variant in the manuscript

| manuscript statement | correct file |
|---|---|
| Table 1 / layer-3 transport and layer-4 benchmark summary | `layer34_dedup_transport_and_benchmark.csv` |
| Table 4 core / std / LASSO unweighted AUC | `layer3_external_validation_AUC.csv` |
| Table 4 MEC-weighted AUC and calibration slopes | `layer3_weighted_sensitivity_MEC.csv` |
| Section 3.4 three-row comparison and the 0.218 / 0.029 decomposition | `layer4_three_row_AUC.csv` |
| Section 3.4 weighted re-check, Brier 0.176 / 0.203, slope -1.52 | `layer4_weighted_full_vs_lipid.csv`, `layer4_weighted_calibration.csv` |
| Section 3.2 complete-case AUC 0.598 / 0.579 and calibration 1.000 / 0.000 | `layer2_complete_case_performance.csv` |
| Section 3.2 multiply-imputed AUC 0.582 and the CC-MI comparison | `layer2_dedup_discrimination_and_missingness.csv` |
| Methods 2.4 grade distribution 1,497 / 49 / 386 / 15 / 24 | `layer2_outcome_rebuild_patient_grades.csv` |
| Methods 2.4 exclusions 511 = 386 + 53 + 72 | `layer2_outcome_rebuild_exclusions.csv` |
| Methods 2.4 / Table S11 record-level cross-tabulation | `layer2_outcome_crosstab_full.csv` |

`_dedup` files are the **repeat-screening sensitivity analysis** (one record per patient)
and must never be quoted as if they were the primary analysis. Note that the primary and
`_dedup` layer-3 files carry near-identical confidence intervals but different point
estimates and very different calibration slopes; always cite the file name alongside the
number.
