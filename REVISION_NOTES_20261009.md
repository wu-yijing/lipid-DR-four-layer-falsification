# Change record — 2026-10-09

Audit-driven changes to the analysis code and its provenance, made in response to an internal
consistency audit of the manuscript and an audit of the publisher screening conditions it must
now pass. No scientific conclusion changes; one manuscript sentence does.

**Pre-change commit:** `6846791` (working tree clean)
**Backup of every file that was edited:** `_bak_before_reportfix_20261009/`
**Roll back everything:** `git checkout 6846791 -- .` and delete the new files listed below.

---

## 1. Substantive fix: silent misalignment in out-of-fold prediction

**Defect.** `oof()` in `4_external/phase4e_local_v2.R` built fold-wise predictions with

```r
p[fold == k] <- plogis(as.numeric(model.matrix(fmla, data = te) %*% as.numeric(coef(f))))
```

`model.matrix()` applies `na.action` and drops incomplete rows. For a model whose predictors
contain missing values the replacement vector is then shorter than `p[fold == k]`, and R recycles
it with a warning instead of failing. In layer 4 this affects only the extended lipid model, whose
terms `tg_log` and `ldl_mmol` are missing in 473 of the 888 reference records.

**Evidence.** `sum(nrow(te))` = 888 versus `sum(nrow(model.matrix(fmla, data = te)))` = 415, with
one length-mismatch warning per fold.

**Consequence.** The committed row `local_lipid_only_std` in `layer4_three_row_AUC.csv` reported
`n = 888` and `AUC = 0.4861`. The quantity exists on **415** records and equals **0.4988**
(95% CI 0.4368-0.5575). The other three rows are unaffected and reproduce exactly.

**Fix.** `oof()` now predicts on the complete rows of each fold and writes the values back to their
own positions, with `stopifnot(length(pr) == sum(okk))` as a guard so the failure mode cannot
recur silently. Applied to:

- `4_external/phase4e_local_v2.R`
- `4_external/phase4e_local_v2_dedup.R`
- `4_external/phase10_local_comparator_increment.R` (new)

**Artefacts.** `output_tables/layer4_three_row_AUC_corrected.csv` holds the corrected table. The
committed file is left in place and flagged, row by row, in `output_tables/PROVENANCE.md`, because
it is cited by the manuscript and should not change under the author's feet.

**Manuscript action.** §3.4's sentence about the extended variable set quotes "AUC 0.486" on an
implied n of 888. Replace with 0.4988 on 415 records, or remove the sentence. The conclusion it
supports — the lipid-only model performs no better — is unchanged.

---

## 2. Reporting gaps closed

New script `3_prediction/phase10_reporting_gaps.py` (runs in ~15 s) supplies:

| quantity | previously | now |
|---|---|---|
| layer-2 complete-case AUC, apparent | 0.598, no interval | 0.5984 (0.5748-0.6398) |
| layer-2 complete-case AUC, optimism-corrected | 0.579, no interval | 0.5793 (0.5702-0.5977) |
| layer-2 lipid IDI | 0.020, no interval | 0.0199 (0.0110-0.0436) |
| collinearity condition number | 8.5, basis unstated | **8.54 column-standardised**, against 75.01 raw and 73.00 for the correlation matrix |
| parsimonious marker classes | ambiguous count (read as 12 models) | explicit `marker_class`, 4 conventional lipids + 2 apolipoproteins + 4 derived indices = 10 |
| second-grader review sample | 251 stated, composition not closable | 251 = strata sum; random stratum **40**, excluded records **26** |
| sample-size requirements | quoted, unverified | every quoted value reproduced from the committed grid, break-even C = 0.63 confirmed |

The optimism-corrected value reproduces the committed 0.5794 and the IDI the committed 0.0199, so
the new intervals are attached to the same estimates the paper already reports.

---

## 3. Increment metrics for the comparison the paper can actually attribute

New script `4_external/phase10_local_comparator_increment.R` (runs in ~9 s). Section 3.4 argues
that the full-model-versus-transported-score contrast cannot be attributed to either cause; that
contrast was the only one carrying NRI and IDI. The attributable contrast now has them, with
intervals:

| contrast | n | dAUC (95% CI) | NRI (95% CI) | IDI (95% CI) |
|---|---|---|---|---|
| full vs locally fitted lipid-only (CORE) | 888 | 0.218 (0.165, 0.266) | 0.375 (0.290, 0.463) | 0.150 (0.120, 0.179) |
| full vs transported Chinese score | 888 | 0.247 (0.189, 0.304) | 0.386 (0.322, 0.449) | 0.154 (0.127, 0.183) |
| locally fitted (CORE) vs transported | 888 | 0.029 (-0.041, 0.096) | 0.060 (-0.022, 0.136) | 0.004 (-0.006, 0.013) |

Two consequences for the manuscript.

- The transport component (0.029) is **indistinguishable from zero**, while the information
  component (0.218) is clearly not. The paper's argument is strengthened by saying so explicitly.
- The quoted NRI 0.410 is from an earlier fold split. On the current split the same contrast gives
  0.386 unweighted (0.421 MEC-weighted). Quote the current value with its weighting, or state the
  provenance.

The script also asserts that the three unaffected three-row AUCs still reproduce, so it doubles as
a regression check.

---

## 4. Provenance made mechanical

The README previously claimed that "the file names in `output_tables/` indicate which script
produces which aggregate output". That was not true: only 8 of 81 files are written under those
names by a committed script; 46 are byte-identical to phase-named outputs, 11 are hand-assembled
(summary tables with no producing script), and 8 are record-level files held locally only.

- `tools/artefact_map.json` — machine-readable register, 81/81 files classified.
- `tools/make_provenance.py` — regenerates `output_tables/PROVENANCE.md`, and `--check` fails if a
  file is unclassified.
- `output_tables/PROVENANCE.md` — per-file register: variant (primary / repeat-screening
  sensitivity / legacy), producing script, and which manuscript element it supports. Flags two
  superseded artefacts and two partially defective files.

The superseded artefacts:

- `layer4_increment_NRI_IDI.csv` — full-model AUC 0.7354 and dAUC 0.2494 from the earlier phase-4
  split; NRI/IDI therefore sit on a different dAUC basis than the dAUC quoted beside them.
- `layer2_outcome_rebuild_crosstab.csv` — its cells are **labels, not counts**, and reading them as
  counts produced the "a handful of cases / several hundred misclassified non-cases" account that
  appeared in an earlier draft of the README. The record-level truth is in
  `layer2_outcome_crosstab_full.csv`: 474 of the 475 miscoded records are true DR cases, and none
  of the 275 records analysed as cases carries a gradeable DR diagnosis.

---

## 5. Self-audit

`tools/audit_repo_consistency.py` — 179 checks, exits non-zero on any failure.

- recomputes every arithmetic claim from the committed tables (cohort accounting, the
  inclusion-exclusion count, the AUC decomposition, events per variable, the Bonferroni threshold,
  the missingness pattern mean);
- checks **every MR table for internal consistency** using the identity that an
  inverse-variance-weighted estimate is the geometric mean of its limits, `OR == sqrt(lo*hi)`;
- gates the README against superseded narratives and internal duplication;
- verifies the alignment guard is present and the defective idiom is gone;
- enforces a **disclosure gate**: no git-tracked aggregate file may contain a record identifier or
  Chinese source-report text;
- optionally cross-checks the manuscript: `--manuscript path.docx`, reporting five manuscript-side
  items that remain open.

`1_mr/phase8_exmore_sensitivity.py` and `3_prediction/phase8_review_revision_analyses.py` gained
columns, not changed values: `ci_p_consistent` (true for all 12 EXMORE rows, which documents that
an earlier audit query about the TG row was unfounded) and `marker_class`.

---

## 6. Publisher-screening package

`POLICY_COMPLIANCE.md` maps each published screening condition — Frontiers' external-validation
requirement for health-dataset MR and for public-data papers, PLOS's STROBE-MR pre-screen and its
2025 database policy, and the reporting items named in the *European Journal of Epidemiology*
editorial — to the concrete artefact that answers it, and names the two analyses still outstanding
with runnable scripts: an outcome-side positive control and replication in a second independent DR
outcome dataset (`1_mr/phase11_positive_control.py`). That script refuses to emit a table when its
inputs are absent, returning exit status 3 with a report of how to obtain each missing file.

---

## 7. New files

```
tools/paths.py                              shared input resolution (LDF_SRC), replaces hard-coded paths
tools/mr_estimators.py                      IVW / MR-Egger / weighted median / heterogeneity, shared
tools/artefact_map.json                     per-file provenance register (81 entries)
tools/make_provenance.py                    regenerates output_tables/PROVENANCE.md
tools/audit_repo_consistency.py             the self-audit (179 checks)
3_prediction/phase10_reporting_gaps.py      intervals, condition-number bases, marker classes
4_external/phase10_local_comparator_increment.R   attributable increments + the alignment fix
1_mr/phase11_positive_control.py            data-gated positive control and outcome replication
POLICY_COMPLIANCE.md                        screening conditions mapped to artefacts
output_tables/PROVENANCE.md                 generated register
```

Modified: `README.md` (rewritten), `4_external/phase4e_local_v2.R`,
`4_external/phase4e_local_v2_dedup.R`, `1_mr/phase8_exmore_sensitivity.py`,
`3_prediction/phase8_review_revision_analyses.py`.

---

## 8. How to verify

```bash
python tools/paths.py                                   # active input configuration
python tools/audit_repo_consistency.py                  # 179 checks, exit 0
python tools/make_provenance.py                         # regenerate + re-check the register
python 3_prediction/phase10_reporting_gaps.py           # the new intervals
Rscript 4_external/phase10_local_comparator_increment.R # corrected three-row table + increments
python 1_mr/phase11_positive_control.py --check         # exit 3 here: inputs are not distributed
```
