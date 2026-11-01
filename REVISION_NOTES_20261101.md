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

---

## Third pass (v6.4): three P1 and five P2 items found by a version-to-version audit

Comparing the submitted manuscript with the v6.3 revision line by line turned up eight
defects that the earlier passes had introduced or left half-closed. None changes a point
estimate, a confidence interval or a conclusion; all are consistency, attribution or
exposition. Three would have been visible to a reviewer who checked the specific place the
corresponding comment pointed at, which is why they were treated as blocking.

| # | where | what was wrong | fix |
|---|---|---|---|
| 27 | response letter, header | carried a date one day before the changes it described | dated to the submission date |
| 28 | response letter, point 12 | quoted 0.467 as the remnant-cholesterol AUC | 0.436, the value in Results 3.5; 0.467 is the transported score's AUC |
| 29 | Discussion 4.5, endpoint scope | the minimum detectable odds ratio was given per 1% HbA1c while the observed upper limit beside it was given per SD | both now per SD, with the per-1% equivalents in parentheses; this is the scale the power calculation of Methods 2.3 uses |
| 30 | Conclusions; Results 3.2 heading; Figure 1 | coefficient *interpretability* was still written as a failure condition, after Table 1 and Supplementary Table S47 had removed stability from the criterion | "unstable", reported but explicitly not a criterion; Figure 1 regenerated from `figures/figure1_roles.py` |
| 31 | Discussion 4.2, closing | the word cut had deleted the sentence that resolved the inventory of discordant reports, leaving the paragraph without a conclusion | one sentence restored: opposite directions for the same exposure under different instrument sets read as instrument dependence, and null results being under-reported means the record over-represents whichever set was significant |
| 32 | Table 5 header | column labelled `p (BH-adjusted)` while the text reports q | relabelled `q (BH-adjusted)` |
| 33 | Discussion 4.1, first paragraph | Supplementary Table S41, a layer-2 robustness table, was cited for the export clause | S41 attached to the source-cohort clause, `Table 4` cited for the export clause |
| 34 | Methods 2.9 | the script producing the three design-based AUC intervals was missing from the list of revision-time scripts | `4_external/revision_20261101_psu_auc_bootstrap.R` added |
| 35 | Discussion 4.1 | stray space before a full stop | removed |

Word budget: the fixes add 54 words, funded by six compressions that remove no estimate,
conclusion or citation — the height rationale in Results 3.1.3 (already stated in
Discussion 4.5), the ApoA1/ApoB and automated-module sentences in 4.5, the
alternative-modalities sentence in 4.4, the freeze sentence in 2.9, and the opening of 4.5.
The declared main-text count moves from 11,979 to 11,978 words.

The three marked copies (clean, tracked, redline) were edited together. Accepting every
revision in the tracked file and removing every struck run from the redline both reproduce
the clean text paragraph for paragraph (510 non-empty paragraphs), which is the same
invariant the earlier passes were checked against.

### What this repository does and does not contain

The manuscript, the response letter, the reviewer correspondence and the packaging manifest
stay out of this repository: `*.docx` and the manuscript-assembly scripts are excluded by
`.gitignore`, and the data-protection rules there also keep individual-level records and
reader packs out of `output_tables/` (see `POLICY_COMPLIANCE.md`). What is here is the
analysis code, the de-identified aggregate tables, and these notes, which describe each
change without reproducing manuscript prose. If a claim in the notes and a table in
`output_tables/` disagree, the table is the authority.

### Layer-2 calibration: stated, not quoted

The layer-2 model is fitted and evaluated on the same records, so its apparent calibration
slope and intercept are 1.00 and 0.00 by construction. Discussion 4.5 now says that
calibration is not reported for this model and why, instead of leaving the gap unexplained.
`3_prediction/revision_20261102_layer2_calibration.R` computes the optimism-corrected slope
and intercept by bootstrap; it is **provided unrun**, because it needs the individual-level
development set, which is not deposited here. No table in `output_tables/` corresponds to
it, and none is claimed.
