# Publisher screening: what the policies ask for, and where this repository answers

This repository is submitted alongside a manuscript that contains a Mendelian randomisation
layer and two NHANES-based layers. Both of those categories are now subject to explicit
editorial screening at the publisher the manuscript targets. This file maps each condition to
the concrete artefact that satisfies it, so that the claim of compliance is checkable rather
than asserted.

Nothing here is a promise of acceptance. It is the evidence package.

---

## The policies, as published

| date | publisher | condition | source |
|---|---|---|---|
| 2024-06/07 | Frontiers | independent external validation required for **all health-dataset MR studies**; MR papers without new validation or institutional data declined at desk review. 5,513 MR submissions rejected since. | [Frontiers, *Cutting through fast-churn science*](https://www.frontiersin.org/news/2025/09/15/cutting-through-fast-churn-science-how-frontiers-raised-the-bar) |
| 2024-04 | PLOS | completed **STROBE-MR checklist required** before an academic editor sees an MR submission; 40% of all submissions rejected at this pre-screen | [Pruller, Mockridge & Vousden, *PLoS One* 2025;20(11):e0337199](https://doi.org/10.1371/journal.pone.0337199) |
| 2025-01/05 | Frontiers | external-validation requirement extended to **all papers based on simple queries of public data, including NHANES**; authors must supply new experimental validation **or additional data from their own institution**. 1,382 NHANES-based submissions rejected since May 2025. | as above |
| 2025-09 | PLOS | updated standards for public health/social-science databases; papers without additional validation rejected automatically. Rejection rate for this class rose from **40% to 94%**. | [PLOS editorial policy update](https://theplosblog.plos.org/2025/09/updates-to-plos-retrospective-health-database-editorial-policy/) |

Supporting commentary: *Science* (doi:10.1126/science.aed1020), *BMJ* 2025;391:r2170.

Frontiers' own wording for the accepted route is explicit: *"authors could use NHANES data
alongside private data from their own institution"* (Elena Vicario, Head of Research Integrity).

---

## Condition-by-condition mapping

| # | screening condition | satisfied by | artefact |
|---|---|---|---|
| 1 | Institutional (non-public) data must be part of the study — public data alone is not enough | Layered 2 is built entirely on the authors' own hospital EMR cohort (2,482 screening records from 2,333 identifiers), held under ethics approval BBYY2023072. The public resources play supporting roles only. | `layer2_*`; manuscript §2.2, §2.4 |
| 2 | A testable hypothesis, not pattern-hunting | Four numbered hypotheses with explicit refutation conditions (H1a/H1b/H1c, H2', H3, H4); the outcome definition was pre-specified before the corrected pipeline was run. | manuscript §1.3, §2.1; `layer34_dedup_transport_and_benchmark.csv` |
| 3 | The role of each public dataset must be stated, and it must not be the source of hypothesis-generating correlations | NHANES is used for exactly two things: external validation of the institutional coefficients (layer 3), and a locally fitted comparator benchmark that makes the transport test interpretable (layer 4). The comparison that carries no transport assumption is the one fitted locally. | `README.md` §"Data roles"; `layer4_three_row_AUC_corrected.csv` |
| 4 | Findings robust to the choice of statistical specification | Layer 1: univariable IVW, MVMR, MR-PRESSO, six cis drug-target analyses, and an EXMORE control-exclusion sensitivity. Layer 2: complete cases versus multiple imputation, inverse-probability-of-selection weighting, a three-row decomposition, and a repeat-screening sensitivity. Every contrast is reported, including the ones that do not flatter the conclusion. | `MR_*.csv`, `MVMR_*.csv`, `cis_*.csv`, `layer2_delta_auc_bootstrap.csv`, `layer2_dedup_*` |
| 5 | Completed STROBE-MR checklist | supplied as a separate submission file | `Checklists_TRIPOD_STROBE_STROBE-MR.docx` |
| 6 | Completed TRIPOD checklist (prediction-model layer) | supplied in the same file | as above |
| 7 | Uncertainty must accompany performance estimates | every headline AUC, IDI and increment metric now carries a 95% interval, with the basis of each interval named | `layer2_reporting_gaps.csv`, `layer4_local_comparator_increment.csv` |
| 8 | Code availability, and code that actually reproduces the numbers | the aggregate register states, per file, which script produces it and which analysis variant it belongs to; a self-audit re-derives every headline quantity from the committed tables and exits non-zero on any failure | `output_tables/PROVENANCE.md`, `tools/make_provenance.py`, `tools/audit_repo_consistency.py` |
| 9 | No unexplained divergence between an artefact in the repository and a number in the paper | the audit gates the manuscript directly: `tools/audit_repo_consistency.py --manuscript path.docx` checks the manuscript for the headline counts and for known-superseded phrasings | `tools/audit_repo_consistency.py --help` |
| 10 | Discussion of gene-environment equivalence | stated explicitly that the estimates are scaled to genetically predicted lipid differences and are not dose-equivalent to a clinical intervention, and that the drug-target results are evidence on the direction of genetic support rather than predictions of trial effects | manuscript §2.3, §4.3 |
| 11 | Results framed as falsification rather than as a claim of causality | the paper claims no causal effect for any factor; its central result is a negative one, and the constructive claim is explicitly marked as corroborative rather than novel | manuscript §4.4, §5 |

---

## The two analyses that are still outstanding

These are the only conditions this repository does **not** yet satisfy, and they are the reason
they are listed here rather than buried in a limitations paragraph.

### A. Outcome-side positive control

**What it is.** Genetically predicted type 2 diabetes, or HbA1c, against the same DR endpoint.
**Why it matters now.** Layer 1's robustness currently rests entirely on internal sensitivity
analyses. EXMORE is a different endpoint *definition within the same FinnGen release*, not a
second dataset, so as it stands "no effect" cannot be separated from "no detectable effect".
A positive control is the direct answer, and it is also the cleanest reply to the charge that
this is an automated pipeline application.
**Script.** `1_mr/phase11_positive_control.py --mode positive`
**Blocked on.** an instrument panel for the positive-control exposure
(`<SRC>/phase1_out/positive_control_instruments.csv`) and the layer-1 primary outcome file
(`LDF_OUTCOME_GZ`). `--check` reports both.

### B. Replication in a second, independent DR outcome dataset

**What it is.** The same lipid instruments against an independent diabetic-retinopathy GWAS.
**Why it matters now.** This is the requirement the policies state in terms. The GWAS Catalog
summary statistics are fully public and require no token, so this is achievable without new
data access permissions.
**Script.** `1_mr/phase11_positive_control.py --mode replicate`
**Blocked on.** a pre-extracted outcome subset (`LDF_OUTCOME_SUBSET`). The lipid instrument
panel is already present (`phase1_out/harmonised_data.csv`).

`phase11` refuses to emit a table when its inputs are absent, and returns exit status 3 with a
report rather than writing a stub. The code path itself has been exercised end to end on the
EXMORE subset purely as a smoke test; the table it produced was deleted, because EXMORE is not
an independent dataset and must not be presented as a replication.

---

## Two things to say explicitly in the cover letter

1. **The public-dataset role, in one sentence, before anything else.** The manuscript carries a
   screening-relevant signal on its first page — "NHANES" appears in the abstract twice — and the
   institutional cohort is not mentioned by name until §2.2. State in the cover letter that the
   hypotheses were developed and tested in the authors' own hospital cohort, and that NHANES is
   used strictly as an independent external validation set and as a locally fitted comparator.
2. **The contribution, against the existing literature, in a table.** Both PLOS One and the
   *Journal of Global Health* checklist format now ask what has already been done on the question
   and what is new. The prior lipid–DR MR studies and DR prediction models are already cited in
   the manuscript; presenting them as a one-page comparison removes the burden of proof.

## What must not be claimed

- Do not describe the EXMORE analysis as external validation. It is a sensitivity analysis on a
  differently defined control set within the same release.
- Do not quote the `local_lipid_only_std` row of `layer4_three_row_AUC.csv`. It is defective; use
  `layer4_three_row_AUC_corrected.csv` (n = 415, AUC 0.4988). See `output_tables/PROVENANCE.md`.
- Do not quote the `_dedup` files as the primary analysis; they are the repeat-screening
  sensitivity analysis.
- Do not quote NRI 0.410 from `layer4_increment_NRI_IDI.csv`. It comes from an earlier fold split;
  the same contrast on the current split is NRI 0.386 (unweighted) or 0.421 (MEC-weighted), with
  intervals, in `layer4_local_comparator_increment.csv`.
