# Revision round 1 (30 October 2026)

Changes made in response to the first reviewer report, and the code that produces them.
The manuscript text and the supplementary material are not part of this repository; only
code and aggregate tables are.

## What changed in the numbers

| Quantity | Before | After |
|---|---|---|
| MVMR, mutually adjusted TG | 0.94 (0.85-1.04), 41 variants | **0.93 (0.83-1.03)**, 133 variants |
| MVMR, LDL-C / HDL-C | 1.01 / 1.04, 61 / 75 variants | 1.01 (0.94-1.07) / 1.04 (0.96-1.13), 133 variants each |
| Conditional F (TG / LDL-C / HDL-C) | not reported | 42.0 / 82.1 / 60.8 |
| Q_A | not reported | 188.8 on 130 df, p = 5.8e-4 |
| Layer-4 reference model | 0.733, n = 888 | 0.786 (0.804 weighted), n = 756 |
| Layer-3 transported score | 0.486, n = 888 | 0.467 (0.424-0.512), n = 756 |
| Temporal validation (full model) | 0.748 / 0.722 | 0.784 / 0.775 |

The layer-3/4 change follows from treating ungradable fundus photographs as a missing
outcome rather than as non-cases (variant C of `NHANES_outcome_conventions.csv`). The
MVMR change follows from refitting the mutually adjusted model on one instrument set.

## Added code

| File | Role |
|---|---|
| `tools/revision_20261030_download_glgc.py` | fetches the three GLGC 2013 summary-statistics files (GCST002216, GCST002222, GCST002223) in parallel byte ranges |
| `tools/revision_20261030_extract_union.sh` | extracts the 134 jointly clumped variants and their FinnGen R13 associations |
| `1_mr/revision_20261030_mvmr_joint_set.R` | builds the union, joint-clumps it, fits the three-exposure MVMR on 133 variants, writes the two MVMR tables |
| `1_mr/revision_20261030_mvmr_diagnostics.R` | conditional F, Q_A and the MVMR-IVW estimates on the joint set |
| `1_mr/revision_20261030_mvmr_validate.R` | reproduces the published univariable IVW on the overlapping subset, as an allele-alignment check |
| `4_external/revision_20261030_nhanes_design.R` | surveys design (strata, PSUs, four-year and fasting weights), design-consistent intervals, the three outcome conventions |
| `4_external/revision_20261030_figure4_redraw.R` | regenerates Figure 4 from the frozen predictions |
| `4_external/revision_20261030_temporal_primary.R` | temporal validation under the primary definition, with the transported score added |

## Added aggregate tables

- `MVMR_core_TG_LDL_HDL.csv` — replaced. The mutually adjusted estimates now come from one
  jointly clumped set of 133 variants, so all three rows carry the same variant count. The
  file keeps its original schema, and `id.exposure` still identifies the exposure.
- `MVMR_joint_set_diagnostics.csv` — new. Conditional F per exposure and the Q_A statistic.
- `NHANES_outcome_conventions.csv` — new. The three conventions side by side, with the
  reference-model and lipid-only discrimination under each.

## How the joint-set MVMR was built

1. The union of the three genome-wide instrument lists contains 198 variants (55 TG,
   73 LDL-C, 84 HDL-C); the pairwise overlaps are 1, 8 and 6, and one variant
   (rs12748152) instruments all three exposures.
2. Joint clumping against the 1000 Genomes Phase 3 European panel (r2 < 0.001, 10 Mb)
   leaves 134 independent variants; 133 carry complete associations for all three
   exposures after harmonisation to a common effect allele (one palindromic variant is
   dropped).
3. Associations are taken from the full GLGC 2013 summary statistics, so every variant
   enters with an association for every exposure rather than for one exposure only.
4. Conditional F and Q_A are computed with `MVMR` 0.4.9 (`strength_mvmr`, `pleiotropy_mvmr`).

## Two things worth knowing when re-running this

One bookkeeping note: the interval printed in Table 4 and Supplementary Table S31
(0.434-0.518) comes from the primary-frame bootstrap in `table4_two_frames.csv`; the phase-20
replicate committed as `NHANES_outcome_conventions.csv` gives 0.517 at the upper limit, a seed
difference in the third decimal.

- **FinnGen reports beta for the ALT allele.** `extract_union.sh` emits `ref`, `alt`,
  `beta`, `sebeta`; `revision_20261030_mvmr_joint_set.R` swaps the two allele columns
  before alignment. Forgetting the swap flips the sign of the outcome association and
  produces a TG estimate of about 1.08 instead of 0.93. `revision_20261030_mvmr_validate.R`
  is the guard: on the 28 TG instruments shared with the archived harmonised table, the
  univariable IVW must return 0.924 (0.839-1.017) from both sources, and it does.
- **EBI throttles single connections to roughly 10 KB/s.** The downloader uses ten
  parallel byte-range requests per file; the last chunk of each file often stalls and has
  to be completed with a manual `curl -r` before the parts are concatenated and checked
  with `gzip -t`.
