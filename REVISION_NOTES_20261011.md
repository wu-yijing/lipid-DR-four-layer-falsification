# Revision notes — 2026-10-11

## Table 4 AUC confidence intervals: unpaired bootstrap corrected

**Status**: the frozen analysis pipeline is **unchanged** (freeze at `647c78b` stands). This note
records one additive script and the manuscript edit it supports.

### The defect

`4_external/phase4_external_v2.R` (line 28), `4_external/phase4_external_v2_dedup.R` (line 31) and
`1_mr/phase4_external.R` (line 25) compute the Table 4 AUC intervals as:

```r
auc_ci <- function(p, y, B = 1000) {
  ok <- !is.na(p) & !is.na(y); p <- p[ok]; y <- y[ok]; n <- length(y)
  out <- replicate(B, auc_fun(p[sample(n, n, TRUE)], y[sample(n, n, TRUE)]))
  quantile(out, c(0.025, 0.975), na.rm = TRUE)
}
```

The two `sample()` calls are **independent**, so the predictor–outcome pairing is destroyed in
every replicate. `auc_fun` then scores a randomly mismatched predictor against the outcome, and the
resulting distribution is centred on AUC = 0.5 whatever the estimate is. The returned interval is
therefore a **null interval**, not a confidence interval.

### Evidence

Midpoints of the six reported intervals (three models x two pipelines), against their point estimates:

| file | model | point estimate | reported interval | midpoint | estimate - midpoint |
|---|---|---|---|---|---|
| `layer3_external_validation_AUC.csv` | core | 0.4979 | 0.4616-0.5379 | **0.4998** | -0.0019 |
| " | std | 0.4928 | 0.4443-0.5573 | **0.5008** | -0.0080 |
| " | lasso | 0.5510 | 0.4152-0.5829 | **0.4990** | **+0.0519** |
| `layer3_dedup_external_validation_AUC.csv` | core | 0.4904 | 0.4615-0.5380 | **0.4998** | -0.0094 |
| " | std | 0.4848 | 0.4453-0.5583 | **0.5018** | -0.0170 |
| " | lasso | 0.5508 | 0.4152-0.5845 | **0.4999** | +0.0510 |
| `layer3_recomputed_AUC_CI_B2000.csv` (**paired**) | core | 0.4979 | 0.4592-0.5368 | **0.4980** | -0.00006 |
| " | std | 0.4928 | 0.4352-0.5497 | **0.4924** | +0.00040 |

All six defective midpoints sit in 0.4990-0.5018 while the point estimates span 0.485-0.551; the
correctly paired intervals agree with their own point estimates to 1e-4. The LASSO row, the only
model whose estimate departs appreciably from 0.5, displayed the defect as an apparently
asymmetric interval (0.415-0.583 around 0.551).

The correct paired implementation was already in use elsewhere in the repository
(`phase4e_local*.R: sample(n, n, TRUE)` applied as a single index to both `p` and `y`;
likewise `phase10_local_comparator_increment.R` and `phase11_lipid_increment_over_full_model.R`).

### The correction

`4_external/phase21_paired_auc_ci.R` (new, additive):

1. rebuilds the three transported predictions with the identical code path;
2. **GATE 1** asserts the rebuilt n / events / AUC reproduce `phase4_out_v2/phase4_external_auc.csv`
   to 1e-12 - passed (max |dAUC| = 3.3e-16);
3. recomputes percentile + BCa intervals with the pairing fixed, using the same seed (20260912),
   B (2000) and operation order as `phase4e_local_v2.R`;
4. **GATE 2** asserts core and std reproduce the archived `phase4_out_v2/phase4e_auc_ci.csv` -
   passed (max |d| = 2.8e-16);
5. writes `output_tables/layer3_paired_AUC_CI.csv`.

The LASSO row had never been recomputed; it is computed here for the first time.

### Result

| model | n | events | AUC | **was (unpaired)** | **now (paired, percentile B=2000)** | BCa |
|---|---|---|---|---|---|---|
| core | 1,041 | 302 | 0.498 | 0.462-0.538 | **0.459-0.537** | 0.460-0.537 |
| std | 494 | 148 | 0.493 | 0.444-0.557 | **0.435-0.550** | 0.438-0.548 |
| LASSO (exploratory) | 190 | 64 | 0.551 | 0.415-0.583 | **0.466-0.638** | 0.458-0.633 |

No point estimate changes; no conclusion changes. The corrected core and std intervals now agree
with Supplementary Table S31a, which already used the paired values - the two were previously
inconsistent (Table 4 versus S31a).

### Manuscript edits made to match

- Table 4, column "AUC unweighted (95% CI)": the three intervals replaced.
- Table 4 footnote: one sentence added stating that the unweighted AUC intervals are
  2,000-replicate percentile bootstrap intervals over records, resampling predictor and outcome
  together.
- Analysis freeze declaration: `4_external/phase21_paired_auc_ci.R` named as the producer of the
  Table 4 intervals.

### Open

The `auc_ci()` functions in the three frozen scripts are left as they were, so that the freeze at
`647c78b` remains a faithful description of what was run. Anyone re-running those scripts will
reproduce the defective intervals; `phase21_paired_auc_ci.R` is the authoritative producer of the
intervals published in Table 4.
