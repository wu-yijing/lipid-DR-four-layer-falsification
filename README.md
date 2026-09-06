# Lipid-based DR risk stratification: a four-layer falsification

Analysis code and **desensitized aggregate results** supporting the manuscript
*"Complete-case analysis inflates discrimination and MR causal signals are
instrument-dependent: a four-layer falsification of lipid-based diabetic
retinopathy risk stratification"* (submitted to *Frontiers in Endocrinology*).

This repository contains **only analysis scripts, figure/table-generation code,
and rendered figures**. It contains **no manuscript text (`.docx`) and no
individual-level data**. Individual data are held externally under approved
ethics and are not publicly shared.

## The four-layer framework

| Layer | Question | Scripts |
|-------|----------|---------|
| 1 — Causal premise | Do lipids have an independent causal effect on DR? | `1_mr/` (two-sample MR, MVMR, MR-PRESSO, temporal, weighted) |
| 1b — cis drug-target MR | Cis-IV (druggable-genome) corroboration | `2_cis_mr/` |
| 2 — Discrimination | Can a lipid score discriminate DR in the screening cohort? | `3_prediction/` (complete-case vs multiple imputation) |
| 3 — Transportability | Does the Chinese lipid score generalise to NHANES? | `4_external/` (NHANES external validation + MI-pooled sensitivity) |
| 4 — Benchmark | What *does* discriminate DR? | `3_prediction/` (NHANES full reference model) |

## Reproducibility

- **MR / prediction / external validation:** R 4.5.2 (`TwoSampleMR`, `MendelianRandomization`, `survey`, `mice`).
- **Figures & tables:** Python 3.13 (`matplotlib`, `python-docx`, `Pillow`).
- Run order follows the layer numbering; each script reads external data from
  paths documented at the top of the file and writes aggregate outputs.

## Key quantitative findings (aggregate)

- Two-sample MR: univariable TG OR 0.86 (p = 0.009) vanishes after MVMR (0.94);
  all six cis drug-target analyses null.
- Discrimination: CC AUC 0.669 vs MI 0.556–0.559 (Δ ≈ 0.11); the gap is driven
  by DR-associated missingness (OR 2.31).
- Transportability: Chinese coefficients applied to NHANES give AUC 0.477–0.572
  (both weighted and unweighted); calibration slope 0.79 / 0.94 / −0.17 under
  MEC weighting — i.e. part of the apparent miscalibration is a survey-design
  artefact, whereas discrimination failure is robust.
- Benchmark: NHANES full model AUC 0.735 vs lipid-only 0.505.

## License

MIT — code only. Manuscript text and data are not covered by this license.
