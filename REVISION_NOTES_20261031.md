# Revision round 1 — addendum, 31 October 2026

Two items closed after the 30 October revision package.

## Linkage-window sensitivity (reviewer 1, point 12)

The manuscript states that restricting the layer-2 analysis to the 1,357 records whose
lipid panel fell within 7 days of the photograph left discrimination unchanged
(apparent 0.594, optimism-corrected 0.573, 350 events). That analysis had been run
ad hoc and no script was committed with the 30 October package. It is now archived as
`3_prediction/phase22_linkage_window.R`, with its output in
`output_tables/layer2_linkage_window.csv`:

| Analysis | n | Events | Apparent | Optimism | Corrected |
|---|---|---|---|---|---|
| All complete cases | 1,417 | 367 | 0.5984 | 0.0191 | 0.5793 |
| Panel within 7 days | 1,357 | 350 | 0.5939 | 0.0212 | **0.5727** |
| Panel outside 7 days | 60 | 17 | 0.7934 | 0.1165 | 0.6769 |

The first two rows reproduce the printed values (0.598 / 0.579 and 0.594 / 0.573) to
four decimals, and the counts match exactly.

The third row was not previously reported and is worth knowing: the 60 records whose
panel falls outside the window behave very differently from the rest (apparent 0.79 on
17 events). That is the expected consequence of a linkage that carries no shared
identifier — those records may belong to another admission or another patient — and it
supports rather than weakens the decision to report the window-restricted estimate. It
is not quoted in the manuscript because 17 events cannot support an estimate.

## How the window was read

`cleaned_dataset.xlsx` (sheet 3, `血脂-糖网关联队列`) carries `关联天数差` in column 49.
The modelling frame `dev_set_v2.csv` carries the corrected outcome but not the interval,
so the two are joined by row order; the row alignment is checked by `stopifnot(nrow(d) == nrow(x))`
plus an exact age match on 2,481 of 2,482 rows (one record has a missing age).

## Note for anyone re-running this

`readxl` hands the workbook path to `utils::unzip`, which cannot translate a UTF-8
filename under a non-UTF-8 code page. The script therefore writes the path as `\u`
escapes and stages a copy into `tempdir()` before reading. Without that step the script
fails with "path does not exist" even though `file.exists()` returns TRUE.

## Supplementary Data File 1 (31 October 2026)

The manuscript declares, in the supplementary-material paragraph, a CSV of per-instrument
strength, allele frequency and effect size for all 662 instrument–exposure pairs. That file
is `output_tables/layer1_instrument_detail_F_R2.csv`. The submission copy,
`Supplementary_Data_File_1.csv`, is identical except that rows are ordered by exposure
(TG, LDL-C, HDL-C, TC, ApoA-I, ApoB) and the exposure labels are harmonised with Table 2
(`LDL`→`LDL-C`, `HDL`→`HDL-C`, `ApoA1`→`ApoA-I`). Column names are unchanged.

Recomputed from that file alone, so that the file can be checked against the tables:

| Exposure | n | Mean F | Cumulative R² | IVW OR (95% CI) | MDE, 80% power |
|---|---|---|---|---|---|
| TG | 55 | 158.85 | 0.0518 | 0.86 (0.77–0.96) | 1.120 |
| LDL-C | 73 | 182.29 | 0.1159 | 0.99 (0.91–1.06) | 1.075 |
| HDL-C | 84 | 152.79 | 0.0793 | 1.03 (0.92–1.15) | 1.092 |
| TC | 86 | 139.67 | 0.0865 | 0.96 (0.86–1.06) | 1.090 |
| ApoA-I | 213 | 154.81 | 0.0065 † | 1.18 (0.86–1.61) | 1.371 |
| ApoB | 151 | 314.69 | 0.0074 † | 0.87 (0.66–1.15) | 1.360 |

Every value matches Supplementary Tables S22 and S23 and Table 2, including the two ranges
quoted in Methods 2.3 (minimum detectable OR 1.08–1.37 on the fixed-effect SE and 1.12–1.56
on the multiplicative random-effects SE), and 582 unique SNPs. Two footnotes were added to
S23 because the file exposes them:

- † The apolipoprotein R² is tabulated on the SD scale. Their betas are per g/L, so
  Σ 2p(1−p)b² lands in g/L² (0.0065, 0.0074); the tabulated 0.105 and 0.139 are ΣF/N.
  Both columns exist in `layer1_instrument_strength_by_exposure.csv`, and S23 printed one
  for the conventional lipids and the other for the apolipoproteins without saying so.
- Median F in S23 is the upper of the two central values when n is even (84 for HDL-C,
  86 for TC), which is 60.36 and 65.79; `median()` gives 59.81 and 65.51. The four
  odd-n exposures are unaffected. The convention is now stated in the table note.

