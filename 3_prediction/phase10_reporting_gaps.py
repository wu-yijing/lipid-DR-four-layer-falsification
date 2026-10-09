# -*- coding: utf-8 -*-
"""phase10 (layer 2): close the reporting gaps identified in the revision audit.

What this adds, and why
-----------------------
An internal audit of the manuscript against this repository found that several
statistics the paper reports carry **no uncertainty interval**, that one
collinearity diagnostic had an unstated basis, and that two supporting tables
could not be traced to a script. None of these change a conclusion, but all of
them are reportable items under TRIPOD and are the kind of omission that inflates
the assessed risk of bias. This script computes them.

Outputs (all new files; nothing already committed is overwritten)
----------------------------------------------------------------
``layer2_reporting_gaps.csv``
    Point estimate, 95% CI and explicit basis for the quantities that previously
    had no interval: complete-case apparent AUC, optimism-corrected AUC, and the
    lipid IDI over the age + sex base model. The basis of every interval is named
    in a column so no reader has to guess which bootstrap was used.
``layer2_design_condition_numbers.csv``
    Condition number of the design matrix under three conventions (raw, column-
    standardised, and the centred/scaled correlation matrix), so the value quoted
    in the manuscript can be attributed to one of them.
``layer2_parsimonious_markers_classified.csv``
    The ten parsimonious single-marker models with an explicit ``marker_class``
    (conventional lipid / apolipoprotein / derived index), which resolves an
    ambiguous count in the Methods.
``layer2_second_grader_composition.csv``
    Manuscript-ready, de-identified composition of the blinding sample, including
    the size of the random population-representative stratum.
``layer2_pmsampsize_verification.csv``
    Verification of every sample-size requirement quoted in the manuscript against
    the committed pmsampsize grid.
``layer2_reporting_gap_inventory.csv``
    Which headline statistics still lack an interval, and what it would take to
    supply one.

Inputs: ``dev_set_v2.csv`` (resolved through ``tools/paths.py``); committed
aggregate tables for the parts that do not need individual-level data.

Usage
-----
    python 3_prediction/phase10_reporting_gaps.py
    python 3_prediction/phase10_reporting_gaps.py --B 400
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import paths  # noqa: E402

SEED = 20261009
PREDICTORS = ["age10", "male", "tg_log", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b"]
BASE = ["age10", "male"]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def auc(y: np.ndarray, s: np.ndarray) -> float:
    """Rank-based (Mann-Whitney) AUC; ties handled by average ranks."""
    y = np.asarray(y, float)
    r = pd.Series(np.asarray(s, float)).rank().values
    n1, n0 = y.sum(), len(y) - y.sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1)


def design(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    return sm.add_constant(df[cols].to_numpy(float), has_constant="add")


def fit_predict(X: np.ndarray, y: np.ndarray, idx_fit, idx_eval):
    """Return predicted probabilities on idx_eval after fitting on idx_fit."""
    try:
        m = sm.Logit(y[idx_fit], X[idx_fit]).fit(disp=0, maxiter=100)
    except Exception:
        return None
    return m.predict(X[idx_eval])


def idi(y: np.ndarray, p_full: np.ndarray, p_base: np.ndarray) -> float:
    """Integrated discrimination improvement (Pencina): change in the mean
    difference of predicted probabilities between cases and non-cases."""
    y = np.asarray(y, bool)
    return float((p_full[y].mean() - p_full[~y].mean())
                 - (p_base[y].mean() - p_base[~y].mean()))


# --------------------------------------------------------------------------
# 1-2. intervals for the layer-2 performance statistics
# --------------------------------------------------------------------------
def performance_intervals(cc: pd.DataFrame, B: int, rng) -> tuple[list[dict], dict]:
    y = cc["any_dr"].to_numpy(float)
    Xa, Xb = design(cc, BASE), design(cc, PREDICTORS)
    n = len(cc)
    idx = np.arange(n)

    pa = fit_predict(Xa, y, idx, idx)
    pb = fit_predict(Xb, y, idx, idx)
    auc_a, auc_b = auc(y, pa), auc(y, pb)
    idi_app = idi(y, pb, pa)

    # Harrell bootstrap optimism, identical in construction to the original
    # phase2 model script: optimism = mean(apparent on replicate - fitted on
    # replicate then evaluated on the original sample).
    opt_a, opt_b = [], []
    boot_a, boot_b, boot_idi = [], [], []
    for _ in range(B):
        s = rng.integers(0, n, n)
        if y[s].sum() in (0, n):
            continue
        aa = fit_predict(Xa, y, s, s)
        ab = fit_predict(Xb, y, s, s)
        ta = fit_predict(Xa, y, s, idx)
        tb = fit_predict(Xb, y, s, idx)
        if any(v is None for v in (aa, ab, ta, tb)):
            continue
        a_app, b_app = auc(y[s], aa), auc(y[s], ab)
        a_tst, b_tst = auc(y, ta), auc(y, tb)
        opt_a.append(a_app - a_tst)
        opt_b.append(b_app - b_tst)
        boot_a.append(a_app)
        boot_b.append(b_app)
        boot_idi.append(idi(y[s], ab, aa))

    opt_a, opt_b = np.array(opt_a), np.array(opt_b)
    boot_a, boot_b, boot_idi = np.array(boot_a), np.array(boot_b), np.array(boot_idi)

    def pct(v, name, point, basis, note=""):
        lo, hi = np.nanpercentile(v, [2.5, 97.5])
        return dict(metric=name, value=round(float(point), 4),
                    lo=round(float(lo), 4), hi=round(float(hi), 4),
                    basis=basis, n=n, events=int(y.sum()), note=note)

    rows = [
        pct(boot_b, "complete-case AUC, apparent (eight-predictor model)", auc_b,
            "percentile bootstrap of the apparent estimate",
            "apparent: fitted and evaluated on the same 1,417 records"),
        pct(boot_b - opt_b, "complete-case AUC, optimism-corrected (eight-predictor model)",
            auc_b - opt_b.mean(),
            "percentile bootstrap of the optimism-corrected estimate",
            "Harrell bootstrap optimism, B=%d; this is the 0.579 quoted in Section 3.2" % B),
        pct(boot_a, "complete-case AUC, apparent (age + sex base model)", auc_a,
            "percentile bootstrap of the apparent estimate",
            "reference for the incremental comparisons"),
        pct(boot_idi, "lipid IDI over the age + sex base model", idi_app,
            "percentile bootstrap of the paired IDI",
            "paired: both models refitted within each replicate"),
    ]

    extra = dict(
        apparent_auc_8pred=round(float(auc_b), 4),
        optimism_8pred=round(float(opt_b.mean()), 4),
        optimism_age_sex=round(float(opt_a.mean()), 4),
        apparent_idi=round(float(idi_app), 4),
        B_effective=int(len(boot_a)),
        B_requested=int(B),
    )
    return rows, extra


def calibration_apparent(cc: pd.DataFrame) -> dict:
    y = cc["any_dr"].to_numpy(float)
    X = design(cc, PREDICTORS)
    p = fit_predict(X, y, np.arange(len(cc)), np.arange(len(cc)))
    lp = np.log(p / (1 - p))
    # recalibration regression on the linear predictor
    m = sm.GLM(y, sm.add_constant(lp), family=sm.families.Binomial()).fit()
    return dict(calib_slope_apparent=round(float(m.params[1]), 4),
                calib_intercept_apparent=round(float(m.params[0]), 6),
                note="apparent: fitted and recalibrated on the same records, so it is not "
                     "evidence of calibration in new data")


# --------------------------------------------------------------------------
# 3. condition number under three conventions
# --------------------------------------------------------------------------
def condition_numbers(cc: pd.DataFrame) -> list[dict]:
    Xr = design(cc, PREDICTORS)                      # raw, with intercept
    Xc = np.column_stack([np.ones(len(cc)),
                          (cc[PREDICTORS] - cc[PREDICTORS].mean())
                          / cc[PREDICTORS].std(ddof=1)])   # centred + unit SD
    corr = np.corrcoef(cc[PREDICTORS].to_numpy(float).T)

    def kappa(M):
        s = np.linalg.svd(M, compute_uv=False)
        return float(s[0] / s[-1])

    return [
        dict(basis="raw (mmol/L, g/L and years on their measured scales, with intercept)",
             condition_number=round(kappa(Xr), 2)),
        dict(basis="column-standardised (centred and scaled to unit SD, with intercept)",
             condition_number=round(kappa(Xc), 2)),
        dict(basis="correlation matrix of the eight predictors (no intercept)",
             condition_number=round(kappa(corr), 2)),
    ]


# --------------------------------------------------------------------------
# 4-6. tables that do not need individual-level data
# --------------------------------------------------------------------------
def marker_classes() -> list[dict]:
    """Attach an explicit class to each parsimonious single-marker model."""
    src = paths.out_dir() / "layer2_parsimonious_models.csv"
    conv = {"TG (log2)", "Total cholesterol", "HDL-C", "LDL-C"}
    apo = {"ApoA-I", "ApoB"}
    out = []
    for r in csv.DictReader(src.open(encoding="utf-8-sig")):
        lab = r["marker"]
        cls = ("conventional lipid" if lab in conv
               else "apolipoprotein" if lab in apo else "derived index")
        out.append(dict(marker=lab, marker_class=cls, n=r["n"], events=r["events"],
                        auc_age_sex=r["auc_age_sex"], auc_plus_marker=r["auc_plus_marker"],
                        dAUC=r["dAUC"], dAUC_lo=r["dAUC_lo"], dAUC_hi=r["dAUC_hi"],
                        dAUC_optimism_corrected=r["dAUC_optimism_corrected"]))
    return out


STRATUM_EN = {
    "A1 all G4 proliferative": "every record graded proliferative DR",
    "A2 all G3 severe NPDR": "every record graded severe NPDR",
    "A3 all G1 mild NPDR": "every record graded mild NPDR (the rarest gradable stratum)",
    "A4 random G2 moderate NPDR": "random sample of moderate NPDR",
    "A5 random G0 no DR": "random sample of grade 0 (no DR)",
    "B1 \"suspected\" stated, no grade given": "records stating suspected retinopathy without a grade",
    "B2 graded G0 but positive lesion flag": "records graded no DR that carry a positive lesion descriptor",
    "C random population-representative": "random sample of the analysis population (reproduces prevalence)",
    "D1 excluded, DR mentioned, no grade": "excluded records that mention DR without a grade",
    "D2 excluded, DR not mentioned": "excluded records that do not mention DR",
    "D3 drove pre-adjudication disagreement": "the records that drove pre-adjudication disagreement",
}


def second_grader_composition() -> tuple[list[dict], dict]:
    sm = json.loads((paths.out_dir() / "layer2_second_grader_sampling_summary.json")
                    .read_text(encoding="utf-8"))
    rows = []
    for k in sorted(sm["by_stratum"]):
        rows.append(dict(stratum=k,
                         description=STRATUM_EN.get(k, "**undescribed**"),
                         n=sm["by_stratum"][k],
                         group=("A" if k.startswith("A") else "B" if k.startswith("B")
                                else "C" if k.startswith("C") else "D")))
    excluded = sum(r["n"] for r in rows if r["group"] == "D")
    random_n = sum(r["n"] for r in rows if r["group"] == "C")
    meta = dict(total=sm["total"], sum_of_strata=sum(r["n"] for r in rows),
                random_stratum=random_n, excluded_records=excluded,
                gradable=sm["gradable"], not_gradable=sm["not_gradable"],
                enriched_dr_prevalence_pct=sm["enriched_DR_prevalence_pct"], seed=sm["seed"])
    return rows, meta


def pmsampsize_verification() -> list[dict]:
    """Check every sample-size requirement quoted in the manuscript against the grid."""
    grid = {float(r["C"]): int(r["n_required"])
            for r in csv.DictReader((paths.out_dir() / "layer2_pmsampsize_grid_fine.csv")
                                    .open(encoding="utf-8-sig"))}
    n_primary = 1971
    rows = []
    for c, quoted in ((0.67, 1058), (0.65, 1369), (0.61, 2578), (0.60, 3132), (0.58, 4900)):
        got = grid.get(c)
        if c == 0.58:
            rows.append(dict(C=0.58, quoted_in_manuscript="exceeds 4,900", grid_value=got,
                             agreement="consistent" if got and got > 4900 else "MISMATCH"))
        else:
            rows.append(dict(C=c, quoted_in_manuscript=quoted, grid_value=got,
                             agreement="consistent" if got == quoted else "MISMATCH"))
    # the requirement falls as assumed discrimination rises, so the study is
    # adequately powered once n_required drops to the cohort size; the break-even
    # is the SMALLEST assumed C at which that happens
    below = [c for c, v in sorted(grid.items()) if v <= n_primary]
    breakeven = min(below) if below else None
    rows.append(dict(C="break-even", quoted_in_manuscript="C ~= 0.63",
                     grid_value=breakeven,
                     agreement="consistent" if breakeven == 0.63 else "CHECK"))
    return rows


def gap_inventory() -> list[dict]:
    """What the manuscript still reports without an interval."""
    return [
        dict(statistic="layer-2 complete-case AUC (apparent and optimism-corrected)",
             interval="supplied here", where="layer2_reporting_gaps.csv",
             action="quote the interval in Section 3.2 / Table 3"),
        dict(statistic="layer-2 lipid IDI", interval="supplied here",
             where="layer2_reporting_gaps.csv", action="quote the interval"),
        dict(statistic="layer-2 NRI at the 25.5% threshold",
             interval="not supplied", where="-",
             action="quote a bootstrap interval, or drop NRI and reason from the IDI"),
        dict(statistic="multiply-imputed pooled AUC (0.582)",
             interval="not supplied by this script", where="3_prediction/phase10b_mi_auc_ci.R",
             action="re-run the converged MI with the settings recorded there; note that the "
                    "interval that carries the paper's claim is already reported, namely the "
                    "paired bootstrap of the complete-case minus imputation difference "
                    "(+0.004, 95% CI -0.004 to 0.012, layer2_delta_auc_bootstrap.csv)"),
        dict(statistic="layer-2 apparent calibration slope 1.000 / intercept 0.000",
             interval="n/a (exact by construction)", where="layer2_reporting_gaps.csv",
             action="label these as apparent wherever they are quoted, so they are not read as "
                    "calibration in new data"),
        dict(statistic="layer-4 NRI 0.410 and IDI 0.152",
             interval="not supplied",
             where="layer4_increment_NRI_IDI.csv (SUPERSEDED; see PROVENANCE.md)",
             action="recompute on the current fold split, and supply intervals"),
        dict(statistic="collinearity: design-matrix condition number",
             interval="n/a", where="layer2_design_condition_numbers.csv",
             action="state which of the three bases the quoted value refers to"),
    ]


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--B", type=int, default=1000, help="bootstrap replicates (default 1000)")
    a = ap.parse_args(argv)
    out = paths.out_dir()
    rng = np.random.default_rng(SEED)

    print("paths:", paths.describe())
    dev = pd.read_csv(paths.require("dev_set_v2"), low_memory=False)
    dev["log2tg"] = np.log2(dev["tg"])
    g = dev[dev["gradeable"] == 1].copy()
    cc = g[["age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr"]].dropna()
    cc = cc.rename(columns={"age": "age10"})
    cc["age10"] = cc["age10"] / 10.0
    cc["tg_log"] = np.log2(cc["tg"])
    print("complete-case analytic set: n=%d, events=%d" % (len(cc), int(cc["any_dr"].sum())))

    rows, extra = performance_intervals(cc, a.B, rng)
    rows.append(dict(metric="complete-case calibration slope (apparent)",
                     value=calibration_apparent(cc)["calib_slope_apparent"],
                     lo="", hi="", basis="recalibration regression on the linear predictor",
                     n=len(cc), events=int(cc["any_dr"].sum()),
                     note="exact by construction when fitted and evaluated on the same records; "
                          "must be labelled 'apparent'"))
    rows.append(dict(metric="complete-case calibration intercept (apparent)",
                     value=calibration_apparent(cc)["calib_intercept_apparent"],
                     lo="", hi="", basis="recalibration regression on the linear predictor",
                     n=len(cc), events=int(cc["any_dr"].sum()),
                     note="must be labelled 'apparent'"))
    pd.DataFrame(rows).to_csv(out / "layer2_reporting_gaps.csv", index=False)
    print("\n" + pd.DataFrame(rows).to_string(index=False))
    print("\nbootstrap detail:", extra)

    pd.DataFrame(condition_numbers(cc)).to_csv(out / "layer2_design_condition_numbers.csv", index=False)
    pd.DataFrame(marker_classes()).to_csv(out / "layer2_parsimonious_markers_classified.csv", index=False)
    comp, meta = second_grader_composition()
    pd.DataFrame(comp).to_csv(out / "layer2_second_grader_composition.csv", index=False)
    pd.DataFrame(pmsampsize_verification()).to_csv(out / "layer2_pmsampsize_verification.csv", index=False)
    pd.DataFrame(gap_inventory()).to_csv(out / "layer2_reporting_gap_inventory.csv", index=False)

    k = condition_numbers(cc)
    for r in k:
        print("condition number [%s]: %.2f" % (r["basis"], r["condition_number"]))
    print("second-grader sample: %d records = %d strata sum, random stratum %d, "
          "excluded %d" % (meta["total"], meta["sum_of_strata"], meta["random_stratum"],
                           meta["excluded_records"]))
    bad = [r for r in pmsampsize_verification() if r["agreement"] != "consistent"]
    print("pmsampsize grid: %s" % ("all quoted values consistent" if not bad else bad))
    print("\nwrote 5 tables to %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
