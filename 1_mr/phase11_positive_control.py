# -*- coding: utf-8 -*-
"""phase11 (layer 1): the two validation analyses the manuscript still owes.

Why these two, and why they are now a precondition rather than a nice-to-have
-----------------------------------------------------------------------------
Several publishers have tightened their screening for submissions built on public
data. Frontiers, since July 2024, requires independent external validation or
additional institutional data for every health-dataset Mendelian randomisation
submission, and since 2025 it applies the same requirement to papers built on simple
queries of public data, including NHANES; PLOS One states that robustness requires
validation of findings "in, at minimum, a second independent dataset"; and the
European Journal of Epidemiology's editorial on the growth of MR literature singles
out, as hallmarks of low-quality MR work, an unjustified exposure-outcome pairing, no
discussion of gene-environment equivalence, and no sensitivity analyses.

This layer is currently the weakest part of that case. Layer 1 rests entirely on
internal sensitivity analyses: multivariable MR, MR-PRESSO, six cis drug-target
analyses and the EXMORE control-exclusion endpoint. EXMORE is a *different endpoint
definition within the same FinnGen release*, not a second independent dataset. So two
analyses are missing:

1. **An outcome-side positive control** (this script, ``--mode positive``).
   Genetically predicted type 2 diabetes, or HbA1c, against the same DR endpoint. If
   the instrument set and dataset can detect an effect that certainly exists, then the
   null for lipids is informative rather than an artefact of a pipeline with no power.
   Without it, "no effect" and "no detectable effect" cannot be separated.

2. **Replication in a second, independent DR outcome dataset** (``--mode replicate``).
   The same lipid instruments against a different DR GWAS. This is the requirement the
   editorial policies actually state.

Both modes are data-gated. Neither can be run from this repository alone, because
neither the positive-control instrument panel nor a second DR outcome dataset is
distributed here. Run ``--check`` first: it reports exactly which inputs are present,
which are missing, and how each missing one is obtained. The script refuses to emit
numbers it cannot compute rather than writing an empty or stub table.

Input contract
--------------
``<SRC>/phase1_out/positive_control_instruments.csv``  (positive mode)
    exposure,SNP,ea,nea,b_exp,se_exp,p_exp,eaf
    One row per instrument-exposure pair. Only the instrument panel is needed: the
    confounder set is handled by MR-Egger and the global test, as in layer 1.

``<SRC>/phase1_out/harmonised_data.csv``  (replicate mode)
    Already produced by phase 1: the lipid instruments with exposure and outcome
    effects, one row per exposure-SNP, column ``exposure`` giving the trait.

outcome summary statistics (either mode)
    Either a FinnGen-format ``.gz``, given by ``LDF_OUTCOME_GZ``, with columns
    rsids, ref, alt, beta, sebeta, pval; or a pre-extracted ``.tsv`` given by
    ``LDF_OUTCOME_SUBSET`` with columns rsid, ref, alt, beta, sebeta.
    For the positive control the sensible outcome is the primary DR endpoint used in
    layer 1 (so that the control speaks to the same analysis); for replication it must
    be an independent DR dataset.

Output
------
``output_tables/MR_positive_control.csv``      (positive mode)
``output_tables/MR_outcome_replication.csv``   (replicate mode)
One row per exposure with n instruments, IVW OR and CI, p-value, the CI-implied
p-value and their agreement, MR-Egger slope and intercept p, weighted median, and
Cochran's Q with I-squared.

Usage
-----
    python 1_mr/phase11_positive_control.py --check
    LDF_OUTCOME_GZ=/path/to/finngen_R13_DM_RETINOPATHY.gz \
        python 1_mr/phase11_positive_control.py --mode positive
    LDF_OUTCOME_SUBSET=/path/to/independent_DR_subset.tsv \
        python 1_mr/phase11_positive_control.py --mode replicate
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import paths  # noqa: E402
from mr_estimators import (align, ci_from_est, ci_p_consistent, egger,  # noqa: E402
                           heterogeneity, ivw, p_from_z, weighted_median)

POS_INSTR = "positive_control_instruments.csv"
REPL_INSTR = "harmonised_data.csv"
NEEDED_GZ_COLS = {"rsids", "ref", "alt", "beta", "sebeta"}
NEEDED_TSV_COLS = {"rsid", "ref", "alt", "beta", "sebeta"}
NEEDED_INSTR_COLS = {"exposure", "SNP", "ea", "nea", "b_exp", "se_exp"}

ACQUISITION = {
    "positive_control_instruments.csv": (
        "Build the instrument panel exactly as in layer 1: genome-wide significant "
        "(p < 5e-8), LD-independent (r2 < 0.001, 10 Mb) variants for the positive-control "
        "exposure. For type 2 diabetes, a suitably powered European-ancestry GWAS "
        "(e.g. DIAMANTE / DIAGRAM) or, for HbA1c, a glycaemic-traits GWAS. Export the "
        "columns listed in the input contract to <SRC>/phase1_out/."
    ),
    "harmonised_data.csv": (
        "Already present: produced by the phase-1 scripts and used by phase 2. If it is "
        "absent, re-run the phase-1 layer."
    ),
    "outcome_gz": (
        "FinnGen release R13, endpoint DM_RETINOPATHY (the layer-1 primary endpoint) for "
        "the positive control; download from https://www.finngen.fi/en/access_results and "
        "set LDF_OUTCOME_GZ. Note that this file is ~800 MB and the script streams it once."
    ),
    "outcome_subset": (
        "For replication, any independent diabetic-retinopathy GWAS. The GWAS Catalog "
        "summary statistics are fully public and need no token: "
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ . Pre-extract the "
        "instrument SNPs to a tsv with columns rsid, ref, alt, beta, sebeta and set "
        "LDF_OUTCOME_SUBSET. A multi-ancestry diabetic eye disease GWAS is one candidate."
    ),
}


# --------------------------------------------------------------------------
#: column aliases accepted for the instrument panel. The layer-1 harmonised data is a
#: TwoSampleMR frame, the drug-target panel uses short names, and either is fine.
ALIASES = {
    "ea": ["ea", "effect_allele.exposure", "effect_allele"],
    "nea": ["nea", "other_allele.exposure", "other_allele"],
    "b_exp": ["b_exp", "beta.exposure", "beta"],
    "se_exp": ["se_exp", "se.exposure", "se"],
    "p_exp": ["p_exp", "pval.exposure", "pval"],
    "eaf": ["eaf", "eaf.exposure", "eaf"],
}


def normalise_instruments(d: pd.DataFrame) -> pd.DataFrame:
    """Map either panel schema onto the internal contract, and apply mr_keep if present."""
    out = d.copy()
    for target, names in ALIASES.items():
        if target in out.columns:
            continue
        for n in names:
            if n in out.columns:
                out[target] = out[n]
                break
    for keep in ("mr_keep", "mr_keep.outcome"):
        if keep in out.columns:
            vals = out[keep]
            if vals.dtype == bool:
                out = out[vals]
            else:
                out = out[vals.astype(str).str.lower().isin(["true", "1", "yes"])]
            break
    return out


def instrument_path(which: str) -> Path:
    return paths.src_dir() / "phase1_out" / (POS_INSTR if which == "positive" else REPL_INSTR)


def outcome_source() -> tuple[Path | None, str]:
    gz = os.environ.get("LDF_OUTCOME_GZ")
    tsv = os.environ.get("LDF_OUTCOME_SUBSET")
    if tsv:
        return Path(tsv).expanduser(), "subset"
    if gz:
        return Path(gz).expanduser(), "gz"
    return None, "none"


def check(mode: str | None = None) -> int:
    """Report input availability. `mode` limits the requirement to the inputs that
    mode actually needs; None requires everything (the bare --check)."""
    print("phase 11 input check")
    print("=" * 72)
    print(paths.describe())
    ok = True

    need_pos = mode in (None, "positive")
    need_rep = mode in (None, "replicate")

    ip = instrument_path("positive")
    print("\n[positive control]%s" % ("" if need_pos else "  (not required for this mode)"))
    print("  instrument panel : %s  %s" % (ip, "PRESENT" if ip.exists() else "MISSING"))
    if ip.exists():
        cols = set(next(csv.reader(ip.open(encoding="utf-8-sig"))))
        miss = NEEDED_INSTR_COLS - cols
        print("  columns          : %s" % ("ok" if not miss else "MISSING %s" % sorted(miss)))
        ok &= not miss
    elif need_pos:
        ok = False
        print("  how to obtain    : %s" % ACQUISITION[POS_INSTR])

    rp = instrument_path("replicate")
    print("\n[outcome replication]%s" % ("" if need_rep else "  (not required for this mode)"))
    print("  lipid instruments: %s  %s" % (rp, "PRESENT" if rp.exists() else "MISSING"))
    if not rp.exists() and need_rep:
        ok = False
        print("  how to obtain    : %s" % ACQUISITION[REPL_INSTR])

    src, kind = outcome_source()
    print("\n[outcome summary statistics]")
    if src is None:
        ok = False
        print("  not configured   : set LDF_OUTCOME_GZ or LDF_OUTCOME_SUBSET")
        print("  positive control : %s" % ACQUISITION["outcome_gz"])
        print("  replication      : %s" % ACQUISITION["outcome_subset"])
    else:
        print("  source (%s)  : %s  %s" % (kind, src, "PRESENT" if src.exists() else "MISSING"))
        if not src.exists():
            ok = False

    print("\n" + "=" * 72)
    if ok:
        print("all inputs present: the requested mode can run")
        return 0
    print("inputs missing: phase 11 cannot be run here.")
    print("This is expected in the published repository, which carries no GWAS extracts;")
    print("it is reported rather than worked around so that no stub result is produced.")
    return 3


# --------------------------------------------------------------------------
def load_outcome(snps: set[str]) -> pd.DataFrame:
    src, kind = outcome_source()
    if src is None or not src.exists():
        raise paths.MissingInput(
            "no outcome summary statistics configured; set LDF_OUTCOME_GZ or "
            "LDF_OUTCOME_SUBSET (run --check for the contract)")
    if kind == "subset":
        d = pd.read_csv(src, sep="\t")
        d = d.rename(columns={"rsid": "SNP"})
        need = {"SNP", "ref", "alt", "beta", "sebeta"}
        if not need <= set(d.columns):
            raise ValueError("subset tsv needs columns %s, found %s"
                             % (sorted(need), sorted(d.columns)))
        return d[["SNP", "ref", "alt", "beta", "sebeta"]].astype({"beta": float, "sebeta": float})

    rows = []
    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        idx = {c: k for k, c in enumerate(hdr)}
        missing = NEEDED_GZ_COLS - set(hdr)
        if missing:
            raise ValueError("outcome .gz needs columns %s, found %s"
                             % (sorted(NEEDED_GZ_COLS), sorted(hdr)))
        remaining = set(snps)
        for r in rd:
            hit = [s for s in r[idx["rsids"]].split(",") if s in remaining]
            if hit:
                rows.append((hit[0], r[idx["ref"]], r[idx["alt"]],
                             r[idx["beta"]], r[idx["sebeta"]]))
                remaining.discard(hit[0])
            if not remaining:
                break
    d = pd.DataFrame(rows, columns=["SNP", "ref", "alt", "beta", "sebeta"])
    return d.astype({"beta": float, "sebeta": float})


def run(instr: pd.DataFrame, out: pd.DataFrame, label: str) -> list[dict]:
    """One IVW / MR-Egger / weighted-median block per exposure."""
    rows = []
    for exposure, g in instr.groupby("exposure"):
        m = g.merge(out, left_on="SNP", right_on="SNP", how="inner")
        b, se, bx = [], [], []
        for _, r in m.iterrows():
            sg = align(r["ea"], r["nea"], r["ref"], r["alt"])
            if sg is None or float(r["b_exp"]) == 0:
                continue
            b.append(sg * float(r["beta"]) / float(r["b_exp"]))
            se.append(float(r["sebeta"]) / abs(float(r["b_exp"])))
            bx.append(float(r["b_exp"]))
        b, se, bx = np.array(b), np.array(se), np.array(bx)
        if len(b) == 0:
            rows.append(dict(analysis=label, exposure=exposure, n=0))
            continue
        if len(b) == 1:
            est, z = float(b[0]), float(b[0]) / float(se[0])
            se_est = float(se[0])
            Q = df_ = i2 = p_het = np.nan
            egg_est = egg_p = ic = ic_p = wm = np.nan
        else:
            est, se_est, Q, df_ = ivw(b, se)
            z = est / se_est
            _, _, p_het, i2 = heterogeneity(b, se)
            egg_est, egg_se, ic, ic_se = egger(b, se, bx)
            egg_p = p_from_z(egg_est / egg_se)
            ic_p = p_from_z(ic / ic_se)
            wm = weighted_median(b, se)
        lo, hi = ci_from_est(est, se_est)
        p = p_from_z(z)
        p_ci, consistent = ci_p_consistent(lo, hi, p)
        rows.append(dict(
            analysis=label, exposure=exposure, n=len(b),
            OR=np.exp(est), lo=lo, hi=hi, p=p, p_from_ci=p_ci, ci_p_consistent=consistent,
            egger_OR=np.exp(egg_est) if np.isfinite(egg_est) else np.nan, egger_slope_p=egg_p,
            egger_intercept_p=ic_p,
            weighted_median_OR=np.exp(wm) if np.isfinite(wm) else np.nan,
            cochran_Q=Q, Q_df=df_, Q_p=p_het, I2=i2))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="report input availability and exit (0 = runnable, 3 = not runnable here)")
    ap.add_argument("--mode", choices=["positive", "replicate"], default="positive")
    a = ap.parse_args(argv)

    if a.check:
        return check(a.mode if a.mode else None)

    if check(a.mode) != 0:
        print("\nRefusing to run: no numbers are emitted for an analysis that cannot be computed.")
        return 3

    ip = instrument_path(a.mode)
    instr = normalise_instruments(pd.read_csv(ip))
    missing = NEEDED_INSTR_COLS - set(instr.columns)
    if missing:
        print("instrument panel is missing columns %s (see --check)" % sorted(missing))
        return 3

    snps = set(instr["SNP"].astype(str))
    print("instruments: %d rows, %d unique SNPs" % (len(instr), len(snps)))
    out = load_outcome(snps)
    print("outcome rows extracted: %d of %d requested" % (len(out), len(snps)))

    rows = run(instr, out, "positive_control" if a.mode == "positive" else "outcome_replication")
    res = pd.DataFrame(rows)
    name = "MR_positive_control.csv" if a.mode == "positive" else "MR_outcome_replication.csv"
    res.to_csv(paths.out_dir() / name, index=False)
    print(res.to_string(index=False, float_format=lambda x: "%.4g" % x))
    print("\nwrote %s" % (paths.out_dir() / name))
    print("\nReminder: a positive control is reported only if it is positive. If this table "
          "shows a null for an exposure with a known effect, the layer-1 null for lipids "
          "cannot be read as informative, and that must be stated in the manuscript.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
