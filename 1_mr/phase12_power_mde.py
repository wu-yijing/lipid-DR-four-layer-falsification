# -*- coding: utf-8 -*-
"""Layer-1 instrument strength, variance explained and minimum detectable effect (MDE).

Computes, from the harmonised two-sample input already used for the main analysis:
  1. per-exposure instrument count, per-SNP F, mean/median/min F
  2. cumulative variance explained (R2 = sum 2p(1-p)b^2) for the SD-scaled GLGC traits
  3. IVW fixed- and multiplicative-random-effects estimates (cross-check vs mr_main.csv)
  4. the smallest OR the design can detect at 80% / 90% power, alpha = 0.05 two-sided
  5. achieved power at OR = 1.10, 1.15, 1.20
  6. a winner's-curse sensitivity: IVW restricted to stricter instrument p thresholds

Pure standard library; reads the harmonised file produced by phase1_mr.R.
"""
from __future__ import annotations
import csv, json, math, os, sys

SRC = os.environ.get("LDF_HARM", r"E:/workbuddy/2026-09-04-18-19-16/phase1_out/harmonised_data.csv")
OUT = os.environ.get("LDF_OUT", r"E:/workbuddy/lipid-DR-four-layer-falsification/output_tables")
Z_A = 1.959963984540054     # two-sided alpha = 0.05
Z_P80 = 0.8416212335729143
Z_P90 = 1.2815515655446004
SD_SCALED = {"HDL", "LDL", "TC", "TG"}   # GLGC exposures are per 1-SD; UKB ApoA1/ApoB are per allele

def norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))

def load(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if str(r.get("mr_keep", "")).upper() != "TRUE":
                continue
            try:
                rows.append({
                    "exposure": r["exposure"], "SNP": r["SNP"],
                    "bx": float(r["beta.exposure"]), "sx": float(r["se.exposure"]),
                    "by": float(r["beta.outcome"]),  "sy": float(r["se.outcome"]),
                    "eaf": float(r["eaf.exposure"]) if r.get("eaf.exposure") not in (None, "", "NA") else None,
                    "px": float(r["pval.exposure"]) if r.get("pval.exposure") not in (None, "", "NA") else None,
                    "ne": float(r["samplesize.exposure"]) if r.get("samplesize.exposure") not in (None, "", "NA") else None,
                    "cases": float(r["ncase.outcome"]), "ctrl": float(r["ncontrol.outcome"]),
                })
            except (TypeError, ValueError):
                continue
    return rows

def ivw(rs):
    num = sum(r["bx"] * r["by"] / r["sy"] ** 2 for r in rs)
    den = sum(r["bx"] ** 2 / r["sy"] ** 2 for r in rs)
    b = num / den
    se = 1.0 / math.sqrt(den)
    k = len(rs)
    q = sum((r["bx"] ** 2 / r["sy"] ** 2) * ((r["by"] / r["bx"]) - b) ** 2 for r in rs)
    se_re = se * math.sqrt(q / (k - 1)) if k > 1 and q > (k - 1) else se
    return {"k": k, "b": b, "se_fe": se, "se_re": se_re, "Q": q, "Q_df": k - 1}

def mde(se, z_power):
    return math.exp((Z_A + z_power) * se)

def power_at(orv, se):
    z = abs(math.log(orv)) / se
    return norm_cdf(z - Z_A) + norm_cdf(-z - Z_A)

def main():
    rows = load(SRC)
    exposures = []
    seen = []
    for r in rows:
        if r["exposure"] not in seen:
            seen.append(r["exposure"])
    strength, mde_rows, wc_rows, detail = [], [], [], []
    for exp in seen:
        rs = [r for r in rows if r["exposure"] == exp]
        fs = [(r["bx"] / r["sx"]) ** 2 for r in rs]
        fs_sorted = sorted(fs)
        r2 = sum(2 * r["eaf"] * (1 - r["eaf"]) * r["bx"] ** 2 for r in rs if r["eaf"] is not None)
        has_eaf = all(r["eaf"] is not None for r in rs)
        est = ivw(rs)
        ne = rs[0]["ne"]
        strength.append({
            "exposure": exp, "n_snp": est["k"], "n_exposure": int(ne) if ne else "",
            "mean_F": round(sum(fs) / len(fs), 3), "median_F": round(fs_sorted[len(fs) // 2], 3),
            "min_F": round(min(fs), 3), "max_F": round(max(fs), 3),
            "prop_F_gt10": round(sum(1 for f in fs if f > 10) / len(fs), 4),
            "cumulative_R2": round(r2, 5) if has_eaf else "",
            "cumulative_R2_basis": "sum 2p(1-p)b^2 (SD-scaled)" if exp in SD_SCALED else "per-allele; needs trait SD",
            "cumulative_R2_F_over_N": round(sum(fs) / ne, 5) if ne else "",
        })
        for r in rs:
            detail.append({"exposure": exp, "SNP": r["SNP"], "beta_exposure": r["bx"], "se_exposure": r["sx"],
                           "F": round((r["bx"] / r["sx"]) ** 2, 3), "eaf": r["eaf"], "pval_exposure": r["px"],
                           "beta_outcome": r["by"], "se_outcome": r["sy"]})
        mde_rows.append({
            "exposure": exp, "n_snp": est["k"],
            "OR_point": round(math.exp(est["b"]), 4),
            "se_IVW_fixed": round(est["se_fe"], 5), "se_IVW_random": round(est["se_re"], 5),
            "MDE_OR_80pct_power": round(mde(est["se_fe"], Z_P80), 3),
            "MDE_OR_90pct_power": round(mde(est["se_fe"], Z_P90), 3),
            "MDE_OR_80pct_power_RE": round(mde(est["se_re"], Z_P80), 3),
            "MDE_OR_90pct_power_RE": round(mde(est["se_re"], Z_P90), 3),
            "power_at_OR_1.10": round(power_at(1.10, est["se_fe"]), 3),
            "power_at_OR_1.15": round(power_at(1.15, est["se_fe"]), 3),
            "power_at_OR_1.20": round(power_at(1.20, est["se_fe"]), 3),
            "power_at_OR_1.30": round(power_at(1.30, est["se_fe"]), 3),
            "power_at_OR_1.20_RE": round(power_at(1.20, est["se_re"]), 3),
        })
        for thr in (5e-8, 1e-10, 5e-10, 1e-12):
            sub = [r for r in rs if r["px"] is not None and r["px"] < thr]
            if len(sub) < 2:
                wc_rows.append({"exposure": exp, "pval_threshold": thr, "n_snp": len(sub), "OR": "", "lci": "", "uci": "", "p": ""})
                continue
            e2 = ivw(sub)
            wc_rows.append({"exposure": exp, "pval_threshold": thr, "n_snp": e2["k"],
                            "OR": round(math.exp(e2["b"]), 4),
                            "lci": round(math.exp(e2["b"] - Z_A * e2["se_re"]), 4),
                            "uci": round(math.exp(e2["b"] + Z_A * e2["se_re"]), 4),
                            "p": round(2 * (1 - norm_cdf(abs(e2["b"] / e2["se_re"]))), 5)})
    os.makedirs(OUT, exist_ok=True)
    def w(name, data, cols):
        with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=cols); wr.writeheader(); wr.writerows(data)
        print("wrote", name, len(data), "rows")
    w("layer1_instrument_strength_by_exposure.csv", strength,
      ["exposure","n_snp","n_exposure","mean_F","median_F","min_F","max_F","prop_F_gt10","cumulative_R2","cumulative_R2_basis","cumulative_R2_F_over_N"])
    w("layer1_MDE_by_exposure.csv", mde_rows,
      ["exposure","n_snp","OR_point","se_IVW_fixed","se_IVW_random","MDE_OR_80pct_power","MDE_OR_90pct_power","MDE_OR_80pct_power_RE","MDE_OR_90pct_power_RE",
       "power_at_OR_1.10","power_at_OR_1.15","power_at_OR_1.20","power_at_OR_1.30","power_at_OR_1.20_RE"])
    w("layer1_winners_curse_threshold_sensitivity.csv", wc_rows,
      ["exposure","pval_threshold","n_snp","OR","lci","uci","p"])
    w("layer1_instrument_detail_F_R2.csv", detail,
      ["exposure","SNP","beta_exposure","se_exposure","F","eaf","pval_exposure","beta_outcome","se_outcome"])
    print(json.dumps(mde_rows, indent=1, ensure_ascii=False))

if __name__ == "__main__":
    main()
