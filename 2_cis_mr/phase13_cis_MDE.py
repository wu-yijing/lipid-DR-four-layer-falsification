# -*- coding: utf-8 -*-
"""Layer-1 cis drug-target arm: the smallest OR each target set could detect.

Reads the cis instrument file (iv_druggable.csv) and the fitted cis results
(druggable_mr_results.csv) produced by phase5a/phase5c, and reports, per target,
the instrument count and the minimum detectable OR at 80% and 90% power
(alpha = 0.05, two-sided) implied by the fitted standard error.
"""
from __future__ import annotations
import csv, math, os

IV = os.environ.get("LDF_CIS_IV", r"E:/workbuddy/2026-09-04-18-19-16/phase1_out/iv_druggable.csv")
RES = os.environ.get("LDF_CIS_RES", r"E:/workbuddy/2026-09-04-18-19-16/phase1_out/druggable_mr_results.csv")
OUT = os.environ.get("LDF_OUT", r"E:/workbuddy/lipid-DR-four-layer-falsification/output_tables")
Z_A, Z80, Z90 = 1.959963984540054, 0.8416212335729143, 1.2815515655446004

ivs = {}
with open(IV, newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        ivs.setdefault(r["target"], []).append(r)

rows = []
with open(RES, newline="", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        se = float(r["se"])
        fs = [(float(x["beta"]) / float(x["se"])) ** 2 for x in ivs.get(r["target"], [])]
        rows.append({
            "target": r["target"], "lipid_trait": r["lipid_trait"],
            "n_iv_extracted": r["n_iv"], "n_iv_harmonised": r["n_harmonised"], "method": r["method"],
            "mean_F": round(sum(fs) / len(fs), 1) if fs else "",
            "min_F": round(min(fs), 1) if fs else "",
            "OR": round(float(r["OR"]), 4),
            "lci": round(float(r["or_lci95"]), 3), "uci": round(float(r["or_uci95"]), 3),
            "p": round(float(r["p"]), 4),
            "MDE_OR_80pct_power": round(math.exp((Z_A + Z80) * se), 3),
            "MDE_OR_90pct_power": round(math.exp((Z_A + Z90) * se), 3),
            "excludes_OR_above": round(math.exp((Z_A + Z80) * se), 3),
        })
cols = ["target","lipid_trait","n_iv_extracted","n_iv_harmonised","method","mean_F","min_F","OR","lci","uci","p",
        "MDE_OR_80pct_power","MDE_OR_90pct_power","excludes_OR_above"]
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "layer1_cis_target_MDE.csv"), "w", newline="", encoding="utf-8") as fh:
    wr = csv.DictWriter(fh, fieldnames=cols); wr.writeheader(); wr.writerows(rows)
for r in rows:
    print(r["target"], "| n_iv=", r["n_iv_harmonised"], "| OR", r["OR"], f'({r["lci"]}-{r["uci"]})',
          "| MDE80", r["MDE_OR_80pct_power"], "| MDE90", r["MDE_OR_90pct_power"])
