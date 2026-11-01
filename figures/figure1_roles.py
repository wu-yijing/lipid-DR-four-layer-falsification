#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1: the four-layer design, with the role each layer plays in the argument.

Redrawn 2026-11-01. The earlier Figure 1 carried the wording of the pre-C version
("No layer answered its own question affirmatively / four eliminations"), which the
manuscript no longer uses: the claim now rests on two constitutive conditions, with the
causal layer explanatory and the benchmark bearing only on the clinical recommendation.
This script prints the role of each layer on its own panel and states the failure twice
over, so the figure and Sections 1.3, 4.1 and 5 say the same thing.

Usage
-----
    python figures/figure1_roles.py --out DIR      # writes Figure1.png and Figure1.tiff
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---- palette -------------------------------------------------------------------------
ROLE = {
    "explanatory": dict(stroke="#2C7FB8", fill="#E8F4FA", text="#1B5E8A"),
    "constitutive": dict(stroke="#B9770E", fill="#FEF9E7", text="#8A5A0B"),
    "clinical": dict(stroke="#6C3FA0", fill="#F3EDF9", text="#4E2C76"),
}
GREEN = {"stroke": "#1E8449", "fill": "#EAF7EE", "text": "#146437"}
RED = "#C0392B"
GREEN_TXT = "#1E8449"

LAYERS = [
    dict(n="Layer 1", name="Causal", role="explanatory",
         lines=[
             "Two-sample MR: six lipids (GLGC / UKB) \u2192 FinnGen R13 DR (15,353 / 88,091)",
             "IVW (main), MR-Egger, weighted median, MVMR, MR-PRESSO, cis drug-target",
             "TG 0.86 (0.77\u20130.96), p = 0.009; mutually adjusted TG 0.93 (0.83\u20131.03); all six cis targets null",
         ],
         verdict=("\u2717", RED, "No lipid met the pre-specified threshold; the TG signal is suggestive but inconclusive")),
    dict(n="Layer 2", name="Discrimination", role="constitutive",
         lines=[
             "Chinese EMR screening cohort (physician-signed grades); n = 1,971; 474 events (24.0%)",
             "Logistic + LASSO, bootstrap 1,000; multiple imputation (m = 50) as primary analysis",
             "Complete cases 0.579 vs multiple imputation 0.582; the earlier 0.06 gap was an outcome-coding artefact",
         ],
         verdict=("\u2717", RED, "Below the margin \u2014 alone sufficient to refute the claim")),
    dict(n="Layer 3", name="Transportability", role="constitutive",
         lines=[
             "Chinese coefficients applied to the NHANES 2005\u20132008 diabetes subgroup",
             "core / std / LASSO (exploratory); unweighted and MEC-weighted; ungradable photographs excluded",
             "AUC 0.475 / 0.478 / 0.554; weighted 0.457 / 0.446; calibration slopes \u22121.92 / \u22122.44",
         ],
         verdict=("\u2717", RED, "Fails to generalise; a score refitted in the target population reaches only 0.540")),
    dict(n="Layer 4", name="Benchmark", role="clinical",
         lines=[
             "NHANES full reference model: HbA1c, duration, insulin use, eGFR, BMI, smoking, blood pressure, sex",
             "Survey-weighted logistic; 10-fold out-of-fold; temporal validation across cycles",
             "AUC 0.786 (0.804 weighted) vs locally fitted lipid-only 0.540 (transported 0.467); lipids add \u22120.003",
         ],
         verdict=("\u2713", GREEN_TXT, "HbA1c, duration and insulin use robust; male sex survives BH adjustment (q = 0.01)")),
]

ROLE_LABEL = {
    "explanatory": "explanatory, not a condition of the claim",
    "constitutive": "constitutive: the claim cannot lack it",
    "clinical": "clinical: bears on the recommendation only",
}


def draw(out_dir):
    fig = plt.figure(figsize=(10.6, 10.35), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    top, gap, h = 98.0, 1.3, 19.0
    for i, L in enumerate(LAYERS):
        y1 = top - i * (h + gap)
        y0 = y1 - h
        c = ROLE[L["role"]]
        # left panel
        ax.add_patch(FancyBboxPatch((0.6, y0), 18.6, h, boxstyle="round,pad=0.2,rounding_size=0.6",
                                    fc=c["fill"], ec=c["stroke"], lw=1.1))
        ax.text(9.9, y0 + h * 0.70, L["n"], ha="center", va="center",
                fontsize=13, fontweight="bold", color="#222222")
        ax.text(9.9, y0 + h * 0.47, L["name"], ha="center", va="center",
                fontsize=11.5, color="#333333")
        ax.text(9.9, y0 + h * 0.22, L["role"], ha="center", va="center",
                fontsize=9.5, style="italic", color=c["text"])
        # right panel
        ax.add_patch(FancyBboxPatch((20.0, y0), 79.4, h, boxstyle="round,pad=0.2,rounding_size=0.6",
                                    fc="white", ec=c["stroke"], lw=1.1))
        for k, line in enumerate(L["lines"]):
            ax.text(22.0, y0 + h - 3.2 - k * 3.2, line, ha="left", va="center",
                    fontsize=10.4, color="#1A1A1A")
        sym, col, txt = L["verdict"]
        ax.text(22.0, y0 + 2.9, sym + "  " + txt, ha="left", va="center",
                fontsize=9.8, color=col, fontweight="bold")

    # ---- banner 1: what the four layers jointly establish
    y1 = 9.6
    ax.add_patch(FancyBboxPatch((4.6, y1), 94.8, 7.2, boxstyle="round,pad=0.2,rounding_size=0.4",
                                fc="#F4F4F2", ec="#555555", lw=0.9))
    ax.text(52.0, y1 + 5.15, "The claim fails twice over on its own definition",
            ha="center", va="center", fontsize=11.2, fontweight="bold", color="#111111")
    ax.text(52.0, y1 + 3.0, "a stand-alone score must discriminate on its own; a transportable score must survive export \u2014 either failure is sufficient",
            ha="center", va="center", fontsize=10.0, color="#222222")
    ax.text(52.0, y1 + 1.05, "replacement factors: glycaemic control (HbA1c), diabetes duration and treatment intensity",
            ha="center", va="center", fontsize=10.0, color="#222222")

    # ---- banner 2: the two transferable findings
    y2 = 1.6
    ax.add_patch(FancyBboxPatch((4.6, y2), 94.8, 7.2, boxstyle="round,pad=0.2,rounding_size=0.4",
                                fc="#FDF3E7", ec="#B9770E", lw=0.9))
    ax.text(52.0, y2 + 5.3, "Two generalisable sources of analytical dependence",
            ha="center", va="center", fontsize=11.2, fontweight="bold", color="#8A5A0B")
    ax.text(52.0, y2 + 3.2, "(i) an outcome-coding error mimicked a missing-data effect; correcting it removed a 0.06 gap",
            ha="center", va="center", fontsize=10.0, color="#222222")
    ax.text(52.0, y2 + 1.5, "(ii) MR causal signals are instrument-dependent (genome-wide TG 0.86 vs LPL-cis 1.27)",
            ha="center", va="center", fontsize=10.0, color="#222222")

    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, "Figure1.png")
    tif = os.path.join(out_dir, "Figure1.tiff")
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(tif, dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    print("wrote", png)
    print("wrote", tif)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__))))
    draw(ap.parse_args().out)
