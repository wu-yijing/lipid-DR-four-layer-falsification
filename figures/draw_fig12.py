# -*- coding: utf-8 -*-
"""Figures 1-2 rebuild: band-layout flowchart + fixed forest layout."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import csv, os

plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.8,
                     "savefig.dpi": 300, "figure.dpi": 150})
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "verm": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9", "grey": "#7F7F7F",
      "red": "#C1272D"}
MM = 1 / 25.4
OUT = r"E:/workbuddy/2026-09-04-18-19-16/figs_frontiers"

def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

# =====================================================================
# FIGURE 1 - band layout
# =====================================================================
fig, ax = plt.subplots(figsize=(180*MM, 120*MM))
ax.set_xlim(0, 180); ax.set_ylim(0, 120); ax.axis("off")

def box(x, y, w, h, fc, ec, lw=0.9, r=1.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw))

def txt(x, y, s, size=7.2, weight="normal", ha="center", va="center",
        color="black", ls=1.35):
    ax.text(x, y, s, fontsize=size, fontweight=weight, ha=ha, va=va,
            color=color, linespacing=ls)

def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                 mutation_scale=9, lw=0.9, color=OI["grey"]))

box(15, 105, 150, 13, "#F2F2F2", "black")
txt(90, 113.6, "Can conventional lipids support a transferable point-of-care DR risk tool?",
    8.3, "bold")
txt(90, 108.6, "premises: (a) causal support    (b) stable, interpretable discrimination    (c) transferability", 7.3)

bands = [
    dict(fc="#EAF3FB", ec=OI["blue"], y=81.5,
         tag="Layer 1\nCausal premise",
         dm="Two-sample MR: 6 lipids (GLGC/UKB) vs FinnGen R13 DR\nIVW, Egger, weighted median, MVMR, PRESSO, cis-target MR",
         res="TG OR 0.86 (p=0.009) -> MVMR 0.94 (p=0.24);\n6 cis targets null; LPL-cis vs genome-wide discordant",
         v="NOT SUPPORTED", vc=OI["red"]),
    dict(fc="#FDF3E3", ec=OI["orange"], y=61.5,
         tag="Layer 2\nDiscrimination",
         dm="Chinese EMR cohort (AI-graded fundus), gradeable\nn=1,966 (275 events); logistic + LASSO, bootstrap 1,000",
         res="CC AUC 0.669 vs MI 0.549 (spec 0.55-0.63); IDI 0.032;\nApoB/LDL-C collinear, uninterpretable",
         v="FRAGILE", vc=OI["orange"]),
    dict(fc="#EAF7F1", ec=OI["green"], y=41.5,
         tag="Layer 3\nTransportability",
         dm="Chinese coefficients applied to NHANES 2005-08\ndiabetes subgroup (n=1,119); core / std / LASSO versions",
         res="AUC 0.477-0.572 (MEC-weighted 0.468-0.572);\ncalibration slope 0.31 / 0.38 / -0.12",
         v="NOT SUPPORTED", vc=OI["red"]),
    dict(fc="#F6ECF5", ec=OI["purple"], y=21.5,
         tag="Layer 4\nBenchmark",
         dm="NHANES full model: HbA1c, duration, insulin, eGFR,\nBMI, smoking, BP; survey-weighted 10-fold OOF",
         res="AUC 0.735 (temporal 0.748/0.722; weighted 0.749);\nBrier 0.176 vs lipid-only 0.505 (Brier 0.236)",
         v="ROBUST", vc=OI["green"]),
]
for b in bands:
    box(3, b["y"], 174, 17, b["fc"], b["ec"])
    txt(4, b["y"]+8.5, b["tag"], 7.1, "bold", ha="left")
    txt(31, b["y"]+8.5, b["dm"], 6.1, ha="left")
    txt(104, b["y"]+8.5, b["res"], 6.1, ha="left")
    txt(173.5, b["y"]+15.2, b["v"], 6.8, "bold", ha="right", color=b["vc"])
for y in (81.5, 61.5, 41.5):
    arrow(90, y, 90, y-3.0)

box(10, 6.0, 160, 12.5, "#FFFFFF", OI["verm"])
txt(90, 15.3, "Two generalisable failure modes quantified", 8, "bold", color=OI["verm"])
txt(90, 3.2, "positive replacement factors: glycaemic control (HbA1c), diabetes duration, treatment intensity", 7.2)
txt(90, 10.6, "(i) complete-case analysis inflates discrimination (ΔAUC 0.03–0.12; missingness DR-associated, OR 2.31)\n"
              "(ii) MR causal signals are instrument-dependent (genome-wide 0.86 vs LPL-cis 1.267 — discordant, neither robust)",
    6.5)
arrow(90, 21.5, 90, 18.8)
ax.set_ylim(0, 122)
fig.savefig(f"{OUT}/Figure1_study_design.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# =====================================================================
# FIGURE 2 - forest, annotations inside axes whitespace
# =====================================================================
pA = read_csv(r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/阶段6_论文整合/fig2_panelA_univariable.csv")
mvmr = read_csv(r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/阶段6_论文整合/fig2_panelB_mvmr.csv")
drug = read_csv(r"E:/workbuddy/2026-09-04-18-19-16/phase1_out/druggable_mr_results.csv")

fig, axes = plt.subplots(1, 2, figsize=(180*MM, 92*MM),
                         gridspec_kw={"width_ratios": [1.08, 1], "wspace": 0.30})

def forest(ax, rows, labels, colors, xlabel, title, annot_x=3.1):
    n = len(rows)
    ys = np.arange(n)[::-1]
    for y, r, c in zip(ys, rows, colors):
        lo, orv, hi = r["lo"], r["or"], r["hi"]
        ax.plot([lo, hi], [y, y], lw=1.1, color=c, solid_capstyle="butt", zorder=2)
        for e in (lo, hi):
            ax.plot([e, e], [y-0.15, y+0.15], lw=1.1, color=c, zorder=2)
        ax.plot(orv, y, "s", ms=4.2, color=c, zorder=3)
        star = " *" if r.get("sig") else ""
        ax.text(annot_x, y, f'{orv:.2f} ({lo:.2f}-{hi:.2f}){star}', fontsize=6.4,
                va="center", ha="left")
    ax.axvline(1, lw=0.8, ls="--", color=OI["grey"], zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=7.2)
    ax.set_xscale("log"); ax.set_xlim(0.45, 8.0)
    ax.set_xticks([0.5, 0.7, 1.0, 1.5, 2.0])
    ax.set_xticklabels(["0.5", "0.7", "1.0", "1.5", "2.0"])
    ax.set_xlabel(xlabel, fontsize=7.5)
    ax.minorticks_off()
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title(title, fontsize=8, pad=5)

rows, labels, colors = [], [], []
for r in pA:
    rows.append({"lo": float(r["or_lci95"]), "or": float(r["or"]),
                 "hi": float(r["or_uci95"]), "sig": float(r["pval"]) < 0.05})
    labels.append(r["lab"]); colors.append(OI["blue"])
for r in mvmr:
    rows.append({"lo": float(r["or_lci95"]), "or": float(r["or"]),
                 "hi": float(r["or_uci95"]), "sig": False})
    labels.append(f'{r["lab"]} (MVMR)'); colors.append(OI["sky"])
forest(axes[0], rows, labels, colors,
       "OR for DR per 1-SD genetically predicted lipid (95% CI, log scale)",
       "Univariable IVW and multivariable MR")

rows, labels, colors = [], [], []
for r in drug:
    rows.append({"lo": float(r["or_lci95"]), "or": float(r["OR"]),
                 "hi": float(r["or_uci95"]), "sig": False})
    labels.append(f'{r["target"]} ({r["lipid_trait"]})'); colors.append(OI["verm"])
forest(axes[1], rows, labels, colors,
       "OR per 1-SD lipid, cis drug-target IVs (95% CI, log scale)",
       "cis drug-target MR (6 lipid-lowering targets)")
axes[1].text(7.8, -1.15, "* p<0.05", fontsize=6.4, ha="right", color=OI["grey"])

for ax, lab in zip(axes, ["(a)", "(b)"]):
    ax.text(0.5, -0.34, lab, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="center", va="top")
fig.savefig(f"{OUT}/Figure2_MR_forest.png", bbox_inches="tight", facecolor="white")
plt.close(fig)
print("fig1+fig2 rebuilt")
