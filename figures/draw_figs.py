# -*- coding: utf-8 -*-
"""Frontiers-compliant figures for the lipid-DR manuscript (v3).
Style: Arial 8-9 pt, 300 dpi, Okabe-Ito colour-blind-safe palette,
no in-figure main title (caption external), lowercase panel letters
centred below each panel, full width 180 mm / single column 85 mm."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import csv, os

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],  # per-glyph fallback: ✓ ✗ △ → from DejaVu
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "savefig.dpi": 300, "figure.dpi": 150,
})
# Okabe-Ito
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "verm": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9",
      "yellow": "#F0E442", "black": "#000000", "grey": "#7F7F7F"}
MM = 1 / 25.4
OUT = r"E:/workbuddy/2026-09-04-18-19-16/figs_frontiers"
os.makedirs(OUT, exist_ok=True)

def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

# =====================================================================
# FIGURE 1 - study design flowchart (full width 180 mm)
# =====================================================================
fig, ax = plt.subplots(figsize=(180*MM, 118*MM))
ax.set_xlim(0, 180); ax.set_ylim(0, 118); ax.axis("off")

def box(x, y, w, h, fc, ec, lw=0.9, r=2):
    b = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw)
    ax.add_patch(b)

def txt(x, y, s, size=7.5, weight="normal", ha="center", va="center", color="black"):
    ax.text(x, y, s, fontsize=size, fontweight=weight, ha=ha, va=va, color=color,
            linespacing=1.45)

def arrow(x0, y0, x1, y1):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                        mutation_scale=9, lw=0.9, color=OI["grey"])
    ax.add_patch(a)

# top question box
box(22, 100, 136, 14, "#F5F5F5", OI["black"])
txt(90, 109.8, "Can conventional lipids support a transferable point-of-care DR risk tool?",
    8.5, "bold")
txt(90, 104.6, "premises: (a) causal support   (b) stable, interpretable discrimination   (c) transferability",
    7.5)

layers = [
    dict(x=1.5,  tag="Layer 1", name="Causal premise", fc="#EAF3FB", ec=OI["blue"],
         l1="Two-sample MR: 6 lipids\n(GLGC/UKB) → FinnGen R13 DR\n(15,353/88,091)",
         l2="IVW, MR-Egger, weighted median,\nMVMR, MR-PRESSO, cis-MR",
         l3="TG OR 0.86 (p=0.009) → MVMR 0.94;\n6 cis drug targets all null",
         v="✗ No independent causal effect"),
    dict(x=45.5, tag="Layer 2", name="Discrimination", fc="#FDF3E3", ec=OI["orange"],
         l1="Chinese EMR screening cohort\n(manually graded fundus), n=1,966",
         l2="Logistic + LASSO, bootstrap 1,000;\nMI (m=50, maxit=50) as primary",
         l3="CC AUC 0.669 vs MI 0.556–0.559\n(spec 0.556–0.63); IDI 0.032",
         v="△ Fragile; uninterpretable"),
    dict(x=89.5, tag="Layer 3", name="Transportability", fc="#EAF7F1", ec=OI["green"],
         l1="Chinese coefficients → NHANES\n2005–2008 diabetes subgroup",
         l2="core / std / LASSO versions;\nunweighted + MEC-weighted",
         l3="AUC 0.513/0.530/0.477 (chance);\nslope 0.31/0.38/−0.12\n(weighted 0.79/0.94/−0.17)",
         v="✗ Fails to generalise"),
    dict(x=133.5, tag="Layer 4", name="Benchmark", fc="#F6ECF5", ec=OI["purple"],
         l1="NHANES full reference model:\nHbA1c, duration, insulin use, eGFR…",
         l2="Weighted logistic, 10-fold OOF,\nbi-directional temporal validation",
         l3="AUC 0.735 (temporal 0.748/0.722);\nvs lipid-only 0.505",
         v="✓ HbA1c, duration, insulin use"),
]
for L in layers:
    box(L["x"], 46, 43, 48, L["fc"], L["ec"])
    txt(L["x"]+21.5, 89.5, f'{L["tag"]} · {L["name"]}', 7.5, "bold")
    txt(L["x"]+21.5, 79.5, L["l1"], 5.6)
    txt(L["x"]+21.5, 68.5, L["l2"], 5.6)
    txt(L["x"]+21.5, 57.0, L["l3"], 5.6)
    txt(L["x"]+21.5, 49.2, L["v"], 6.1, "bold")
    arrow(L["x"]+21.5, 46, L["x"]+21.5, 40.5)

# convergence box
box(20, 26, 140, 14, "#F5F5F5", OI["black"])
txt(90, 35.8, "All three premises fail — lipids cannot support a transferable bedside DR risk tool",
    7.5, "bold")
txt(90, 30.3, "positive replacement factors: glycaemic control (HbA1c), diabetes duration, treatment intensity",
    6.8)
arrow(90, 26, 90, 24.5)

box(10, 5, 160, 19, "#FFFFFF", OI["verm"])
txt(90, 20.5, "Two generalisable failure modes quantified", 8, "bold", color=OI["verm"])
txt(90, 14.5, "(i) complete-case analysis inflates discrimination (ΔAUC 0.03–0.12; missingness DR-associated, OR 2.31)",
    6.5)
txt(90, 9.5, "(ii) MR causal signals are instrument-dependent (genome-wide 0.86 vs LPL-cis 1.267 — discordant, neither robust)",
    6.5)
# arrows from top to layers
for L in layers:
    arrow(90, 100, L["x"]+21.5, 94.2)
fig.savefig(f"{OUT}/Figure1_study_design.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# =====================================================================
# FIGURE 2 - MR forest, 2 panels (full width)
# =====================================================================
pA = read_csv(r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/阶段6_论文整合/fig2_panelA_univariable.csv")
mvmr = read_csv(r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/阶段6_论文整合/fig2_panelB_mvmr.csv")
drug = read_csv(r"E:/workbuddy/2026-09-04-18-19-16/phase1_out/druggable_mr_results.csv")

fig, axes = plt.subplots(1, 2, figsize=(180*MM, 95*MM),
                         gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.52})

def forest(ax, rows, labels, colors, xlabel, xrefs=(0.5, 1.0, 2.0)):
    n = len(rows)
    ys = np.arange(n)[::-1]
    for y, r, c, lb in zip(ys, rows, colors, labels):
        lo, orv, hi = r["lo"], r["or"], r["hi"]
        sig = r.get("sig", False)
        ax.plot([lo, hi], [y, y], lw=1.1, color=c, solid_capstyle="butt", zorder=2)
        for e in (lo, hi):
            ax.plot([e, e], [y-0.14, y+0.14], lw=1.1, color=c, zorder=2)
        ax.plot(orv, y, "s", ms=4.2, color=c, zorder=3)
        star = " *" if sig else ""
        # right-aligned INSIDE the axes so OR labels never collide with the
        # neighbouring panel's y-tick labels (xlim extended to 3.5 for headroom)
        ax.text(0.99, y, f'{orv:.2f} ({lo:.2f}–{hi:.2f}){star}', fontsize=6.4,
                va="center", ha="right", transform=ax.get_yaxis_transform(),
                family="Arial")
    ax.axvline(1, lw=0.8, ls="--", color=OI["grey"], zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=7.2)
    ax.set_xscale("log"); ax.set_xlim(0.45, 5.5)
    ax.set_xticks([0.5, 0.7, 1.0, 1.5, 2.0])
    ax.set_xticklabels(["0.5", "0.7", "1.0", "1.5", "2.0"])
    ax.set_xlabel(xlabel, fontsize=7.5)
    ax.minorticks_off()
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)

# panel a: univariable IVW (6) + MVMR (3)
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
       "OR for DR per 1-SD genetically predicted lipid (95% CI, log scale)")
axes[0].set_title("Univariable IVW and multivariable MR", fontsize=8, pad=5)

# panel b: cis drug-target MR
rows, labels, colors = [], [], []
for r in drug:
    rows.append({"lo": float(r["or_lci95"]), "or": float(r["OR"]),
                     "hi": float(r["or_uci95"]), "sig": False})
    labels.append(f'{r["target"]} ({r["lipid_trait"]})'); colors.append(OI["verm"])
forest(axes[1], rows, labels, colors,
       "OR for DR per 1-SD lipid via cis drug-target IVs")
axes[1].set_title("cis drug-target MR (6 lipid-lowering targets)", fontsize=8, pad=5)
axes[1].text(0.99, -0.16, "* p<0.05", transform=axes[1].transAxes,
             fontsize=6.4, ha="right", color=OI["grey"])

# panel letters centred below panels
for ax, lab in zip(axes, ["(a)", "(b)"]):
    ax.text(0.5, -0.30, lab, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="center", va="top")
fig.savefig(f"{OUT}/Figure2_MR_forest.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# =====================================================================
# FIGURE 3 - calibration, 3 panels (full width)
# =====================================================================
cal_files = {
    "core": r"E:/workbuddy/2026-09-04-18-19-16/phase4_out/phase4_calib_core.csv",
    "std": r"E:/workbuddy/2026-09-04-18-19-16/phase4_out/phase4_calib_std.csv",
    "lasso": r"E:/workbuddy/2026-09-04-18-19-16/phase4_out/phase4_calib_lasso.csv",
}
ext = {r["model"]: r for r in read_csv(r"E:/workbuddy/2026-09-04-18-19-16/phase4_out/phase4_external_auc.csv")}

fig, axes = plt.subplots(1, 3, figsize=(180*MM, 66*MM), sharey=True,
                         gridspec_kw={"wspace": 0.18})
cols = [OI["blue"], OI["orange"], OI["verm"]]
for ax, (nm, path), c, lab in zip(axes, cal_files.items(), cols,
                                  ["(a) core", "(b) std", "(c) LASSO"]):
    rows = read_csv(path)
    pm = np.array([float(r["pred_mean"]) for r in rows])
    ob = np.array([float(r["obs_rate"]) for r in rows])
    nn = np.array([float(r["n"]) for r in rows])
    lim = 0.65
    ax.plot([0, lim], [0, lim], ls="--", lw=0.9, color=OI["grey"], zorder=1)
    ax.scatter(pm, ob, s=14, color=c, zorder=3)
    ax.plot(pm, ob, lw=0.9, color=c, zorder=2)
    r = ext[nm]
    ax.text(0.04, 0.93, f'AUC {float(r["auc"]):.3f}\nslope {float(r["calib_slope"]):.2f}\nintercept {float(r["calib_intercept"]):.2f}',
            transform=ax.transAxes, fontsize=6.8, va="top",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Predicted risk", fontsize=7.5)
    if nm == "core":
        ax.set_ylabel("Observed any-DR rate", fontsize=7.5)
    ax.set_title({"core": "core (age, sex, TC, HDL-C)",
                  "std": "std (+TG, LDL-C)",
                  "lasso": "LASSO (+ApoB)"}[nm], fontsize=8, pad=5)
    ax.text(0.5, -0.34, lab, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="center", va="top")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.savefig(f"{OUT}/Figure3_calibration.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

# =====================================================================
# FIGURE 4 - ROC (single column 85 mm)
# =====================================================================
roc = read_csv(r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/阶段6_论文整合/fig4_roc_curves.csv")
fpr = np.array([float(r["fpr"]) for r in roc])
tf = np.array([float(r["tpr_full"]) for r in roc])
tl = np.array([float(r["tpr_lipid"]) for r in roc])

fig, ax = plt.subplots(figsize=(100*MM, 100*MM))
ax.plot(fpr, tf, lw=1.4, color=OI["purple"], label="Full model (AUC 0.735)")
ax.plot(fpr, tl, lw=1.4, color=OI["blue"], label="Lipid-only score (AUC 0.505)")
ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color=OI["grey"])
ax.set_xlabel("1 − specificity", fontsize=8)
ax.set_ylabel("Sensitivity", fontsize=8)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.005)
ax.legend(loc="lower right", fontsize=7, frameon=False)
ax.text(0.5, -0.30, "(a)", transform=ax.transAxes, fontsize=9,
        fontweight="bold", ha="center", va="top")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.savefig(f"{OUT}/Figure4_ROC.png", bbox_inches="tight", facecolor="white")
plt.close(fig)

print("figures written:", os.listdir(OUT))
