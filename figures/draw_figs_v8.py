# v8 (2026-09-28): reads the corrected-outcome (option A) outputs in phase4_out_v2.
# Changes vs v7: Fig3/Fig4 inputs, Fig3 weighted slope annotations, Fig4 transported AUC label,
# Fig1 schematic figures (layer 2/3/4 numbers and the reframed methodological warning).
# -*- coding: utf-8 -*-
"""Frontiers-compliant figures v6 (final) for the lipid-DR manuscript.

Figure 1 : vertical 4-row layout, no in-figure main title, >=8 pt, MI 0.556-0.559.
Figure 2 : OR labels INSIDE each panel right-aligned (CI shown as error bars,
           full values in Table 2); unified log axis 0.45-4.0; short exposure
           labels; no stray "* p<0.05" text.
Figure 3 : MEC-weighted calibration slope annotated next to the unweighted one.
Figure 4 : legend states CHINESE coefficients, split over two lines so it fits
           the single-column width; no redundant panel letter.
Exact widths: 180 mm (double column) / 85 mm (single column), 300 dpi, Arial 8-9 pt.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import csv, os, shutil

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "savefig.dpi": 300, "figure.dpi": 150,
})
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "verm": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9",
      "yellow": "#F0E442", "black": "#000000", "grey": "#7F7F7F"}

MM = 1 / 25.4
P6 = r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/阶段6_论文整合"
P1 = r"E:/workbuddy/2026-09-04-18-19-16/phase1_out"
P4 = r"E:/workbuddy/2026-09-04-18-19-16/phase4_out_v2"
OUT = r"E:/workbuddy/2026-09-04-18-19-16/figs_frontiers"
DEST_PNG = r"E:/workbuddy/Front  in  Endocrinology  投稿"
DEST_TIFF = [r"E:/workbuddy/Front  in  Endocrinology  投稿",
             r"E:/workbuddy/Front  in  Endocrinology  投稿/手稿1_图片"]
os.makedirs(OUT, exist_ok=True)


def read_csv(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


# =====================================================================
# FIGURE 1 - study design, vertical 4-row layout (180 x 152 mm)
# =====================================================================
fig = plt.figure(figsize=(180 * MM, 152 * MM))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 180)
ax.set_ylim(0, 152)
ax.axis("off")


def box(x, y, w, h, fc, ec, lw=0.9, r=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=%d" % r,
                                fc=fc, ec=ec, lw=lw))


def txt(x, y, s, size=8, weight="normal", ha="center", va="center", color="black"):
    ax.text(x, y, s, fontsize=size, fontweight=weight, ha=ha, va=va,
            color=color, linespacing=1.35)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=9, lw=0.9, color=OI["grey"]))


ROW_H = 22.0
rows = [
    dict(y=127.0, tag="Layer 1", name="Causal premise", fc="#EAF3FB", ec=OI["blue"],
         l1="Two-sample MR: 6 lipids (GLGC/UKB) \u2192 FinnGen R13 DR (15,353/88,091)",
         l2="IVW (main), MR-Egger, weighted median, MVMR, MR-PRESSO, cis-MR",
         l3="TG OR 0.86 (p=0.009) \u2192 MVMR 0.94 (p=0.24); all 6 cis targets null",
         v="\u2717  No independent causal effect"),
    dict(y=102.0, tag="Layer 2", name="Discrimination", fc="#FDF3E3", ec=OI["orange"],
         l1="Chinese EMR screening cohort (graded fundus); n=1,971; 474 events (24.0%)",
         l2="Logistic + LASSO, bootstrap 1,000; multiple imputation (m=50, maxit=50)",
         l3="CC 0.579 vs MI 0.582; an earlier 0.06 gap was an outcome-coding artefact",
         v="\u25b3  Fragile; coefficients uninterpretable"),
    dict(y=77.0, tag="Layer 3", name="Transportability", fc="#EAF7F1", ec=OI["green"],
         l1="Chinese coefficients applied to the NHANES 2005\u20132008 diabetes subgroup",
         l2="core / std / LASSO versions; unweighted and MEC-weighted analyses",
         l3="AUC 0.498 / 0.493 / 0.551 (at chance); weighted slope \u22120.91 / \u22121.81 / +2.02",
         v="\u2717  Fails to generalise"),
    dict(y=52.0, tag="Layer 4", name="Benchmark", fc="#F6ECF5", ec=OI["purple"],
         l1="NHANES full reference model: HbA1c, duration, insulin use, eGFR\u2026",
         l2="Survey-weighted logistic, 10-fold out-of-fold; temporal validation",
         l3="AUC 0.733 vs locally fitted lipid-only 0.515 (transported 0.486)",
         v="\u2713  HbA1c, duration, insulin use robust"),
]
for L in rows:
    y = L["y"]
    box(2, y, 40, ROW_H, L["fc"], L["ec"])
    txt(22, y + 13.8, L["tag"], 8.5, "bold")
    txt(22, y + 8.0, L["name"], 8)
    box(44, y, 134, ROW_H, "#FFFFFF", L["ec"])
    txt(47, y + 17.8, L["l1"], 8, ha="left")
    txt(47, y + 12.6, L["l2"], 8, ha="left")
    txt(47, y + 7.4, L["l3"], 8, ha="left")
    txt(47, y + 2.6, L["v"], 8, "bold", ha="left", color=L["ec"])

box(18, 35, 144, 13, "#F5F5F5", OI["black"])
txt(90, 44.3, "All three premises fail \u2014 lipids cannot support a transferable bedside DR risk tool",
    8, "bold")
txt(90, 38.8, "replacement factors: glycaemic control (HbA1c), diabetes duration, treatment intensity", 8)
arrow(90, 51.6, 90, 48.4)

box(8, 2, 164, 30, "#FFFFFF", OI["verm"])
txt(90, 27.6, "Two generalisable sources of analytical dependence", 8, "bold", color=OI["verm"])
txt(90, 21.0, "(i) an outcome-coding error mimicked a missing-data effect; correcting it removed a 0.06 gap", 8)
txt(90, 14.6, "(ii) MR causal signals are instrument-dependent (genome-wide 0.86 vs LPL-cis 1.267)", 8)
arrow(90, 34.6, 90, 32.4)

fig.savefig("%s/Figure1_study_design.png" % OUT, facecolor="white")
plt.close(fig)

# =====================================================================
# FIGURE 2 - MR forest, 2 panels (180 x 95 mm)
# =====================================================================
pA = read_csv("%s/fig2_panelA_univariable.csv" % P6)
mvmr = read_csv("%s/fig2_panelB_mvmr.csv" % P6)
drug = read_csv("%s/druggable_mr_results.csv" % P1)
SHORT = {"Total cholesterol": "TC"}

fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 95 * MM))
fig.subplots_adjust(left=0.135, right=0.975, top=0.855, bottom=0.215, wspace=0.36)


def forest(ax, rows_, labels, colors, xlabel):
    n = len(rows_)
    ys = np.arange(n)[::-1]
    for y, r, c in zip(ys, rows_, colors):
        ax.plot([r["lo"], r["hi"]], [y, y], lw=1.1, color=c, solid_capstyle="butt", zorder=2)
        for e in (r["lo"], r["hi"]):
            ax.plot([e, e], [y - 0.14, y + 0.14], lw=1.1, color=c, zorder=2)
        ax.plot(r["or"], y, "s", ms=4.2, color=c, zorder=3)
        star = " *" if r.get("sig") else ""
        ax.text(0.975, y, "%.2f%s" % (r["or"], star),
                fontsize=8, va="center", ha="right", transform=ax.get_yaxis_transform())
    ax.axvline(1, lw=0.8, ls="--", color=OI["grey"], zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(0.45, 4.0)
    ax.set_xticks([0.5, 0.7, 1.0, 1.5, 2.0, 3.0])
    ax.set_xticklabels(["0.5", "0.7", "1.0", "1.5", "2.0", "3.0"])
    ax.set_xlabel(xlabel, fontsize=8)
    ax.minorticks_off()
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)


rows_, labels, colors = [], [], []
for r in pA:
    rows_.append({"lo": float(r["or_lci95"]), "or": float(r["or"]),
                  "hi": float(r["or_uci95"]), "sig": float(r["pval"]) < 0.05})
    labels.append(SHORT.get(r["lab"], r["lab"]))
    colors.append(OI["blue"])
for r in mvmr:
    rows_.append({"lo": float(r["or_lci95"]), "or": float(r["or"]),
                  "hi": float(r["or_uci95"]), "sig": False})
    labels.append("%s (MVMR)" % r["lab"])
    colors.append(OI["sky"])
forest(axes[0], rows_, labels, colors,
       "OR for DR per 1-SD (GLGC) or per allele\n(UKB ApoA-I/ApoB); 95% CI")
axes[0].set_title("Univariable IVW and multivariable MR", fontsize=9, pad=5)

rows_, labels, colors = [], [], []
for r in drug:
    rows_.append({"lo": float(r["or_lci95"]), "or": float(r["OR"]),
                  "hi": float(r["or_uci95"]), "sig": False})
    labels.append("%s (%s)" % (r["target"], r["lipid_trait"]))
    colors.append(OI["verm"])
forest(axes[1], rows_, labels, colors,
       "OR for DR per 1-SD TG via cis\ndrug-target IVs; 95% CI")
axes[1].set_title("cis drug-target MR (6 lipid-lowering targets)", fontsize=9, pad=5)

for a, lab in zip(axes, ["(A)", "(B)"]):
    a.text(0.5, -0.30, lab, transform=a.transAxes, fontsize=9,
           fontweight="bold", ha="center", va="top")
fig.savefig("%s/Figure2_MR_forest.png" % OUT, facecolor="white")
plt.close(fig)

# =====================================================================
# FIGURE 3 - calibration, 3 panels (180 x 68 mm)
# =====================================================================
cal_files = {"core": "%s/phase4_calib_core.csv" % P4,
             "std": "%s/phase4_calib_std.csv" % P4,
             "lasso": "%s/phase4_calib_lasso.csv" % P4}
ext = {r["model"]: r for r in read_csv("%s/phase4_external_auc.csv" % P4)}
WS = {"core": -0.91, "std": -1.81, "lasso": 2.02}

fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 68 * MM), sharey=True)
fig.subplots_adjust(left=0.10, right=0.985, top=0.845, bottom=0.225, wspace=0.17)
cols = [OI["blue"], OI["orange"], OI["verm"]]
for ax, (nm, path), c, lab in zip(axes, cal_files.items(), cols, ["(A)", "(B)", "(C)"]):
    rws = read_csv(path)
    pm = np.array([float(r["pred_mean"]) for r in rws])
    ob = np.array([float(r["obs_rate"]) for r in rws])
    lim = 0.65
    ax.plot([0, lim], [0, lim], ls="--", lw=0.9, color=OI["grey"], zorder=1)
    ax.scatter(pm, ob, s=14, color=c, zorder=3)
    ax.plot(pm, ob, lw=0.9, color=c, zorder=2)
    r = ext[nm]
    ax.text(0.965, 0.035,
            "AUC %.3f\nslope %.2f  (MEC-weighted %.2f)\nintercept %.2f"
            % (float(r["auc"]), float(r["calib_slope"]), WS[nm], float(r["calib_intercept"])),
            transform=ax.transAxes, fontsize=8, va="bottom", ha="right",
            bbox=dict(fc="white", ec="none", alpha=0.88, pad=1.2))
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Predicted risk", fontsize=8)
    if nm == "core":
        ax.set_ylabel("Observed any-DR rate", fontsize=8)
    ax.set_title({"core": "core (age, sex, TC, HDL-C)",
                  "std": "std (+TG, LDL-C)",
                  "lasso": "LASSO (+ApoB)"}[nm], fontsize=9, pad=5)
    ax.text(0.5, -0.30, lab, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="center", va="top")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.savefig("%s/Figure3_calibration.png" % OUT, facecolor="white")
plt.close(fig)

# =====================================================================
# FIGURE 4 - ROC (85 x 88 mm, single column)
# =====================================================================
roc = read_csv("%s/phase4f_roc3.csv" % P4)
fpr = np.array([float(r["fpr"]) for r in roc])
tf = np.array([float(r["tpr_full"]) for r in roc])
tl = np.array([float(r["tpr_local"]) for r in roc])
tc = np.array([float(r["tpr_chinese"]) for r in roc])

fig, ax = plt.subplots(figsize=(85 * MM, 88 * MM))
fig.subplots_adjust(left=0.205, right=0.975, top=0.975, bottom=0.135)
ax.plot(fpr, tf, lw=1.4, color=OI["purple"], label="Full reference model\n(AUC 0.733)")
ax.plot(fpr, tl, lw=1.4, color=OI["orange"], label="Local lipid-only model\n(AUC 0.515)")
ax.plot(fpr, tc, lw=1.4, ls="--", color=OI["blue"],
        label="Transported Chinese score\n(AUC 0.486)")
ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color=OI["grey"])
ax.set_xlabel("1 \u2212 specificity", fontsize=8)
ax.set_ylabel("Sensitivity", fontsize=8)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.005)
ax.legend(loc="lower right", fontsize=8, frameon=False)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.savefig("%s/Figure4_ROC.png" % OUT, facecolor="white")
plt.close(fig)

# =====================================================================
# copy PNGs + build 300-dpi RGB LZW TIFFs
# =====================================================================
NAMES = {"Figure1_study_design.png": "Figure1",
         "Figure2_MR_forest.png": "Figure2",
         "Figure3_calibration.png": "Figure3",
         "Figure4_ROC.png": "Figure4"}
log = []
from PIL import Image
for png, base in NAMES.items():
    src = "%s/%s" % (OUT, png)
    shutil.copy2(src, "%s/%s" % (DEST_PNG, png))
    im = Image.open(src)
    log.append("%s : %dx%d px = %.1f x %.1f mm @300dpi"
               % (base, im.size[0], im.size[1], im.size[0] / 300 * 25.4, im.size[1] / 300 * 25.4))
    im = im.convert("RGB")
    for d in DEST_TIFF:
        im.save("%s/%s.tiff" % (d, base), format="TIFF",
                compression="tiff_lzw", dpi=(300, 300))

with open("%s/_render.log" % OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("OK")
