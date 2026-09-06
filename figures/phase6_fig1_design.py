# Figure 1: four-layer falsification design flowchart (publication style)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = "E:/workbuddy/2026-09-04-18-19-16/phase6_out/fig1_study_design.png"

BLUE, RED, GREEN, AMBER, GREY = "#2C7FB8", "#C0392B", "#1E8449", "#B9770E", "#555555"
FILL = {"r": "#FDEDEC", "g": "#EAF7EE", "y": "#FEF9E7", "b": "#E8F4FA"}

rows = [
    ("1  CAUSAL\nMendelian\nrandomization",
     "Two-sample MR (IVW / MR-Egger / weighted median)\n"
     "+ multivariable MR (TG + LDL-C + HDL-C mutually adjusted)\n"
     "6 lipid GWAS (GLGC + UKB Neale) \u2192 FinnGen R13 DR\n"
     "(15,353 cases / 88,091 controls)",
     "\u2717  No independent causal effect of any lipid\n"
     "TG: OR 0.86 (p=0.009) \u2192 0.94 (p=0.24) after MVMR", "r"),
    ("2  DISCRIMINATION\nChinese DR\nscreening cohort",
     "Bengbu cohort 2019\u20132021; n=2,482 \u2192 gradable 1,966 \u2192\n"
     "complete cases 1,207 (121 events)\n"
     "Logistic + LASSO; bootstrap 1,000 optimism correction",
     "\u25b3  AUC 0.669 (corrected) but severe collinearity\n"
     "(ApoB +5.03 vs LDL-C \u22121.25); lipid IDI only 0.032", "y"),
    ("3  TRANSPORTABILITY\nCross-ethnicity\nexternal validation",
     "NHANES 2005\u20132008, US diabetes subgroup n=1,119\n"
     "Chinese coefficients applied (units converted)\n"
     "lipid-only score re-fitted in NHANES as comparator",
     "\u2717  AUC 0.48\u20130.53 (\u2248 random); calibration slope\n"
     "0.31\u20130.38 (one version negative) \u2014 fails to generalize", "r"),
    ("4  BENCHMARK\nFull reference\nmodel",
     "NHANES diabetes subgroup n=888; survey-weighted;\n"
     "10-fold out-of-fold predictions\n"
     "HbA1c + diabetes duration + insulin therapy + covariates",
     "\u2713  AUC 0.735 (temporal validation 0.748 / 0.722)\n"
     "lipid-only AUC 0.505 \u2192 \u0394AUC = 0.230", "g"),
]

fig, ax = plt.subplots(figsize=(13.5, 7.6), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

ys = [85.5, 62.5, 39.5, 16.5]
H = 19.5
for (tag, mid, res, key), y in zip(rows, ys):
    y0 = y - H / 2
    # left tag bar
    ax.add_patch(FancyBboxPatch((1, y0), 15.5, H, boxstyle="round,pad=0.25",
                                fc=BLUE, ec="none", alpha=0.92))
    ax.text(8.75, y, tag, ha="center", va="center", fontsize=8.6,
            color="white", weight="bold", linespacing=1.5)
    # middle design box
    ax.add_patch(FancyBboxPatch((18.5, y0), 43.5, H, boxstyle="round,pad=0.25",
                                fc=FILL["b"], ec=BLUE, lw=0.9))
    ax.text(40.25, y, mid, ha="center", va="center", fontsize=7.6,
            color="#1a1a1a", linespacing=1.55)
    # arrow
    ax.annotate("", xy=(64.6, y), xytext=(62.6, y),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color=GREY))
    # result box
    rfc = {"r": FILL["r"], "y": FILL["y"], "g": FILL["g"]}[key]
    rec = {"r": RED, "y": AMBER, "g": GREEN}[key]
    ax.add_patch(FancyBboxPatch((65.2, y0), 33.8, H, boxstyle="round,pad=0.25",
                                fc=rfc, ec=rec, lw=1.1))
    ax.text(82.1, y, res, ha="center", va="center", fontsize=7.7,
            color=rec, weight="bold", linespacing=1.6)

ax.text(50, 97.5,
        "Lipid-based bedside DR risk scores: a four-layer falsification design",
        ha="center", va="center", fontsize=12.5, weight="bold")
ax.text(50, 3.2,
        "Layers 1\u20133 jointly reject the premise that lipids can support a generalizable "
        "bedside DR risk tool; Layer 4 identifies the factors that do (glycemic control, "
        "duration, treatment intensity).",
        ha="center", va="center", fontsize=8, color=GREY, style="italic")

plt.savefig(OUT, bbox_inches="tight", facecolor="white")
print("fig1_study_design.png written")
