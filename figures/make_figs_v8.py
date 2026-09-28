# Generate draw_figs_v8.py from draw_figs_v7.py:
#   - point Fig3/Fig4 at the corrected-outcome outputs (phase4_out_v2)
#   - update the MEC-weighted calibration slopes annotated on Fig3
#   - update the Fig4 transported-score AUC label
#   - update the numbers baked into the Fig1 schematic (including the reframed warning)
import io, sys
sys.stdout.reconfigure(encoding="utf-8")

SRC = r"E:/workbuddy/2026-09-04-18-19-16/figs_frontiers/draw_figs_v7.py"
DST = r"E:/workbuddy/2026-09-04-18-19-16/figs_frontiers/draw_figs_v8.py"

s = io.open(SRC, encoding="utf-8").read()
n = 0


def rep(old, new):
    global s, n
    assert old in s, "ANCHOR NOT FOUND: " + old[:90]
    s = s.replace(old, new)
    n += 1


# ---- paths: Fig3 / Fig4 read the corrected-outcome outputs ----
rep('P4 = r"E:/workbuddy/2026-09-04-18-19-16/phase4_out"',
    'P4 = r"E:/workbuddy/2026-09-04-18-19-16/phase4_out_v2"')

# ---- Fig3 annotated weighted slopes ----
rep('WS = {"core": 0.79, "std": 0.94, "lasso": -0.17}',
    'WS = {"core": -0.91, "std": -1.81, "lasso": 2.02}')

# ---- Fig1 schematic: layer 2 ----
rep('l1="Chinese EMR screening cohort (graded fundus); n=1,966; 275 events (14.0%)",',
    'l1="Chinese EMR screening cohort (graded fundus); n=1,971; 474 events (24.0%)",')
rep('l3="CC AUC 0.668 vs MI 0.608 (same correction); IPW: not selection bias",',
    'l3="CC 0.579 vs MI 0.582; an earlier 0.06 gap was an outcome-coding artefact",')

# ---- Fig1 schematic: layer 3 ----
rep('l3="AUC 0.513 / 0.530 / 0.477 (all at chance); slope 0.31 / 0.38 / \\u22120.12",',
    'l3="AUC 0.498 / 0.493 / 0.551 (at chance); weighted slope \\u22120.91 / \\u22121.81 / +2.02",')

# ---- Fig1 schematic: layer 4 ----
rep('l3="AUC 0.733 vs locally fitted lipid-only 0.515 (transported 0.505)",',
    'l3="AUC 0.733 vs locally fitted lipid-only 0.515 (transported 0.486)",')

# ---- Fig1 schematic: bottom warning box ----
rep('txt(90, 21.0, "(i) missing-data handling shifts discrimination (dAUC ~ 0.06); IPW shows it is not selection", 8)',
    'txt(90, 21.0, "(i) an outcome-coding error mimicked a missing-data effect; correcting it removed a 0.06 gap", 8)')

# ---- Fig4 ROC legend label ----
rep('"Transported Chinese score\\n(AUC 0.505)")',
    '"Transported Chinese score\\n(AUC 0.486)")')

s = ("# v8 (2026-09-28): reads the corrected-outcome (option A) outputs in phase4_out_v2.\n"
     "# Changes vs v7: Fig3/Fig4 inputs, Fig3 weighted slope annotations, Fig4 transported AUC label,\n"
     "# Fig1 schematic figures (layer 2/3/4 numbers and the reframed methodological warning).\n"
     + s)

io.open(DST, "w", encoding="utf-8").write(s)
print("wrote", DST, "| replacements:", n)
