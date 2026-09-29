# -*- coding: utf-8 -*-
"""Build dev_set_v2_dedup.csv: one record per patient identifier, retaining the most recent
   screening year (sensitivity analysis for repeat screening), then generate the _dedup R
   scripts by repointing the _v2 pipeline at that input.
"""
import sys, os, io, re
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
SRC = r"E:/workbuddy/2026-09-04-18-19-16"

aud = pd.read_csv(f"{SRC}/phase2_dev/severity_rebuild_audit.csv", low_memory=False)
dev = pd.read_csv(f"{SRC}/phase2_dev/dev_set_v2.csv", low_memory=False)
assert (dev["id"].values == aud["pid"].values).all(), "row order no longer aligned"

aud = aud.reset_index(drop=True)
d = aud.sort_values(["pid", "year"], kind="mergesort")
keep = sorted(d.groupby("pid", sort=False).tail(1).index)
print("records %d -> %d (patients %d)" % (len(aud), len(keep), aud.pid.nunique()))

out = dev.loc[keep].reset_index(drop=True)
out.to_csv(f"{SRC}/phase2_dev/dev_set_v2_dedup.csv", index=False, encoding="utf-8")

g = out[out.gradeable == 1]
cc = g[["age", "male", "tg", "tc", "hdl_c", "ldl_c", "apo_a1", "apo_b", "any_dr"]].dropna()
print("gradeable %d | events %d (%.1f%%) | PDR %d"
      % (len(g), int(g.any_dr.sum()), 100 * g.any_dr.mean(), int((g.severity == 4).sum())))
print("severity:", g.severity.value_counts().sort_index().to_dict())
print("complete cases %d | events %d (%.1f%%)" % (len(cc), int(cc.any_dr.sum()), 100 * cc.any_dr.mean()))
print("dedup file written:", out.shape)

# ---------------------------------------------------------------- generate _dedup scripts
STAGE = [("phase2_out_v2", "@@P2OV@@"), ("phase4_out_v2", "@@P4OV@@"), ("phase7_out_v2", "@@P7OV@@"),
         ("phase2_out", "@@P2O@@"), ("phase4_out", "@@P4O@@"), ("phase7_out", "@@P7O@@")]
BACK = {"@@P2OV@@": "phase2_out_v2d", "@@P4OV@@": "phase4_out_v2d", "@@P7OV@@": "phase7_out_v2d",
        "@@P2O@@": "phase2_out_v2d", "@@P4O@@": "phase4_out_v2d", "@@P7O@@": "phase7_out_v2d"}

for src, dst in [("phase2_v2.R", "phase2_v2_dedup.R"),
                 ("phase2_model_v2.R", "phase2_model_v2_dedup.R"),
                 ("phase4_external_v2.R", "phase4_external_v2_dedup.R"),
                 ("phase4e_local_v2.R", "phase4e_local_v2_dedup.R"),
                 ("phase7h_weighted_v2.R", "phase7h_weighted_v2_dedup.R")]:
    s = io.open(f"{SRC}/{src}", encoding="utf-8", errors="replace").read()
    for a, b in STAGE:
        s = s.replace(a, b)
    for a, b in BACK.items():
        s = s.replace(a, b)
    s = s.replace("dev_set_v2.csv", "dev_set_v2_dedup.csv")
    s = s.replace("phase_v2_", "phase_v2d_")
    s = s.replace("v2_phase7h_", "v2d_phase7h_")
    io.open(f"{SRC}/{dst}", "w", encoding="utf-8").write(
        "## AUTO-GENERATED from %s -- sensitivity analysis: one record per patient\n"
        "## (most recent screening year), input dev_set_v2_dedup.csv.\n\n" % src + s)
    print("generated", dst)
