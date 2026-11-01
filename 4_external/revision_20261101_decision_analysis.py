#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Decision analysis of the layer-4 models (Supplementary Table S48, Figure S2).

Reads the out-of-fold predictions written by
``4_external/revision_20261030_figure4_redraw.R`` and writes three aggregate tables into
``output_tables/``. The predictions themselves are individual-level and are never
committed (see POLICY_COMPLIANCE.md); pass their path with ``LDF_OOF``.

    LDF_OOF=/path/to/nhanes_oof_predictions.csv python revision_20261101_decision_analysis.py

Input columns: ``row, any_dr, wt4, fold, p_full, p_loc, p_tr`` on the primary frame
(n = 756, 267 events; ungradable photographs carry a missing outcome, so they are absent
from the frame rather than labelled non-cases).

Outputs (all survey-weighted by the 4-year MEC weight, WTMEC2YR/2):
  output_tables/dca_net_benefit_curve.csv   net benefit over thresholds 0.05-0.60,
                                            three models + refer-all, weighted and unweighted
  output_tables/dca_operating_points.csv    sensitivity, specificity, PPV and referral
                                            burden at fixed thresholds
  output_tables/dca_fixed_operating.csv     fixed fractions of the weighted population
                                            referred, and the referral needed for 90%
                                            sensitivity

Reading rule used in the manuscript: net benefit may not be compared with "refer all"
above the prevalence, because there the default is itself harmful; the lipid scores are
described as unable to beat that default within the range where it is reasonable.
"""
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "output_tables")


def net_benefit(p, y, w, t):
    pos = p >= t
    tp = w[pos & (y == 1)].sum() / w.sum()
    fp = w[pos & (y == 0)].sum() / w.sum()
    return tp - fp * t / (1 - t)


def operating_point(p, y, w, t):
    pos = p >= t
    wtp = w[pos & (y == 1)].sum() / w.sum()
    wfp = w[pos & (y == 0)].sum() / w.sum()
    wfn = w[~pos & (y == 1)].sum() / w.sum()
    wtn = w[~pos & (y == 0)].sum() / w.sum()
    tp = int((pos & (y == 1)).sum())
    fp = int((pos & (y == 0)).sum())
    fn = int((~pos & (y == 1)).sum())
    tn = int((~pos & (y == 0)).sum())
    return {
        "sens": tp / (tp + fn) if (tp + fn) > 0 else np.nan,
        "spec": tn / (tn + fp) if (tn + fp) > 0 else np.nan,
        "ppv": tp / (tp + fp) if (tp + fp) > 0 else np.nan,
        "referred": pos.mean(),
        "wsens": wtp / (wtp + wfn) if (wtp + wfn) > 0 else np.nan,
        "wspec": wtn / (wtn + wfp) if (wtn + wfp) > 0 else np.nan,
        "wreferred": wtp + wfp,
    }


def main():
    path = os.environ.get("LDF_OOF", "")
    if not path or not os.path.exists(path):
        sys.exit("LDF_OOF must point at nhanes_oof_predictions.csv written by "
                 "4_external/revision_20261030_figure4_redraw.R")

    d = pd.read_csv(path)
    y = d["any_dr"].to_numpy().astype(int)
    w = d["wt4"].to_numpy().astype(float)
    models = {
        "full_reference": d["p_full"].to_numpy(),
        "local_lipid": d["p_loc"].to_numpy(),
        "transported": d["p_tr"].to_numpy(),
    }
    os.makedirs(OUT, exist_ok=True)
    print(f"n = {len(y)}, events = {int(y.sum())}, weighted prevalence = "
          f"{w[y == 1].sum() / w.sum():.4f}")

    # ---- net-benefit curves
    ths = np.round(np.arange(0.05, 0.605, 0.01), 2)
    rows = []
    for name, p in models.items():
        for t in ths:
            rows.append({
                "model": name, "threshold": t,
                "nb_w": net_benefit(p, y, w, t),
                "nb_u": net_benefit(p, y, np.ones_like(w), t),
            })
    curve = pd.DataFrame(rows)
    allw = pd.DataFrame([{
        "threshold": t,
        "nb_all_w": net_benefit(np.ones_like(y, dtype=float), y, w, t),
        "nb_all_u": net_benefit(np.ones_like(y, dtype=float), y, np.ones_like(w), t),
    } for t in ths])
    curve = curve.merge(allw, on="threshold")
    curve.to_csv(os.path.join(OUT, "dca_net_benefit_curve.csv"), index=False)

    # ---- fixed thresholds
    fixed_t = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
    rows = []
    for name, p in models.items():
        for t in fixed_t:
            r = operating_point(p, y, w, t)
            r.update(model=name, threshold=t, nb_w=net_benefit(p, y, w, t))
            rows.append(r)
    pd.DataFrame(rows)[["model", "threshold", "sens", "spec", "ppv", "referred",
                        "wsens", "wspec", "wreferred", "nb_w"]] \
        .to_csv(os.path.join(OUT, "dca_operating_points.csv"), index=False)

    # ---- fixed referral burden (weighted) and the burden needed for 90% sensitivity
    W = w.sum()
    rows = []
    for name, p in models.items():
        order = np.argsort(-p)
        for frac in (0.25, 0.40):
            k = int(np.searchsorted(np.cumsum(w[order]) / W, frac) + 1)
            sel = order[:k]
            comp = np.setdiff1d(np.arange(len(p)), sel)
            t = p[order][k - 1]
            rows.append({
                "model": name, "scenario": f"top {int(frac*100)}% (weighted) referred",
                "threshold": round(float(t), 3),
                "nb_w": round(net_benefit(p, y, w, t), 4),
                "referred": round(w[sel].sum() / W, 4),
                "sens": round(w[sel][y[sel] == 1].sum() / w[y == 1].sum(), 4),
                "spec": round(w[comp][y[comp] == 0].sum() / w[y == 0].sum(), 4),
                "k": k,
            })
        cum = np.cumsum(w[order] * (y[order] == 1)) / w[y == 1].sum()
        k = int(np.searchsorted(cum, 0.90) + 1)
        sel = order[:k]
        comp = np.setdiff1d(np.arange(len(p)), sel)
        t = p[order][k - 1]
        rows.append({
            "model": name, "scenario": "referral needed for 90% sensitivity",
            "threshold": round(float(t), 3),
            "nb_w": round(net_benefit(p, y, w, t), 4),
            "referred": round(w[sel].sum() / W, 4),
            "sens": round(w[sel][y[sel] == 1].sum() / w[y == 1].sum(), 4),
            "spec": round(w[comp][y[comp] == 0].sum() / w[y == 0].sum(), 4),
            "k": k,
        })
    pd.DataFrame(rows)[["model", "scenario", "threshold", "nb_w", "referred",
                        "sens", "spec", "k"]] \
        .to_csv(os.path.join(OUT, "dca_fixed_operating.csv"), index=False)

    print("wrote dca_net_benefit_curve.csv, dca_operating_points.csv, dca_fixed_operating.csv")

    if "--figure" in sys.argv:
        make_figure(curve, models, y, w, sys.argv[sys.argv.index("--figure") + 1])


def make_figure(curve, models, y, w, path):
    """Two-panel Figure S2. Not committed as a binary (the disclosure gate blocks new
    binaries); regenerate with ``--figure figures/FigureS2_decision_analysis.png``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    col = {"full_reference": "#1F6FB2", "local_lipid": "#E1962B", "transported": "#B4639B"}
    lab = {"full_reference": "Full reference model",
           "local_lipid": "Local lipid-only model",
           "transported": "Transported Chinese score"}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.linewidth": 0.8, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"wspace": 0.26})

    a = ax[0]
    for m in col:
        s = curve[curve.model == m]
        a.plot(s.threshold, s.nb_w, color=col[m], lw=1.6, label=lab[m])
    s = curve[curve.model == "full_reference"]
    a.plot(s.threshold, s.nb_all_w, color="0.45", lw=1.2, ls="--", label="Refer all")
    a.plot(s.threshold, np.zeros(len(s)), color="0.15", lw=0.9, ls=":", label="Refer none")
    a.axvline(w[y == 1].sum() / w.sum(), color="0.82", lw=0.8)
    a.set_xlim(0.05, 0.60)
    a.set_ylim(-0.10, 0.32)
    a.set_xlabel("Threshold probability")
    a.set_ylabel("Net benefit")
    a.set_title("A  Decision curves (survey-weighted)", loc="left", fontsize=8.5)
    a.text(0.40, -0.075, "vertical line: prevalence", fontsize=6, color="0.55")
    a.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.7,
             borderpad=0.1, labelspacing=0.35)

    b = ax[1]
    for m, p in models.items():
        order = np.argsort(-p)
        cumy = np.cumsum(w[order] * (y[order] == 1)) / w[y == 1].sum()
        cumw = np.cumsum(w[order]) / w.sum()
        b.plot(cumw, cumy, color=col[m], lw=1.6, label=lab[m])
    b.plot([0, 1], [0, 1], color="0.55", lw=1.0, ls="--", label="No ranking")
    b.set_xlim(0, 1)
    b.set_ylim(0, 1)
    b.set_xlabel("Proportion referred")
    b.set_ylabel("Proportion of cases detected")
    b.set_title("B  Detection yield at a given referral burden", loc="left", fontsize=8.5)
    for f in (0.25, 0.40):
        b.axvline(f, color="0.85", lw=0.8, zorder=0)
    b.legend(frameon=False, fontsize=6.5, loc="lower right", handlelength=1.7,
             borderpad=0.1, labelspacing=0.35)

    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print("wrote figure:", path)


if __name__ == "__main__":
    main()
