# -*- coding: utf-8 -*-
"""phase8: peer-review revision analyses for layer 2
   1) Parsimonious lipid models (age + sex + one marker) with optimism-corrected dAUC
   2) Lipid missing-data pattern distribution
   3) Record-level cross-tabulation of the earlier outcome coding vs the rebuilt outcome
   Inputs :  <SRC>/phase2_dev/dev_set.csv (earlier outcome), dev_set_v2.csv (rebuilt outcome)
   Outputs:  output_tables/layer2_parsimonious_models.csv
             output_tables/layer2_missingness_patterns.csv
             output_tables/layer2_outcome_crosstab_full.csv
   Note   :  dev_set_v2.csv is row-aligned with dev_set.csv (2,482 screening records).
"""
import sys, os, numpy as np, pandas as pd, statsmodels.api as sm

SRC = r"E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output_tables")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20260929)
LIP = ['tg', 'tc', 'hdl_c', 'ldl_c', 'apo_a1', 'apo_b']

dev2 = pd.read_csv(f"{SRC}/dev_set_v2.csv", low_memory=False)
dev1 = pd.read_csv(f"{SRC}/dev_set.csv", low_memory=False)
assert (dev1['id'].values == dev2['id'].values).all(), "dev sets are no longer row-aligned"
g = dev2[dev2['gradeable'] == 1].copy()
g['log2tg'] = np.log2(g.tg)
g['log2tg_hdl'] = np.log2(g.tg / g.hdl_c)
g['non_hdl'] = g.tc - g.hdl_c
g['rc'] = g.tc - g.ldl_c - g.hdl_c
g['apo_ratio'] = g.apo_b / g.apo_a1

# ---------------- 1) parsimonious models ----------------
def auc(y, s):
    y = np.asarray(y, float); r = pd.Series(np.asarray(s, float)).rank().values
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1)

def fit_auc(X, y, idx_fit, idx_eval):
    try:
        m = sm.Logit(y.iloc[idx_fit], sm.add_constant(X.iloc[idx_fit], has_constant='add')).fit(disp=0, maxiter=100)
    except Exception:
        return np.nan
    p = m.predict(sm.add_constant(X.iloc[idx_eval], has_constant='add'))
    yy = y.iloc[idx_eval].values
    if yy.sum() in (0, len(yy)): return np.nan
    return auc(yy, p.values)

MARKERS = [('TG (log2)', 'log2tg'), ('Total cholesterol', 'tc'), ('HDL-C', 'hdl_c'), ('LDL-C', 'ldl_c'),
           ('ApoA-I', 'apo_a1'), ('ApoB', 'apo_b'), ('Non-HDL-C', 'non_hdl'),
           ('Remnant cholesterol', 'rc'), ('TG/HDL-C (log2)', 'log2tg_hdl'), ('ApoB/ApoA-I', 'apo_ratio')]
B = 400
rows = []
for label, col in MARKERS:
    sub = g[['age', 'male', col, 'any_dr']].dropna().reset_index(drop=True)
    Xa, Xb, y = sub[['age', 'male']], sub[['age', 'male', col]], sub['any_dr'].astype(float)
    n = len(sub); idx = np.arange(n)
    a_app, b_app = fit_auc(Xa, y, idx, idx), fit_auc(Xb, y, idx, idx)
    boot = []
    for _ in range(B):
        s = rng.integers(0, n, n)
        if y.iloc[s].sum() in (0, n): continue
        aa, bb = fit_auc(Xa, y, s, s), fit_auc(Xb, y, s, s)
        ao, bo = fit_auc(Xa, y, s, idx), fit_auc(Xb, y, s, idx)
        boot.append((aa, bb, bb - aa, ao, bo))
    boot = np.array(boot)
    lo, hi = np.nanpercentile(boot[:, 2], [2.5, 97.5])
    opt_a = np.nanmean(boot[:, 0] - boot[:, 3]); opt_b = np.nanmean(boot[:, 1] - boot[:, 4])
    rows.append(dict(marker=label, n=n, events=int(y.sum()), auc_age_sex=round(a_app, 3),
                     auc_plus_marker=round(b_app, 3), dAUC=round(b_app - a_app, 3),
                     dAUC_lo=round(lo, 3), dAUC_hi=round(hi, 3),
                     dAUC_optimism_corrected=round((b_app - opt_b) - (a_app - opt_a), 3)))
    print("done", label, flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/layer2_parsimonious_models.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))

# ---------------- 2) missingness patterns ----------------
nm = g[LIP].isna().sum(axis=1)
pat = g[LIP].isna().apply(lambda r: "".join('1' if x else '0' for x in r), axis=1)
dist = nm.value_counts().sort_index().rename_axis('lipids_missing').reset_index(name='records')
dist['pct'] = (100 * dist.records / len(g)).round(1)
dist.to_csv(f"{OUT}/layer2_missingness_patterns.csv", index=False)
print("\npatterns:", dist.to_dict('records'))
print("any missing:", int((nm > 0).sum()), "mean missing among those: %.2f" % nm[nm > 0].mean())
print(pat.value_counts().head(8).to_dict())

# ---------------- 3) outcome cross-tabulation ----------------
m = pd.DataFrame({
    'sev_old': dev1['severity'].values, 'grade_old': dev1['gradeable'].values, 'dr_old': dev1['any_dr'].values,
    'sev_new': dev2['severity'].values, 'grade_new': dev2['gradeable'].values, 'dr_new': dev2['any_dr'].values})
m['new_class'] = np.where(m.grade_new == 0, 'not_gradable', np.where(m.dr_new == 1, 'DR', 'no_DR'))
m['old_class'] = np.where(m.grade_old == 1, np.where(m.dr_old == 1, 'analysed_as_case', 'analysed_as_non_case'),
                          'excluded_' + m.sev_old.astype(str))
ct = pd.crosstab(m.old_class, m.new_class, margins=True)
ct.to_csv(f"{OUT}/layer2_outcome_crosstab_full.csv")
print("\n", ct.to_string())
