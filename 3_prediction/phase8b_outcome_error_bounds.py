# -*- coding: utf-8 -*-
"""phase8b: how much could unverified outcome coding change the layer-2 conclusion?
   1) Internal-consistency audit of the reconstructed outcome against the source record
   2) Probabilistic bias analysis over assumed sensitivity / specificity of the outcome
   3) Extreme bounds: assigning every not-gradable record to non-cases, then to cases
   Inputs : <SRC>/phase2_dev/dev_set_v2.csv, severity_rebuild_audit.csv
   Outputs: output_tables/layer2_outcome_audit_and_bounds.csv
"""
import sys, os, numpy as np, pandas as pd, statsmodels.api as sm

SRC = r"E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output_tables")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20260929)

dev = pd.read_csv(f"{SRC}/dev_set_v2.csv", low_memory=False)
aud = pd.read_csv(f"{SRC}/severity_rebuild_audit.csv", low_memory=False)
assert len(dev) == len(aud) == 2482, "row alignment lost"
txt = (aud['Ltxt'].fillna('') + ' ' + aud['Rtxt'].fillna(''))
les = ['L_hem_num', 'L_exu_num', 'R_hem_num', 'R_exu_num']
flags = ['L_micro', 'L_hardex', 'L_macular', 'L_neov', 'R_micro', 'R_hardex', 'R_macular', 'R_neov']
les_pos = (dev[les].fillna(0) > 0).any(axis=1) | (dev[flags].fillna(0) > 0).any(axis=1)
dr_kw = txt.str.contains('非增殖性|增殖性|视网膜病变|微血管瘤', regex=True).values
nodr_kw = txt.str.contains('未见明显糖尿病视网膜病变|未见微血管瘤|未见出血', regex=True).values

rows = []
rows.append(dict(block='A_consistency', item='graded no-DR but positive lesion count/flag',
                 n=int((dev.gradeable.eq(1) & dev.any_dr.eq(0) & les_pos).sum()), denom=int((dev.gradeable.eq(1) & dev.any_dr.eq(0)).sum())))
rows.append(dict(block='A_consistency', item='rebuilt as DR but no DR keyword',
                 n=int((dev.gradeable.eq(1) & dev.any_dr.eq(1) & ~dr_kw).sum()), denom=int((dev.gradeable.eq(1) & dev.any_dr.eq(1)).sum())))
rows.append(dict(block='A_consistency', item='explicit absent-DR text but rebuilt as DR',
                 n=int((dev.any_dr.eq(1) & nodr_kw & ~dr_kw).sum()), denom=int(dev.any_dr.sum())))
rows.append(dict(block='A_consistency', item='not-gradable records carrying DR keywords',
                 n=int((dev.gradeable.eq(0) & dr_kw).sum()), denom=int((dev.gradeable.eq(0)).sum())))

# ---- 2) probabilistic bias analysis ----
g = dev[dev.gradeable == 1].copy()
g['log2tg'] = np.log2(g.tg)
CC = g[['age', 'male', 'log2tg', 'tc', 'hdl_c', 'ldl_c', 'apo_a1', 'apo_b', 'any_dr']].dropna().reset_index(drop=True)
X = CC[['age', 'male', 'log2tg', 'tc', 'hdl_c', 'ldl_c', 'apo_a1', 'apo_b']]
d_obs = CC['any_dr'].astype(int).values
p_obs = d_obs.mean()

def auc(y, s):
    y = np.asarray(y, float); r = pd.Series(np.asarray(s, float)).rank().values
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n0 * n1)

def fit_auc(yv):
    try:
        m = sm.Logit(yv, sm.add_constant(X, has_constant='add')).fit(disp=0, maxiter=60)
    except Exception:
        return np.nan
    return auc(yv, m.predict(sm.add_constant(X, has_constant='add')).values)

print("observed AUC (complete cases) = %.3f" % fit_auc(d_obs))
B = 100
for Se in [1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.40]:
    for Sp in [0.99, 0.95]:
        P1 = p_obs * Se + (1 - p_obs) * (1 - Sp)
        PPV = p_obs * Se / P1
        NPV = (1 - p_obs) * Sp / (1 - P1)
        prob = np.where(d_obs == 1, PPV, 1 - NPV)
        aucs = []
        for _ in range(B):
            y = (rng.random(len(d_obs)) < prob).astype(float)
            if y.sum() in (0, len(y)): continue
            aucs.append(fit_auc(y))
        rows.append(dict(block='B_bias_analysis', item='Se=%.2f, Sp=%.2f' % (Se, Sp),
                         n=len(aucs), auc_mean=round(np.nanmean(aucs), 3),
                         auc_lo=round(np.nanpercentile(aucs, 2.5), 3),
                         auc_hi=round(np.nanpercentile(aucs, 97.5), 3),
                         true_prev=round(p_obs * Se / P1 * p_obs + (1 - NPV * 0) * 0, 4)))
        print("  Se=%.2f Sp=%.2f -> AUC %.3f" % (Se, Sp, np.nanmean(aucs)), flush=True)

# ---- 3) extreme bounds on the 511 not-gradable records ----
exc = dev[dev.gradeable == 0]
for lbl, val in [('all not-gradable as non-cases', 0), ('all not-gradable as cases', 1),
                 ('primary analysis (excluded)', None)]:
    if val is None:
        d = g[['tg', 'tc', 'hdl_c', 'ldl_c', 'any_dr']].copy()
    else:
        d = pd.concat([g[['tg', 'tc', 'hdl_c', 'ldl_c', 'any_dr']],
                       pd.DataFrame({'tg': exc.tg, 'tc': exc.tc, 'hdl_c': exc.hdl_c, 'ldl_c': exc.ldl_c, 'any_dr': val})],
                      ignore_index=True)
    res = {}
    for c in ['tg', 'tc', 'hdl_c', 'ldl_c']:
        s = d[[c, 'any_dr']].dropna()
        res[c] = round(auc(s.any_dr, s[c]), 3)
    rows.append(dict(block='C_extreme_bounds', item=lbl, n=len(d), prevalence=round(100 * d.any_dr.mean(), 1), **res))
    print("  %-32s prevalence %.1f%%  %s" % (lbl, 100 * d.any_dr.mean(), res), flush=True)

# ---- 4) selective-exclusion check ----
inc = g[g.any_dr == 1]
print("\nmean lipids - included DR (n=%d): TG %.2f TC %.2f | not-gradable with DR keywords (n=%d): TG %.2f TC %.2f"
      % (len(inc), inc.tg.mean(), inc.tc.mean(), int((dev.gradeable.eq(0) & dr_kw).sum()),
         dev[(dev.gradeable.eq(0)) & dr_kw].tg.mean(), dev[(dev.gradeable.eq(0)) & dr_kw].tc.mean()))

pd.DataFrame(rows).to_csv(f"{OUT}/layer2_outcome_audit_and_bounds.csv", index=False)
print("\nwritten:", f"{OUT}/layer2_outcome_audit_and_bounds.csv")
