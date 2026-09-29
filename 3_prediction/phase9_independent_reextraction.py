# -*- coding: utf-8 -*-
"""phase9: independent re-extraction of the DR outcome from the source ophthalmology text (Route D1)

Structural difference from the primary pipeline (rebuild_outcome.py):
  primary  : a first-match keyword cascade over the whole eye text, with a fixed priority order
             (ungradable -> non-proliferative grades -> bare proliferative -> explicit absence).
  this one : (1) every eye text is normalised (病变/改变 unified, full-width punctuation unified);
             (2) the text is split into clauses;
             (3) EVERY clause is scored independently for (negation scope, severity term, target term);
             (4) the eye grade is the maximum over contributing clauses, and any logical conflict
                 between clauses is reported rather than silently resolved by order;
             (5) the record grade is the worst gradable eye.
No clause exits early, so the substring-ordering failure mode that produced the original coding
error cannot occur by construction.

Outputs
  output_tables/layer2_independent_reextraction.csv   (per record: both parsers, agreement flags)
  output_tables/layer2_independent_reextraction_summary.csv
"""
import sys, os, re, json
import numpy as np
import pandas as pd

SRC = r"E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output_tables")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- normalisation
NEG = r'(?:未见|未|无|没有|排除|不考虑|不存在)'
QUANT_RE = re.compile(r'出血数量|出血总面积|出血最大面积|渗出数量|渗出总面积|渗出最大面积')
# NOTE: the proliferative pattern must not match the 增殖性 inside 非增殖性 — that is precisely the
# substring trap that produced the original coding error. A negative lookbehind makes it explicit
# rather than relying on the order in which patterns happen to be tested.
SEV = [('proliferative DR', 4, r'(?<!非)增殖性'),
       ('severe NPDR', 3, r'重度\s*非增殖性'),
       ('moderate NPDR', 2, r'中度\s*非增殖性'),
       ('mild NPDR', 1, r'轻度\s*非增殖性')]
BARE_GRADE = [('severe', 3, r'重度'), ('moderate', 2, r'中度'), ('mild', 1, r'轻度')]
# adjudication switch: True reproduces the pre-adjudication variant that mapped 可疑 DR to grade 0
SUSPECT_AS_GRADE0 = bool(int(os.environ.get('SUSPECT_AS_GRADE0', '0')))
TARGET = r'(?:糖尿病)?视网膜(?:病变|改变)|糖网'
UNGRADABLE = r'图片不可读|不可读|眼底看不清|看不清|屈光介质混浊|模糊|无法判读|图像不清'

def normalise(s):
    s = str(s)
    s = s.replace('，', ',').replace('；', ';').replace('：', ':').replace('。', '.')
    s = s.replace('改变', '病变')          # 病变 / 改变 are used interchangeably in this EMR
    return s

def clauses(s):
    return [c.strip() for c in re.split(r'[,;.\n]+', normalise(s)) if c.strip()]

def score_eye(raw):
    """returns (grade or None, evidence list, flags)"""
    txt = normalise(raw)
    if not txt.strip():
        return None, ['empty text'], ['ungradable:empty']
    cl = clauses(raw)
    evidence, flags, grades = [], [], []

    for c in cl:
        # quantitative lesion-count clause, carries no diagnostic grade
        if QUANT_RE.search(c):
            evidence.append('quant-only:%s' % c[:28]); continue
        # explicit unreadable / ungradable statement
        if re.search(UNGRADABLE, c):
            flags.append('ungradable-text'); evidence.append('ungradable:%s' % c[:28]); continue
        # target condition present?
        if not re.search(TARGET, c):
            continue
        negated = bool(re.search(NEG + r'.{0,12}' + TARGET, c)) or bool(re.search(NEG, c))
        qualifier = 'suspected' if re.search(r'可疑|疑似|考虑', c) else ('visible' if re.search(r'可见|见', c) else 'plain')
        g = None
        for name, val, pat in SEV:
            if re.search(pat, c):
                g = val; break
        if g is None:
            for name, val, pat in BARE_GRADE:          # bare 轻/中/重度 without 非增殖性
                if re.search(pat, c):
                    g = val; flags.append('bare-grade-term'); break
        if g is None:
            if negated and qualifier in ('plain', 'visible'):
                grades.append((0, 'explicit-absence:%s' % c[:24]))
            elif qualifier == 'suspected':
                # ADJUDICATION (see Table S14 note): '可疑' (suspected) DR permits no ICDR grade and
                # is NOT equivalent to "no apparent retinopathy", so under the prespecified rule
                # ("an explicit ICDR grade is required") it is indeterminate. An earlier run of this
                # same parser mapped it to 0; that variant is reported as the pre-adjudication stage
                # because it is what surfaced the ambiguity.
                if SUSPECT_AS_GRADE0:
                    grades.append((0, 'suspected->0:%s' % c[:24])); flags.append('suspected')
                else:
                    evidence.append('suspected-no-grade:%s' % c[:24])
                    flags.append('suspected->indeterminate')
            else:
                # DR stated but no severity grade -> not gradeable
                evidence.append('no-grade-stated:%s' % c[:28]); flags.append('dr-without-grade')
            continue
        if negated:
            grades.append((0, 'negated-severity:%s' % c[:24]))
        else:
            grades.append((g, '%s:%s' % (name, c[:24])))

    if not grades:
        return None, evidence, flags + ['ungradable:no-grade-clause']
    vals = sorted({v for v, _ in grades})
    if len(vals) > 1:
        flags.append('conflict')
        evidence.append('CONFLICT clauses=%s' % vals)
    return max(vals), evidence, flags

# ---------------------------------------------------------------- run over the cohort
aud = pd.read_csv(f"{SRC}/severity_rebuild_audit.csv", low_memory=False)
dev = pd.read_csv(f"{SRC}/dev_set_v2.csv", low_memory=False)
assert len(aud) == len(dev) == 2482 and (aud['pid'].values == dev['id'].values).all()

rows = []
for i, r in aud.iterrows():
    gl, el, fl = score_eye(r['Ltxt']); gr, er, fr = score_eye(r['Rtxt'])
    flagged = fl + fr
    eyes = [g for g in (gl, gr) if g is not None]
    if 'ungradable-text' in flagged and eyes:
        gradable = 1; grade = max(eyes); flags = ';'.join(sorted(set(flagged)) + ['mixed:ungradable-and-grade'])
    elif eyes and 'dr-without-grade' in flagged:
        gradable = 1; grade = max(eyes); flags = ';'.join(sorted(set(flagged)))
    elif eyes:
        gradable = 1; grade = max(eyes); flags = ';'.join(sorted(set(flagged)))
    else:
        gradable = 0; grade = np.nan; flags = ';'.join(sorted(set(flagged)))
    rows.append(dict(pid=r['pid'], sev_new=dev.loc[i, 'severity'], gradeable_new=dev.loc[i, 'gradeable'],
                     any_dr_new=dev.loc[i, 'any_dr'],
                     sev_ind=grade, gradeable_ind=gradable, any_dr_ind=(1 if (gradable and grade >= 1) else 0),
                     flags=str(flags), evidence=' | '.join((el + er)[:4])))
R = pd.DataFrame(rows)
R['grade_match'] = (R.sev_new.fillna(-1) == R.sev_ind.fillna(-1)).astype(int)
R['anydr_match'] = (R.any_dr_new.astype(int) == R.any_dr_ind.astype(int)).astype(int)
R['gradable_match'] = (R.gradeable_new.astype(int) == R.gradeable_ind.astype(int)).astype(int)
R.to_csv(f"{OUT}/layer2_independent_reextraction.csv", index=False)

from sklearn.metrics import cohen_kappa_score
g = R[(R.gradeable_new == 1) & (R.gradeable_ind == 1)]
both = R[R.gradeable_new.eq(1) & R.gradeable_ind.eq(1)]
kw = cohen_kappa_score(both.sev_new.astype(int), both.sev_ind.astype(int), weights='quadratic')
ku = cohen_kappa_score(R.any_dr_new.astype(int), R.any_dr_ind.astype(int))
kg = cohen_kappa_score(R.gradeable_new.astype(int), R.gradeable_ind.astype(int))

summary = {
 'records': int(len(R)),
 'gradable_agree': int(R.gradable_match.sum()), 'gradable_total': int(len(R)),
 'anydr_agree': int(R.anydr_match.sum()),
 'anydr_n': int(len(R)),
 'five_level_agree_both_gradable': int(R.loc[R.gradeable_new.eq(1) & R.gradeable_ind.eq(1), 'grade_match'].sum()),
 'five_level_n': int((R.gradeable_new.eq(1) & R.gradeable_ind.eq(1)).sum()),
 'kappa_gradable': round(float(kg), 3),
 'kappa_any_dr': round(float(ku), 3),
 'kappa_five_level_quadratic_weighted': round(float(kw), 3),
 'discordant_records': int((R.anydr_match == 0).sum() + (R.gradable_match == 0).sum()),
 'suspect_as_grade0_variant': SUSPECT_AS_GRADE0,
}

# ---- bootstrap CIs for the three kappa statistics ----
rng = np.random.default_rng(20260929)
bs = {'gradable': [], 'any_dr': [], 'five_level': []}
idx = np.arange(len(R))
for _ in range(2000):
    s = rng.integers(0, len(R), len(R))
    b = R.iloc[s]
    bs['gradable'].append(cohen_kappa_score(b.gradeable_new.astype(int), b.gradeable_ind.astype(int)))
    bs['any_dr'].append(cohen_kappa_score(b.any_dr_new.astype(int), b.any_dr_ind.astype(int)))
    bb = b[b.gradeable_new.eq(1) & b.gradeable_ind.eq(1)]
    if bb.sev_new.notna().all() and bb.sev_ind.notna().all() and bb.sev_new.nunique() > 1:
        bs['five_level'].append(cohen_kappa_score(bb.sev_new.astype(int), bb.sev_ind.astype(int), weights='quadratic'))
for k, v in bs.items():
    if v:
        summary['kappa_%s_lo' % k] = round(float(np.nanpercentile(v, 2.5)), 3)
        summary['kappa_%s_hi' % k] = round(float(np.nanpercentile(v, 97.5)), 3)

# ---- the discordant records, listed for adjudication ----
bad = R[(R.gradable_match == 0) | (R.anydr_match == 0) | (R.grade_match == 0)].copy()
src_aud = aud.set_index('pid')
bad['Ltxt'] = [str(src_aud.loc[p, 'Ltxt'])[:200] for p in bad.pid]
bad['Rtxt'] = [str(src_aud.loc[p, 'Rtxt'])[:200] for p in bad.pid]
bad.to_csv(f"{OUT}/layer2_independent_reextraction_discordances.csv", index=False)

# ---- sensitivity: if the six suspected records are counted as non-cases ----
if len(bad):
    g2 = dev[dev.gradeable == 1]
    n2 = len(g2) + len(bad); e2 = int(g2.any_dr.sum())
    summary['sensitivity_n'] = n2
    summary['sensitivity_prevalence_pct'] = round(100 * e2 / n2, 2)
    summary['primary_prevalence_pct'] = round(100 * e2 / len(g2), 2)

pd.DataFrame([summary]).to_csv(f"{OUT}/layer2_independent_reextraction_summary.csv", index=False)
print(json.dumps(summary, ensure_ascii=False, indent=2))
print("\n=== 五级分布对照 ===")
print(pd.crosstab(R.sev_new, R.sev_ind, dropna=False))
print("\n=== any-DR 交叉表 ===")
print(pd.crosstab(R.any_dr_new, R.any_dr_ind))
print("\n=== 不一致记录（前 40 条）===")
bad = R[(R.gradable_match == 0) | (R.anydr_match == 0) | (R.grade_match == 0)]
for _, r in bad.head(40).iterrows():
    a = aud[aud.pid == r.pid].iloc[0]
    print("pid %s | new %s/%s | ind %s/%s | flags=%s\n    L=%s\n    R=%s"
          % (r.pid, r.sev_new, r.any_dr_new, r.sev_ind, r.any_dr_ind, r.flags,
             str(a['Ltxt'])[:110], str(a['Rtxt'])[:110]))
print("\n不一致记录总数:", len(bad))
