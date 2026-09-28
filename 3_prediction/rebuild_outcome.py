# -*- coding: utf-8 -*-
"""Rebuild dev_set with a CORRECTED DR outcome, following the original clean_data.py
   parsing but with the substring trap fixed (非增殖性 must be tested before 增殖性).

   Per eye priority:
     1  ungradable text (屈光介质混浊 / 图片不可读 / 眼底看不清)      -> ungradable
     2  重度非增殖性 ->3 ; 中度非增殖性 ->2 ; 轻度非增殖性 ->1
        bare 增殖性糖尿病视网膜 (NOT preceded by 非) ->4
     3  未见明显糖尿病视网膜病变 / 无糖尿病视网膜病变 / 未见明显糖尿病视网膜改变 ->0
     4  可见糖尿病视网膜病变 (DR explicitly present, grade not stated)  -> feature-based
     5  可疑糖尿病视网膜病变                                          -> suspect
     6  出血数量… (pure AI quantitation, no diagnosis)                -> quant
     7  otherwise                                                     -> other

   Patient level (primary, "option A"): grade = max over eyes with an explicit severity
   statement (categories 2-3 above); DR = grade >= 1; patients with no such eye are excluded.
   Sensitivity variant A' additionally admits category 4 (grade inferred from features).
"""
import sys, os, re, glob, collections
import pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')

BASE = r'E:/workbuddy/血脂，糖网等数据清洗'
OUTDIR = r'E:/workbuddy/2026-09-04-18-19-16/phase2_dev'
FEAT = ["微血管瘤", "出血", "硬性渗出", "棉絮斑", "静脉串珠", "微血管异常", "新生血管",
        "玻璃体积血", "视网膜前出血", "纤维增殖膜", "视网膜脱离", "黄斑水肿"]


def parse_eye(text):
    t = '' if pd.isna(text) else str(text).strip()
    if t == '':
        return dict(cat='missing', grade=np.nan)
    if ('屈光介质混浊' in t) or ('图片不可读' in t) or ('眼底看不清' in t):
        return dict(cat='ungradable', grade=np.nan)
    # ---- explicit severity text: 非增殖性 MUST be tested before bare 增殖性 ----
    if '重度非增殖性' in t:
        return dict(cat='sev3', grade=3)
    if '中度非增殖性' in t:
        return dict(cat='sev2', grade=2)
    if '轻度非增殖性' in t:
        return dict(cat='sev1', grade=1)
    if re.search(r'(?<!非)增殖性糖尿病视网膜', t):
        return dict(cat='sev4_pdr', grade=4)
    if ('未见明显糖尿病视网膜病变' in t) or ('无糖尿病视网膜病变' in t) or ('未见明显糖尿病视网膜改变' in t):
        return dict(cat='nodr_text', grade=0)
    if '可见糖尿病视网膜病变' in t:
        return dict(cat='visible_dr', grade=np.nan)
    if '可疑糖尿病视网膜病变' in t:
        return dict(cat='suspect', grade=np.nan)
    if '出血数量' in t:
        return dict(cat='quant', grade=np.nan)
    return dict(cat='other', grade=np.nan)


def feature_grade(L, R):
    """grade from structured feature flags, used only for the A' sensitivity."""
    g = []
    for feats in (L, R):
        if any(feats.get(x) for x in ['f_新生血管', 'f_玻璃体积血', 'f_视网膜前出血',
                                      'f_纤维增殖膜', 'f_视网膜脱离']):
            g.append(4)
        elif feats.get('f_黄斑水肿'):
            g.append(3)
        elif any(feats.get(x) for x in ['f_棉絮斑', 'f_静脉串珠', 'f_微血管异常']):
            g.append(2)
        else:
            g.append(1)
    return max(g)


# ---------- read raw screening files ----------
files = sorted(glob.glob(os.path.join(BASE, '2019-2023糖网数据', '*.xls')))
rows = []
for f in files:
    raw = pd.read_excel(f, header=None)
    hdr = raw[raw[0] == '序号'].index[0]
    d = pd.read_excel(f, header=hdr)
    d.columns = [str(c).strip() for c in d.columns]
    yr = os.path.basename(f)[:4]
    for _, r in d.iterrows():
        if pd.isna(r.get('患者编号')) and pd.isna(r.get('患者姓名')):
            continue
        Ltxt = r.get('左眼糖网筛查评估结果（仅供参考）')
        Rtxt = r.get('右眼糖网筛查评估结果（仅供参考）')
        L = parse_eye(Ltxt)
        R = parse_eye(Rtxt)
        Lf = {f'f_{x}': (f'可见{x}' in str(Ltxt)) or (f'可疑{x}' in str(Ltxt)) for x in FEAT}
        Rf = {f'f_{x}': (f'可见{x}' in str(Rtxt)) or (f'可疑{x}' in str(Rtxt)) for x in FEAT}
        rows.append(dict(
            year=yr, pid=str(r.get('患者编号')).strip(), hosp=str(r.get('医院编号')),
            name=r.get('患者姓名'), sex=r.get('性别'), age=r.get('年龄'),
            Lcat=L['cat'], Lgrade=L['grade'], Rcat=R['cat'], Rgrade=R['grade'],
            Ltxt=str(Ltxt), Rtxt=str(Rtxt), **{k: v for k, v in Lf.items()},
            **{k + '_R': v for k, v in Rf.items()},
            L_hem_num=r.get('左眼出血数') if '左眼出血数' in d.columns else np.nan,
        ))
scr = pd.DataFrame(rows)
print('原始筛查记录:', len(scr), ' 唯一患者编号:', scr.pid.nunique())

# ---------- patient-level ----------
def agg(r):
    cats = [r.Lcat, r.Rcat]
    grades = [g for g in (r.Lgrade, r.Rgrade) if pd.notna(g)]
    sev = max(grades) if grades else np.nan
    return pd.Series(dict(
        sev_text=sev,
        graded=bool(grades),
        any_dr=max(grades) >= 1 if grades else False,
        cat_pair='|'.join(sorted(cats)),
    ))

scr = pd.concat([scr, scr.apply(agg, axis=1)], axis=1)
print()
print('=== 患者层面：逐眼类别组合 Top 12 ===')
for k, v in scr.cat_pair.value_counts().head(12).items():
    print('   %-26s %4d' % (k, v))
print()
print('=== 方案 A（仅显式分级文字）===')
print('  可评级（graded）            :', int(scr.graded.sum()))
print('    DR (grade>=1)            :', int(scr.any_dr.sum()))
print('    其中 PDR (grade 4)        :', int((scr.sev_text == 4).sum()))
print('    无DR (grade 0)            :', int((scr.sev_text == 0).sum()))
print('  不可评级（需排除）           :', int((~scr.graded).sum()))
print()
print('  -- 不可评级者的逐眼类别构成 --')
print(scr[~scr.graded].cat_pair.value_counts().head(10).to_string())
print()
prev = 100 * scr.any_dr.mean()
print('  患病率（可评级者中）        : %.1f%%' % prev)

scr.to_csv(os.path.join(OUTDIR, 'severity_rebuild_audit.csv'), index=False, encoding='utf-8')
print('\n已写出 severity_rebuild_audit.csv')
