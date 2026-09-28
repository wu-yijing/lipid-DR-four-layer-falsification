# -*- coding: utf-8 -*-
"""Write dev_set_v2.csv (corrected DR outcome) with the SAME schema/order as dev_set.csv,
   so every downstream R script runs unchanged against a corrected input.

   Primary  (option A) : DR = explicit severity text (轻/中/重 NPDR, bare PDR) on either eye
   Sensitivity (option A'): additionally admits '可见糖尿病视网膜病变' (DR stated, grade from
                          features) and pure AI-quantitation records (grade from haem/ exudate
                          counts, the convention already used in clean_data.py line 128)
"""
import sys, os, re, glob
import pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')

BASE = r'E:/workbuddy/血脂，糖网等数据清洗'
OUT = r'E:/workbuddy/2026-09-04-18-19-16/phase2_dev'
FEAT = ["微血管瘤", "出血", "硬性渗出", "棉絮斑", "静脉串珠", "微血管异常", "新生血管",
        "玻璃体积血", "视网膜前出血", "纤维增殖膜", "视网膜脱离", "黄斑水肿"]


def parse_eye(text):
    t = '' if pd.isna(text) else str(text).strip()
    if t == '':
        return ('missing', np.nan)
    if ('屈光介质混浊' in t) or ('图片不可读' in t) or ('眼底看不清' in t):
        return ('ungradable', np.nan)
    if '重度非增殖性' in t:
        return ('sev', 3)
    if '中度非增殖性' in t:
        return ('sev', 2)
    if '轻度非增殖性' in t:
        return ('sev', 1)
    if re.search(r'(?<!非)增殖性糖尿病视网膜', t):
        return ('sev', 4)
    if ('未见明显糖尿病视网膜病变' in t) or ('无糖尿病视网膜病变' in t) or ('未见明显糖尿病视网膜改变' in t):
        return ('nodr', 0)
    if '可见糖尿病视网膜病变' in t:
        return ('visible', np.nan)
    if '可疑糖尿病视网膜病变' in t:
        return ('suspect', np.nan)
    if '出血数量' in t:
        return ('quant', np.nan)
    return ('other', np.nan)


def quant_grade(t):
    m = re.search(r'出血数量：(\d+).*?渗出数量：(\d+)', str(t))
    if not m:
        return np.nan
    return 1 if (int(m.group(1)) > 0 or int(m.group(2)) > 0) else 0


def fgrade(t):
    g = {}
    for x in FEAT:
        g[x] = (f'可见{x}' in str(t)) or (f'可疑{x}' in str(t))
    if any(g[x] for x in ['新生血管', '玻璃体积血', '视网膜前出血', '纤维增殖膜', '视网膜脱离']):
        return 4
    if g['黄斑水肿']:
        return 3
    if g['棉絮斑'] or g['静脉串珠'] or g['微血管异常']:
        return 2
    return 1


# ---------- read raw screening files (same order/skip rule as clean_data.py) ----------
files = sorted(glob.glob(os.path.join(BASE, '2019-2023糖网数据', '*.xls')))
recs = []
for f in files:
    raw = pd.read_excel(f, header=None)
    hdr = raw[raw[0] == '序号'].index[0]
    d = pd.read_excel(f, header=hdr)
    d.columns = [str(c).strip() for c in d.columns]
    for _, r in d.iterrows():
        if pd.isna(r.get('患者编号')) and pd.isna(r.get('患者姓名')):
            continue
        recs.append(r)
scr = pd.DataFrame(recs).reset_index(drop=True)
print('重建筛查行数:', len(scr))

coh = pd.read_excel(os.path.join(BASE, '清洗结果', 'cleaned_dataset.xlsx'),
                    sheet_name='血脂-糖网关联队列')
print('清洗后队列行数:', len(coh))

# ---------- verify row alignment before using positional merge ----------
pid_s = scr['患者编号'].astype(str).str.strip().tolist()
pid_c = coh['患者编号'].astype(str).str.strip().tolist()
same = sum(a == b for a, b in zip(pid_s, pid_c))
print('按行序 患者编号 一致: %d / %d (%.1f%%)' % (same, len(pid_c), 100 * same / len(pid_c)))
assert same == len(pid_c), 'row order mismatch -> must join, not align positionally'
print('行序校验通过，可按位置合并')

# ---------- outcomes ----------
Lc, Lg = zip(*[parse_eye(t) for t in scr['左眼糖网筛查评估结果（仅供参考）']])
Rc, Rg = zip(*[parse_eye(t) for t in scr['右眼糖网筛查评估结果（仅供参考）']])
Lg, Rg = np.array(Lg, float), np.array(Rg, float)
Lt = scr['左眼糖网筛查评估结果（仅供参考）'].tolist()
Rt = scr['右眼糖网筛查评估结果（仅供参考）'].tolist()

graded = ~np.isnan(Lg) | ~np.isnan(Rg)
sevA = np.nanmax(np.vstack([Lg, Rg]), axis=0)
gradeable = graded.astype(float)
any_dr = ((sevA >= 1) & graded).astype(float)
pdr = ((sevA == 4) & graded).astype(float)
vtpr = ((sevA >= 2) & graded).astype(float)

# ---- sensitivity A': fill ungraded eyes from 'visible DR' features / quantitation ----
Lg2 = Lg.copy(); Rg2 = Rg.copy()
for i in range(len(scr)):
    if np.isnan(Lg2[i]):
        Lg2[i] = fgrade(Lt[i]) if Lc[i] == 'visible' else (quant_grade(Lt[i]) if Lc[i] == 'quant' else np.nan)
    if np.isnan(Rg2[i]):
        Rg2[i] = fgrade(Rt[i]) if Rc[i] == 'visible' else (quant_grade(Rt[i]) if Rc[i] == 'quant' else np.nan)
graded2 = ~np.isnan(Lg2) | ~np.isnan(Rg2)
sevB = np.nanmax(np.vstack([Lg2, Rg2]), axis=0)
gradeable_b = graded2.astype(float)
any_dr_b = ((sevB >= 1) & graded2).astype(float)
pdr_b = ((sevB == 4) & graded2).astype(float)

# ---------- assemble dev_set_v2.csv ----------
out = pd.DataFrame()
out['id'] = coh['患者编号']
age = pd.to_numeric(coh['年龄'], errors='coerce')
out['age'] = age.where(age.between(14, 120))
out['male'] = (coh['性别'].astype(str).str.strip() == '男').astype('float')
for src, dst in [('TG', 'tg'), ('TC', 'tc'), ('HDL_C', 'hdl_c'), ('LDL_C', 'ldl_c'),
                 ('APO_A1', 'apo_a1'), ('APO_B', 'apo_b')]:
    out[dst] = pd.to_numeric(coh[src], errors='coerce')
out['dr_raw'] = coh['患者DR'].astype('float')
out['severity'] = sevA                    # corrected max-eye grade (0-4, NaN if ungraded)
out['gradeable'] = gradeable
out['any_dr'] = any_dr
out['pdr'] = pdr
out['vtpr_dr'] = vtpr
out['sev_Aprime'] = sevB
out['gradeable_Aprime'] = gradeable_b
out['any_dr_Aprime'] = any_dr_b
les = [('左眼出血数', 'L_hem_num'), ('左眼出血面积', 'L_hem_area'), ('左眼渗出数', 'L_exu_num'),
       ('左眼渗出面积', 'L_exu_area'), ('右眼出血数', 'R_hem_num'), ('右眼出血面积', 'R_hem_area'),
       ('右眼渗出数', 'R_exu_num'), ('右眼渗出面积', 'R_exu_area'), ('左眼微血管瘤', 'L_micro'),
       ('左眼硬性渗出', 'L_hardex'), ('左眼黄斑水肿', 'L_macular'), ('左眼新生血管', 'L_neov'),
       ('右眼微血管瘤', 'R_micro'), ('右眼硬性渗出', 'R_hardex'), ('右眼黄斑水肿', 'R_macular'),
       ('右眼新生血管', 'R_neov')]
for s, d_ in les:
    out[d_] = coh[s].astype('float') if s in coh.columns else np.nan

out.to_csv(os.path.join(OUT, 'dev_set_v2.csv'), index=False, encoding='utf-8')

# ---------- QC ----------
def qc(name, g, a, sevv, p):
    gg = g.astype(bool)
    print('\n--- %s ---' % name)
    print('  gradeable n=%d  events=%d  prevalence=%.1f%%  PDR=%d' % (gg.sum(), a[gg].sum(), 100 * a[gg].mean(), int(p[gg].sum())))
    print('  gradeable severity distribution:', {int(k): int(v) for k, v in
          pd.Series(sevv[gg]).value_counts().sort_index().items()})
    print('  6-lipid complete in gradeable: %d' % int(out.loc[gg, ['tg', 'tc', 'hdl_c', 'ldl_c', 'apo_a1', 'apo_b']].notna().all(axis=1).sum()))
    cc = gg & out[['age', 'male'] + ['tg', 'tc', 'hdl_c', 'ldl_c', 'apo_a1', 'apo_b']].notna().all(axis=1).values
    print('  complete cases (age/sex/6 lipids): n=%d events=%d (%.1f%%)' % (cc.sum(), int(a[cc].sum()), 100 * a[cc].mean()))
    return gg.sum(), a[gg].sum()

print('\n===== QC =====')
qc('option A  (primary)', gradeable, any_dr, sevA, pdr)
qc("option A' (sensitivity)", gradeable_b, any_dr_b, sevB, pdr_b)
qc('OLD (for reference)', out['severity'].notna().values * 0 + pd.read_csv(os.path.join(OUT, 'dev_set.csv'))['gradeable'].values,
   pd.read_csv(os.path.join(OUT, 'dev_set.csv'))['any_dr'].values,
   pd.read_csv(os.path.join(OUT, 'dev_set.csv'))['severity'].values,
   pd.read_csv(os.path.join(OUT, 'dev_set.csv'))['pdr'].values)
print('\ndev_set_v2.csv written:', out.shape)
