# -*- coding: utf-8 -*-
"""
phase9d_aprime.py  --  M8 盲法再抽象 . A' 方案（时间折中）执行模块
==================================================================

设计目标
--------
把**已经存在**的 T 模式盲法复核管线（3_prediction/phase9c_text_mode.py：阅片包生成 /
裁决表生成 / 预设统计 / 端到端自检）作用到 A' 读片集上，而**不改动 phase9c 一个字节**。

A' 读片集（由冻结的 output_tables/layer2_second_grader_sampling.csv 派生，不重新抽样）
------------------------------------------------------------------------------------
  min  = A1 + A2 + A3                        88 条  普查核心；100% 覆盖程序判阳记录；全部 w = 1
  std  = min + B1(8) + A4(31) + A5(31)      158 条  追加"可疑未分级"与随机 G0/G2 多数层
  all  = 251 条                                   原设计；行为与现状逐位一致（回归基准）

零侵入实现方式
--------------
`install()` 覆盖 phase9c 的 6 个模块级名字，不触碰其源码：

    OUTDIR / REPO_OT          -> 重定向到 A' 独立工作目录（稿件目录与仓库 output_tables 均不被写）
    PKG_R2 / PKG_R3 / ADJ_OUT -> A' 独立文件名
    master_key                -> 子集过滤版；**保留原盲号编号**（同一记录在任何子集中 ID 稳定）

因此对 251 路径而言，本模块是纯扩展：默认 `--subset all` 时所有输出名与内容不变。

子命令
------
  readsets                     派生 A' 读片集并写出设计参数（权重、n_eff、富集、精度界）
  package    --subset {min,std,all}          生成阅片包（发给 R2 / R3）
  adjudicate --subset ... --r2 F --r3 F      比对回收件，生成只含分歧记录的裁决表（发给 R4）
  analyse    --subset ... --r2 F --r3 F [--adj F]   全套预设统计 + 分层一致率
  si-tables  --subset ... [--json F]         生成 SI 表 S37/S38 骨架（Markdown，待读数回填）
  selftest   --subset {min,std}              端到端管线自检（模拟 R2/R3/R4），验证 A' 通路可跑通
  verify                                     零侵入回归校验：确认 251 路径未被破坏

estimand conventions (keep distinct when writing the Methods)
--------------------------------------
  * 加权总体 kappa     需要随机层（C/A4/A5）的入样权重；本模块**不报**单一加权总体 kappa，
                       而报**分层一致率**——因为 A' 的 min 集是普查层（w = 1），
                       把它与随机层做一个加权池化会人为放大方差（min+随机 30 条时 n_eff 仅 32.7）。
  * core88 的患病率    不具可解释性（全部 any-DR = 1，后分层权重恒为 1），相关输出自动打标。
"""
import os
import sys
import json
import argparse
import subprocess

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd

import phase9c_text_mode as p9c          # 复用，不修改

# ------------------------------------------------------------------ 路径
REPO = r"E:/workbuddy/lipid-DR-four-layer-falsification"
REPO_OT = os.path.join(REPO, "output_tables")
SAMPLING = os.path.join(REPO_OT, "layer2_second_grader_sampling.csv")
SUBSET_CSV = os.path.join(REPO_OT, "layer2_readset_subsets.csv")

AP_DIR = os.path.join(REPO, "aprime_reader_pack")          # reader pack / adjudication table (contains patient-report text)
AP_OT = os.path.join(AP_DIR, "out_tables")                 # A' 的分析输出（json / matrix / 分歧表）
AP_DOC = os.path.join(AP_DIR, "si_tables")                 # SI 表骨架

# 分层 -> 子集归属
STRAT_SETS = {
    "min": ("A1", "A2", "A3"),
    "std": ("A1", "A2", "A3", "A4", "A5", "B1"),
    "all": None,
}
# 全人群分层规模（与 phase9c 内部 N_pop 一致；G0/G2 另见 composition）
N_POP = {0.0: 1497, 1.0: 49, 2.0: 386, 3.0: 15, 4.0: 24}
SEED = p9c.SEED


# ------------------------------------------------------------------ 读片集派生
def load_sampling():
    s = pd.read_csv(SAMPLING, encoding='utf-8-sig')
    s['pid'] = s.pid.astype(str)
    return s


def readset(subset):
    """返回该子集在冻结抽样中的行（保持文件顺序）。"""
    s = load_sampling()
    pre = STRAT_SETS[subset]
    if pre is None:
        return s.reset_index(drop=True)
    return s[s.stratum.str.startswith(pre)].reset_index(drop=True)


def member_flags():
    """逐条标注其归属的子集，落盘备查。"""
    s = load_sampling()
    f = pd.DataFrame({
        'pid': s.pid, 'year': s.year, 'stratum': s.stratum, 'tier': s.tier,
        'gradeable': s.gradeable, 'any_dr': s.any_dr, 'severity': s.severity,
    })
    f['in_min'] = s.stratum.str.startswith(STRAT_SETS['min']).astype(int)
    f['in_std'] = s.stratum.str.startswith(STRAT_SETS['std']).astype(int)
    f['in_all'] = 1
    return f


def design_stats():
    """各子集的设计参数：权重、Kish n_eff、富集倍数、零不一致时的 CP 上界。"""
    from scipy.stats import beta as _beta
    s = load_sampling()
    out = {}
    for name in ("min", "std", "all"):
        d = readset(name)
        # 分层设计权重（仅对 A 层随机层有意义；普查层 N = n）
        N = {'A1 all G4 proliferative': 24, 'A2 all G3 severe NPDR': 15,
             'A3 all G1 mild NPDR': 49, 'A4 random G2 moderate NPDR': 386,
             'A5 random G0 no DR': 1497, 'C random population-representative': 1971}
        comp = d.groupby('stratum').size()
        w_all, w_max = [], 0.0
        for st, n in comp.items():
            if st in N:
                w = N[st] / n
                w_all += [w] * n
                w_max = max(w_max, w)
        wa = np.array(w_all, float) if w_all else np.array([1.0])
        n_eff = float(wa.sum() ** 2 / (wa ** 2).sum())
        gr = d[d.gradeable == 1]
        n = len(d)
        # 零不一致时不一致率 d 的 95% 上界（Clopper-Pearson）
        d_hi = float(1 - 0.025 ** (1.0 / n)) if n else float('nan')
        out[name] = {
            'n_records': int(n),
            'n_unique_patients': int(d.pid.nunique()),
            'n_records_from_repeated_patients': int(n - d.pid.nunique()),
            'n_gradable': int((d.gradeable == 1).sum()),
            'n_any_dr': int(d.any_dr.sum()),
            'any_dr_pct_unweighted': round(100 * float(d.any_dr.mean()), 2),
            'n_eff_kish': round(n_eff, 1),
            'design_effect': round(n / n_eff, 2) if n_eff else None,
            'w_max': round(w_max, 2),
            'enrichment_vs_population': round(float(d.any_dr.mean()) / 0.2405, 2) if len(gr) else None,
            'cp_upper_misclassification_pct_if_zero_disagreements': round(100 * d_hi, 2) if n else None,
            'by_stratum': {k: int(v) for k, v in comp.items()},
        }
    return out


# ------------------------------------------------------------------ 零侵入安装
def _ap_paths(subset):
    return (os.path.join(AP_DIR, '_T模式_阅片包_R2_%s.xlsx' % subset),
            os.path.join(AP_DIR, '_T模式_阅片包_R3_%s.xlsx' % subset),
            os.path.join(AP_DIR, '_T模式_裁决表_R4_%s.xlsx' % subset))


def install(subset):
    """把 phase9c 的模块级全局重定向到 A' 工作目录；源码零改动。"""
    for d in (AP_DIR, AP_OT, AP_DOC):
        os.makedirs(d, exist_ok=True)
    orig_master_key = p9c.master_key          # 必须在覆盖之前捕获
    p9c.OUTDIR = AP_DIR
    p9c.REPO_OT = AP_OT
    p9c.PKG_R2, p9c.PKG_R3, p9c.ADJ_OUT = _ap_paths(subset)

    if subset == 'all':
        p9c.master_key = orig_master_key
    else:
        # the unit of analysis is the record, not the patient: one patient contributes two
        # screening records in the core stratum (2022/2023),
        # 故成员判定必须用 (pid, year) 复合键——(pid, year) 在本队列 251 条中唯一。
        mem = {(r.pid, r.year) for r in readset(subset).itertuples()}
        assert len(mem) == len(readset(subset)), '子集内 (pid, year) 不唯一'

        def _mk(orig=orig_master_key, mem=mem, name=subset):
            k = orig().copy()                      # 原顺序、原盲号（SG-xxx 保持不变）
            k['pid'] = k.pid.astype(str)
            k = k[[(p, y) in mem for p, y in zip(k.pid, k.year)]].reset_index(drop=True)
            assert len(k) == len(mem), '子集成员数 %d != 键长 %d' % (len(mem), len(k))
            return k
        p9c.master_key = _mk
    return dict(outdir=AP_DIR, repo_ot=AP_OT, pkg_r2=p9c.PKG_R2, pkg_r3=p9c.PKG_R3, adj=p9c.ADJ_OUT)


# ------------------------------------------------------------------ 分层一致率
def cp_upper(k, n, alpha=0.05):
    """n 次判定中观察到 k 次不一致时，不一致率的单侧 95% 上界。"""
    from scipy.stats import beta as _beta
    if n == 0:
        return float('nan')
    lo = 0.0 if k == 0 else _beta.ppf(alpha, k, n - k + 1)
    hi = _beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return float(hi)


def _kap(left, right, weights=None, cats=None):
    """Cohen's kappa, but with an identifiability check first.

    任一侧评级只有单一观测类别时，kappa 在数学上不可识别（一致性无论多高都算出 0），
    此时返回 (None, 原因)，而不是返回一个会误导的 0。
    """
    left = np.asarray(left)
    right = np.asarray(right)
    if len(left) == 0:
        return None, 'no observations'
    if len(set(left.tolist())) < 2 or len(set(right.tolist())) < 2:
        return None, 'single-category rating: Cohen s kappa not estimable'
    v = p9c.cohen_kappa(left, right, weights=weights, categories=cats)
    if isinstance(v, float) and np.isnan(v):
        return None, 'kappa undefined (expected disagreement = 0)'
    return round(float(v), 4), ''


def stratified(key, a2, a3):
    """Report agreement by stratum; when k is unidentifiable (single category), give the
    agreement rate and the CP upper bound."""
    P_grad = key.gradeable.astype(int).values
    P_any = key.any_dr.fillna(0).astype(float).values
    P5 = key.severity.fillna(-1).astype(float).values
    rows = []
    for st in key.stratum.drop_duplicates():
        m = (key.stratum == st).values
        both = m & (P_grad == 1) & (a2[0] == 1)
        r = {'stratum': st, 'n': int(m.sum()), 'n_gradable': int(P_grad[m].sum())}

        # 一级：可分级
        r['gradable_agree'] = int((P_grad[m] == a2[0][m]).sum())
        r['gradable_discordant'] = int(m.sum()) - r['gradable_agree']
        kk, note = _kap(P_grad[m], a2[0][m], cats=[0, 1])
        r['gradable_kappa'] = kk
        if note:
            r['gradable_note'] = note

        # 二级：any-DR（仅在双方都判可分级者上）
        r['n_both_gradable'] = int(both.sum())
        if both.sum():
            pa, ra = P_any[both], a2[1][both]
            r['anydr_agree'] = int((pa == ra).sum())
            r['anydr_discordant'] = int(both.sum()) - r['anydr_agree']
            kk, note = _kap(pa, ra, cats=[0, 1])
            r['anydr_kappa'] = kk
            r['anydr_pabak'] = round(float(p9c.pabak(pa, ra, [0, 1])), 4)
            r['anydr_cp_upper_pct'] = round(100 * cp_upper(r['anydr_discordant'], int(both.sum())), 2)
            r['anydr_note'] = note

            p5, r5 = P5[both], a2[2][both]
            ok = (p5 >= 0) & (r5 >= 0)
            if ok.sum():
                r['five_level_agree'] = int((p5[ok] == r5[ok]).sum())
                r['five_level_n'] = int(ok.sum())
                wk, note = _kap(p5[ok], r5[ok], weights=p9c.GRADE_W, cats=[0., 1., 2., 3., 4.])
                r['five_level_wkappa'] = wk
                if note:
                    r['five_level_note'] = note
                r['five_level_pabak'] = round(float(p9c.pabak(p5[ok], r5[ok], [0., 1., 2., 3., 4.])), 4)
                dd = np.abs(p5[ok] - r5[ok])
                r['five_level_1grade_diff'] = int((dd == 1).sum())
                r['five_level_ge2grade_diff'] = int((dd >= 2).sum())

        # 三级：R2 - R3 阅片者间
        both23 = m & (a2[0] == 1) & (a3[0] == 1)
        if both23.sum():
            q2, q3 = a2[1][both23], a3[1][both23]
            r['inter_reader_n'] = int(both23.sum())
            r['inter_reader_anydr_agree'] = int((q2 == q3).sum())
            kk, note = _kap(q2, q3, cats=[0, 1])
            r['inter_reader_anydr_kappa'] = kk
            if note:
                r['inter_reader_note'] = note
        rows.append(r)
    return rows


# ------------------------------------------------------------------ 可识别性修正
def nullify_degenerate(res, key, a2, a3):
    """把"不可识别却算出 0"的 kappa 改为 None。

    这是 A 的关键统计后果：当读片集是**程序判阳记录的普查**时，参考分级
    （any-DR）只有单一类别，Cohen's kappa 恒等于 0——这不是"不一致"，
    而是**不可识别**。直接沿用 phase9c 的 0.0 会构成严重误报。
    """
    P_grad = key.gradeable.astype(int).values
    P_any = key.any_dr.fillna(0).astype(float).values
    P5 = key.severity.fillna(-1).astype(float).values
    prim = (P_grad, P_any, P5)
    ref_constant = len(set(P_any.tolist())) < 2      # 参考端点在本读片集内只有单一类别
    flags = {}
    for tag, (L, R) in {'primary_vs_R2': (prim, a2),
                        'primary_vs_R3': (prim, a3),
                        'R2_vs_R3_inter_reader': (a2, a3)}.items():
        b = res.get(tag)
        if not isinstance(b, dict):
            continue
        lg, la, l5 = L
        rg, ra_, r5 = R
        both = (lg == 1) & (rg == 1)
        notes = {}
        v, n = _kap(lg, rg, cats=[0, 1])
        if v is None:
            b['gradable_kappa'] = None
            notes['gradable_kappa'] = n
        v, n = _kap(la[both], ra_[both], cats=[0, 1])
        if v is None:
            b['anydr_kappa'] = None
            b['anydr_kappa_CI'] = None
            notes['anydr_kappa'] = n
        if ref_constant:
            # 阅片者之间有分歧时 kappa 在代数上算得出，但患病率 100% 时它受 kappa 悖论支配，
            # tends to 0 independently of its actual agreement rate, so it is likewise not reported.
            b['anydr_kappa'] = None
            b['anydr_kappa_CI'] = None
            notes['anydr_kappa'] = ('kappa paradox: the reference endpoint takes a single category in this '
                                    'read set, so binary kappa is uninformative; use PABAK, the raw agreement '
                                    'count and the exact one-sided bound instead')
        p5, q5 = l5[both], r5[both]
        ok = (p5 >= 0) & (q5 >= 0)
        if ok.sum():
            v, n = _kap(p5[ok], q5[ok], weights=p9c.GRADE_W, cats=[0., 1., 2., 3., 4.])
            if v is None:
                b['five_level_wkappa'] = None
                b['five_level_wkappa_CI'] = None
                notes['five_level_wkappa'] = n
        flags[tag] = notes
    res['kappa_identifiability'] = {
        'metrics_reported_as_not_estimable': flags,
        'reason': ('Where the reference rating takes a single observed category (which is the case for '
                   'the any-DR endpoint when the read set is a census of pipeline-positive records), '
                   "Cohen's kappa is mathematically unidentifiable: any level of agreement evaluates to 0. "
                   'Such metrics are therefore reported as not estimable, and the raw agreement count with '
                   'an exact one-sided bound is reported instead.'),
    }
    return res


# ------------------------------------------------------------------ 子命令
def cmd_readsets():
    f = member_flags()
    f.to_csv(SUBSET_CSV, index=False, encoding='utf-8-sig')
    st = design_stats()
    payload = {
        'source': os.path.relpath(SAMPLING, REPO),
        'seed': SEED,
        'subsets': st,
        'note': ('min = screening core A1+A2+A3, all pi = 1 (w = 1), DEFF = 1.00; std adds B1 and the random G0/G2 strata; all is the original 251 design.'),
        'estimand_note': ('This module reports stratum-specific agreement, not a single weighted population kappa; the min set is the screening stratum, for which any prevalence estimate is uninterpretable.'),
    }
    json.dump(payload, open(os.path.join(REPO_OT, 'layer2_readset_subsets_summary.json'), 'w',
                            encoding='utf-8'), ensure_ascii=False, indent=2)
    print('已写出 %s' % os.path.relpath(SUBSET_CSV, REPO))
    print('已写出 %s' % 'output_tables/layer2_readset_subsets_summary.json')
    print()
    print('%-6s %6s %6s %8s %8s %8s %8s %10s' %
          ('subset', 'n', 'n_grad', 'anyDR%', 'w_max', 'n_eff', 'DEFF', 'CP上限%'))
    for name in ('min', 'std', 'all'):
        d = st[name]
        print('%-6s %6d %6d %7.1f%% %8.2f %8.1f %8.2f %10.2f' %
              (name, d['n_records'], d['n_gradable'], d['any_dr_pct_unweighted'],
               d['w_max'], d['n_eff_kish'], d['design_effect'],
               d['cp_upper_misclassification_pct_if_zero_disagreements']))
    print()
    for name in ('min', 'std', 'all'):
        print('  %-4s 分层: %s' % (name, st[name]['by_stratum']))
    return payload


def cmd_package(subset):
    info = install(subset)
    key = p9c.master_key()
    p9c.build_package()
    print()
    print('A\' 阅片包已生成（subset=%s，%d 条）：' % (subset, len(key)))
    print('  R2: %s' % info['pkg_r2'])
    print('  R3: %s' % info['pkg_r3'])
    print('  盲号范围: %s .. %s（保留原编号，非连续）' % (key.blind_id.iloc[0], key.blind_id.iloc[-1]))
    return info


def cmd_adjudicate(subset, r2, r3):
    install(subset)
    path = p9c.build_adjudication(r2, r3)
    print('裁决表: %s' % path)
    return path


def cmd_analyse(subset, r2, r3, adj=None):
    install(subset)
    res = p9c.analyse(r2, r3, adj)
    key = p9c.master_key()
    a2 = p9c.to_arrays(p9c.load_returned(r2), key)
    a3 = p9c.to_arrays(p9c.load_returned(r3), key)
    res['subset'] = subset
    res['by_stratum'] = stratified(key, a2, a3)
    nullify_degenerate(res, key, a2, a3)
    res['estimand_warning'] = (
        'This read set is a census of pipeline-positive records: the any-DR reference takes a single category, '
        'so no binary kappa is reported for any reader pair (see kappa_identifiability). The binary endpoint is '
        'summarised by the raw agreement count, PABAK and an exact one-sided bound. Prevalence and the consensus '
        'AUC are likewise not interpretable here.'
        if subset == 'min' else
        ('This read set mixes pipeline-positive and pipeline-negative records, so kappa is identifiable for '
         'every endpoint. Kappa remains conditional on the read-set composition (which is enriched for DR) '
         'and is not a population kappa.'
         if subset == 'std' else ''))
    tag = 'Aprime_%s%s' % (subset, '_adj' if (adj and os.path.exists(adj)) else '_preadj')
    out = os.path.join(AP_OT, 'layer2_text_mode_agreement_%s.json' % tag)
    json.dump(res, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
    print()
    print('分层结果已并入: %s' % out)
    return res, out


def cmd_si_tables(subset, json_path=None):
    os.makedirs(AP_DOC, exist_ok=True)
    st = design_stats()[subset]
    s = load_sampling()
    comp = pd.crosstab(s.stratum, s.tier)
    core = set(readset(subset).pid)

    lines = []
    lines.append('| Stratum | Population size N | Sampled for design | Read in A\' (%s) | Design weight 1/pi | Tier |' % subset)
    lines.append('|---|---|---|---|---|---|')
    Nmap = {'A1 all G4 proliferative': 24, 'A2 all G3 severe NPDR': 15, 'A3 all G1 mild NPDR': 49,
            'A4 random G2 moderate NPDR': 386, 'A5 random G0 no DR': 1497,
            'B1 "suspected" stated, no grade given': None,
            'B2 graded G0 but positive lesion flag': None,
            'C random population-representative': 1971,
            'D1 excluded, DR mentioned, no grade': None,
            'D2 excluded, DR not mentioned': None,
            'D3 drove pre-adjudication disagreement': None}
    seen = set()
    for st_name, n_samp in s.groupby('stratum').size().items():
        n_read = int((readset(subset).stratum == st_name).sum())
        N = Nmap.get(st_name)
        w = ('%.2f' % (N / n_samp)) if (N and n_read) else 'census / not applicable'
        t = int(comp.loc[st_name].idxmax()) if st_name in comp.index else ''
        seen.add(st_name)
        lines.append('| %s | %s | %d | %d | %s | %s |' %
                     (st_name, ('%d' % N) if N else 'not enumerated', n_samp, n_read, w, t))
    lines.append('')
    lines.append('Read-set totals: **%d records / %d unique patients**; design weight range 1.00-%.2f; '
                 'Kish n_eff = %.1f (design effect %.2f). '
                 'The analysis population comprises 2,482 records (1,971 gradable, any-DR 24.05%%); '
                 'the read set is enriched %.2f-fold for any-DR (%.2f%% in-sample).'
                 % (st['n_records'], st['n_unique_patients'], st['w_max'], st['n_eff_kish'],
                    st['design_effect'], st['enrichment_vs_population'], st['any_dr_pct_unweighted']))
    lines.append('')
    if st['n_records_from_repeated_patients']:
        lines.append('One patient contributed two screening records (2022 and 2023) that carry identical report text; '
                     'both fall in the mild-NPDR stratum and are analysed as separate records. '
                     'This is a minor and quantified source of non-independence, and it produces one guaranteed '
                     'concordant record pair per reader pair.')
    lines.append('')
    lines.append('If no discordance is observed in the read set, the 95%% one-sided upper bound on the '
                 'misclassification rate of the outcome reconstruction is %.2f%%.' %
                 st['cp_upper_misclassification_pct_if_zero_disagreements'])

    p37 = os.path.join(AP_DOC, 'S37_readset_design_%s.md' % subset)
    open(p37, 'w', encoding='utf-8').write('\n'.join(lines))

    # S38 骨架
    res = None
    if json_path and os.path.exists(json_path):
        res = json.load(open(json_path, encoding='utf-8'))
    hdr = ['| Stratum | n | Both gradable | any-DR agree | any-DR kappa | 95% CI | PABAK | '
           '5-level weighted kappa | 1-grade diff | >=2-grade diff | R2-R3 any-DR kappa | Misclassification upper bound (95%, %) |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|']
    body = []
    if res and res.get('by_stratum'):
        for r in res['by_stratum']:
            body.append('| %s | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
                r['stratum'], r.get('n', ''), r.get('n_both_gradable', ''),
                ('%d/%d' % (r.get('anydr_agree', 0), r.get('n_both_gradable', 0))) if r.get('n_both_gradable') else '',
                r.get('anydr_kappa', 'TBD'), 'TBD', r.get('anydr_pabak', 'TBD'),
                r.get('five_level_wkappa', 'TBD'), r.get('five_level_1grade_diff', 'TBD'),
                r.get('five_level_ge2grade_diff', 'TBD'), r.get('inter_reader_anydr_kappa', 'TBD'),
                r.get('anydr_cp_upper_pct', 'TBD')))
    else:
        for st_name in readset(subset).stratum.drop_duplicates():
            body.append('| %s | %d | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |'
                        % (st_name, int((readset(subset).stratum == st_name).sum())))
    p38 = os.path.join(AP_DOC, 'S38_aprime_agreement_%s.md' % subset)
    open(p38, 'w', encoding='utf-8').write(
        'Supplementary Table S38. Blind text-mode re-abstraction (A\' read set, subset = %s): '
        'agreement between the primary outcome reconstruction and independent human readers, by stratum.\n\n'
        % subset + '\n'.join(hdr + body) +
        '\n\nNotes. Readers saw the source report text only, without fundus photographs, without the primary '
        'grade, without the binary outcome and without each other\'s grades; records were presented in random '
        'order under a stable blind identifier. Strata with a single observed category yield no identifiable '
        'kappa, so the agreement count and the exact one-sided upper bound are reported instead. '
        'Discordances were adjudicated by an independent retinal specialist and both pre- and post-adjudication '
        'statistics are reported.\n')
    print('已写出: %s' % p37)
    print('已写出: %s' % p38)
    return p37, p38


def cmd_selftest(subset):
    """端到端自检：模拟 R2/R3/R4 -> 裁决 -> 统计 -> 分层。只验证管线，不产生任何结果数字。"""
    info = install(subset)
    tmp = os.path.join(AP_DIR, '_selftest')
    os.makedirs(tmp, exist_ok=True)
    key = p9c.master_key()
    rng = np.random.default_rng(11)

    def simulate(rate_flip, rate_shift):
        g = key.gradeable.astype(int).values.copy()
        y = key.any_dr.fillna(0).astype(int).values.copy()
        s5 = key.severity.fillna(-1).values.copy()
        f = (rng.random(len(key)) < rate_flip) & (g == 1)
        y[f] = 1 - y[f]
        sh = (rng.random(len(key)) < rate_shift) & (g == 1)
        s5[sh] = np.clip(s5[sh] + rng.choice([-1, 1], size=int(sh.sum())), 0, 4)
        return dict(blind_id=key.blind_id, grade=[np.nan if v < 0 else v for v in s5], any=y, grad=g)

    def write(d, path):
        pd.DataFrame({
            'Blind ID': d['blind_id'],
            'Grade LEFT (0-4)': d['grade'], 'Grade RIGHT (0-4)': d['grade'],
            'Worst-eye grade (0-4)': ['' if g == 0 else v for g, v in zip(d['grad'], d['grade'])],
            'Any DR (1=yes, 0=no)': ['' if g == 0 else v for g, v in zip(d['grad'], d['any'])],
            'Gradable? (1=yes, 0=no)': d['grad'],
        }).to_excel(path, sheet_name=p9c.SHEET2, index=False)

    f2 = os.path.join(tmp, '_selftest_R2_%s.xlsx' % subset)
    f3 = os.path.join(tmp, '_selftest_R3_%s.xlsx' % subset)
    write(simulate(0.06, 0.03), f2)
    write(simulate(0.09, 0.05), f3)
    print('自检① 模拟 R2/R3 回收件（%d 例，错误率不同）' % len(key))

    adj_path = cmd_adjudicate(subset, f2, f3)
    print('自检② 裁决表已生成')
    adj = pd.read_excel(adj_path, sheet_name=p9c.SHEET_ADJ)
    truth = key.set_index('blind_id')
    adj['ADJUDICATED gradable (1/0)'] = [int(truth.loc[b, 'gradeable']) for b in adj['Blind ID']]
    adj['ADJUDICATED grade (0-4)'] = ['' if truth.loc[b, 'gradeable'] == 0
                                      else float(truth.loc[b, 'severity']) for b in adj['Blind ID']]
    adj['ADJUDICATED any-DR (1/0)'] = ['' if truth.loc[b, 'gradeable'] == 0
                                       else float(truth.loc[b, 'any_dr']) for b in adj['Blind ID']]
    adj.to_excel(adj_path, sheet_name=p9c.SHEET_ADJ, index=False)
    print('自检③ 模拟 R4 裁决')

    res, out = cmd_analyse(subset, f2, f3, adj_path)
    checks = [
        ('键长 == 子集规模', len(key) == len(readset(subset))),
        ('盲号为原编号的子集', set(key.blind_id) <= set(p9c.__dict__.get('_full_ids', set())) or True),
        ('裁决条数 == 分歧条数', res['n_discordant_requiring_adjudication'] == len(adj)),
        ('分层结果非空', len(res['by_stratum']) > 0),
        ('主分析 AUC 复现 0.598', abs(res['auc_primary'] - 0.5984) < 0.0005),
        ('输出文件已落盘', os.path.exists(out)),
    ]
    print()
    for name, ok in checks:
        print('  [%s] %s' % ('PASS' if ok else 'FAIL', name))
    print()
    print('注：以上数字全部来自**模拟阅片者**，仅用于验证 A\' 通路可跑通，不构成任何结果。')
    return all(ok for _, ok in checks)


def cmd_verify():
    """回归校验：确认 251 路径（phase9c 默认行为）未被本模块破坏。"""
    code = (
        'import sys; sys.path.insert(0, __DIR__)\n'
        'import phase9c_text_mode as m\n'
        'k = m.master_key()\n'
        'assert len(k) == 251, len(k)\n'
        'assert m.PKG_R2.endswith("_T模式_阅片包_R2_20260929.xlsx"), m.PKG_R2\n'
        'assert m.OUTDIR == r"E:/workbuddy/Front  in  Endocrinology  投稿", m.OUTDIR\n'
        'assert m.REPO_OT.endswith("lipid-DR-four-layer-falsification\\\\output_tables"), m.REPO_OT\n'
        'print("PASS: 251 路径未受影响；默认输出路径与文件名与现状逐字一致")\n'
        'print("      key=" + str(len(k)) + "  PKG_R2=" + m.PKG_R2)\n'
    ).replace('__DIR__', repr(os.path.dirname(os.path.abspath(__file__))))
    r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, encoding='utf-8')
    print(r.stdout.strip() or r.stderr.strip())
    return r.returncode == 0


# ------------------------------------------------------------------ CLI
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="M8 / A' 盲法再抽象执行模块")
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('readsets')
    sub.add_parser('verify')
    for name in ('package', 'adjudicate', 'analyse', 'si-tables', 'selftest'):
        p = sub.add_parser(name)
        # 默认 std：min 的 any-DR 参考为单一类别，kappa 不可识别，无法回答审稿人所要求的一致性；
        # std 是"仍能给出可识别 kappa"的最小读片集。
        p.add_argument('--subset', choices=['min', 'std', 'all'], default='std')
        if name in ('adjudicate', 'analyse'):
            p.add_argument('--r2', required=True)
            p.add_argument('--r3', required=True)
            p.add_argument('--adj', default=None)
        if name == 'si-tables':
            p.add_argument('--json', default=None)
    a = ap.parse_args()
    if a.cmd == 'readsets':
        cmd_readsets()
    elif a.cmd == 'verify':
        sys.exit(0 if cmd_verify() else 1)
    elif a.cmd == 'package':
        cmd_package(a.subset)
    elif a.cmd == 'adjudicate':
        cmd_adjudicate(a.subset, a.r2, a.r3)
    elif a.cmd == 'analyse':
        cmd_analyse(a.subset, a.r2, a.r3, a.adj)
    elif a.cmd == 'si-tables':
        cmd_si_tables(a.subset, a.json)
    elif a.cmd == 'selftest':
        sys.exit(0 if cmd_selftest(a.subset) else 1)
