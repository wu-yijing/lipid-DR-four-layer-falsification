# -*- coding: utf-8 -*-
"""
phase9b_second_grader_sampling.py
---------------------------------
为"第二阅片者盲法复核"生成抽样方案与工作簿（路线 D1 的确认性加固件）。

分层设计（定向层优先抽取，保证决策关键记录全部入样）
  Tier 1（必做，一致性核心）
    A1 全部 G4 增殖期            A2 全部 G3 重度 NPDR
    A3 全部 G1 轻度 NPDR         A4 随机 G2 中度 NPDR (31)
    A5 随机 G0 无 DR (31)
    B1 文本"可疑"但未给分级（全部）
    B2 判为 G0 但结构化病灶标志阳性（全部）
  Tier 2（建议）
    C  随机代表性样本 (40) —— 富集样本无法估患病率，须由它复现
    D3 驱动 D1 裁决前不一致的 6 条（全部）
  Tier 3（完整）
    D1 排除记录·文本提及 DR 但无分级 (10)
    D2 排除记录·文本未提及 DR (10)

关键事实（脚本内断言，防止定义漂移）
  * 病灶"计数"字段（L_hem_num 等）只有 386 条有值且全部落在不可分级记录上，
    不能用于内部一致性审计；必须用二值病灶"标志"字段（1,971 条可分级者全部有值）。
  * 病灶标志阳性但判为 G0 者恰为 27 条（= Supplementary Table S13 的 27/1,497）。
  * 不可分级记录中带明确分级者为 0 条（排除规则无内部矛盾）。
  * 文本含"可疑"者 17 条：6 条可分级·G0，11 条被排除；全部未给分级。

输出
  <投稿目录>/_作者待办_二次阅片抽样清单_v2.xlsx
    1_grader_worksheet_blinded  去盲工作表（随机顺序，无原分级，交给阅片者）
    2_master_key_internal       盲号对照（保留在研究团队，不得外发）
    3_design_notes              抽样设计说明
    4_analysis_template         预设统计与报告模板
  代码仓库 output_tables/ 落盘抽样明细与汇总 json。

运行：直接执行；随机种子固定，可复现。
"""
import sys, os, json
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

SRC = r"E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
REPO_OT = r"E:/workbuddy/lipid-DR-four-layer-falsification/output_tables"
OUTDIR = r"E:/workbuddy/Front  in  Endocrinology  投稿"
OUT = os.path.join(OUTDIR, "_作者待办_二次阅片抽样清单_v2.xlsx")
PRE_ADJ = os.path.join(REPO_OT, "layer2_independent_reextraction_discordances_preadjudication.csv")
SEED = 20260929

aud = pd.read_csv(f"{SRC}/severity_rebuild_audit.csv", low_memory=False)
dev = pd.read_csv(f"{SRC}/dev_set_v2.csv", low_memory=False)
assert len(aud) == len(dev), "audit 与 dev 行数不一致"
assert (aud['any_dr'].astype(float).fillna(-1).values ==
        dev['any_dr'].astype(float).fillna(-1).values).all(), "audit 与 dev 行序未对齐"

FLG = ['L_micro', 'L_hardex', 'L_macular', 'L_neov',
       'R_micro', 'R_hardex', 'R_macular', 'R_neov']
assert dev[FLG][dev.gradeable.eq(1)].notna().all().all(), "病灶标志列在可分级者中存在缺失"

D = pd.DataFrame({
    'pid': aud['pid'].astype(str),
    'year': aud['year'],
    'sex': aud['sex'],
    'age': aud['age'],
    'Ltxt': aud['Ltxt'].fillna(''),
    'Rtxt': aud['Rtxt'].fillna(''),
    'gradeable': dev['gradeable'].values,
    'any_dr': dev['any_dr'].values,
    'severity': dev['severity'].values,
})
D['les_flag'] = (dev[FLG].fillna(0) > 0).any(axis=1).values
D['txt'] = D.Ltxt + ' ' + D.Rtxt
D['suspected'] = D.txt.str.contains('可疑')
D['grade_words'] = D.txt.str.contains('轻度|中度|重度|非增殖性|增殖性')
D['mentioned'] = D.txt.str.contains('视网膜病变|微血管瘤|非增殖性|增殖性')
GRADE = {0: 'G0 no DR', 1: 'G1 mild NPDR', 2: 'G2 moderate NPDR',
         3: 'G3 severe NPDR', 4: 'G4 proliferative DR'}
D['grade_label'] = D.severity.map(GRADE)
G = D.gradeable.eq(1)
NG = ~G

# ---- 脚本内断言：定义漂移保护 ----
n_flag_g0 = int((G & D.any_dr.eq(0) & D.les_flag).sum())
assert n_flag_g0 == 27, "判为 G0 但病灶标志阳性者应为 27，实为 %d" % n_flag_g0
assert int((NG & D.grade_words).sum()) == 0, "存在被排除却带明确分级者"
assert int(D.suspected.sum()) == 17, "含'可疑'记录应为 17 条"
assert int((NG & D.suspected).sum()) == 11, "被排除的可疑记录应为 11 条"
n_ng_ment = int((NG & D.mentioned & ~D.grade_words).sum())
assert n_ng_ment == 83, "被排除但提及 DR 者应为 83，实为 %d" % n_ng_ment

pre_adj_pids = list(pd.read_csv(PRE_ADJ)['pid'].astype(str))
assert len(pre_adj_pids) == 6, "裁决前不一致记录应为 6 条"

rng = np.random.default_rng(SEED)
picked = {}
def take(mask, n, label):
    idx = [i for i in D.index[mask] if i not in picked]
    if n is not None and len(idx) > n:
        idx = list(rng.choice(idx, size=n, replace=False))
    for i in sorted(idx):
        picked[i] = label
    return idx

# ================= 定向层（先抽） =================
L_B2 = 'B2 graded G0 but positive lesion flag'
L_B1 = 'B1 "suspected" stated, no grade given'
L_D3 = 'D3 drove pre-adjudication disagreement'
take(G & D.any_dr.eq(0) & D.les_flag, None, L_B2)
take(D.suspected & ~D.grade_words, None, L_B1)
for i in D.index[D.pid.isin(pre_adj_pids)]:
    picked[i] = L_D3

# ================= Tier 1 核心分层 =================
for lab, sev, n in [('A1 all G4 proliferative', 4, None), ('A2 all G3 severe NPDR', 3, None),
                    ('A3 all G1 mild NPDR', 1, None), ('A4 random G2 moderate NPDR', 2, 31),
                    ('A5 random G0 no DR', 0, 31)]:
    take(G & D.severity.eq(sev), n, lab)

# ================= Tier 2 =================
take(G, 40, 'C random population-representative')

# ================= Tier 3 =================
take(NG & D.mentioned & ~D.grade_words, 10, 'D1 excluded, DR mentioned, no grade')
take(NG & ~D.mentioned, 10, 'D2 excluded, DR not mentioned')

sel = D.loc[sorted(picked)].copy()
sel['stratum'] = [picked[i] for i in sel.index]
TIER = {'A1': 1, 'A2': 1, 'A3': 1, 'A4': 1, 'A5': 1, 'B1': 1, 'B2': 1,
        'C': 2, 'D3': 2, 'D1': 3, 'D2': 3}
sel['tier'] = sel.stratum.str.split().str[0].map(TIER)
assert sel.tier.notna().all(), '分层标签未映射到 Tier：%s' % sorted(sel.loc[sel.tier.isna(),'stratum'].unique())

print("=== 分层结果 ===")
print(sel.groupby(['tier', 'stratum']).size().to_string())
print("\n总例数: %d | 可分级 %d | 不可分级 %d"
      % (len(sel), int(sel.gradeable.sum()), int(sel.gradeable.eq(0).sum())))
for t in [1, 2, 3]:
    print("  Tier %d: %3d 例  (累计 %d)" % (t, int((sel.tier == t).sum()), int((sel.tier <= t).sum())))
print("\n=== 年份覆盖 ===")
print(sel.year.value_counts().sort_index().to_string())
print("\n=== 核心集 A 五级覆盖 vs 全人群 ===")
for sev, lab in GRADE.items():
    na = int(((sel.severity == sev) & sel.stratum.str.startswith('A')).sum())
    npop = int((D.loc[G, 'severity'] == sev).sum())
    print("  %-20s 抽 %3d / 人群 %4d (%.1f%%)" % (lab, na, npop, 100 * na / max(npop, 1)))
n_empty = int((sel.Ltxt.str.strip().eq('') & sel.Rtxt.str.strip().eq('')).sum())
n_half = int(((sel.Ltxt.str.strip().eq('') | sel.Rtxt.str.strip().eq('')) & ~(sel.Ltxt.str.strip().eq('') & sel.Rtxt.str.strip().eq(''))).sum())
print("\n报告中文本为空的记录（文本模式无法判读，需按影像模式补）: 双眼皆空 %d 条 / 单眼为空 %d 条" % (n_empty, n_half))
print("\n富集后 DR 比例（用于说明 kappa 的患病率依赖性）: %.1f%%"
      % (100 * sel.loc[sel.gradeable.eq(1), 'any_dr'].mean()))

# ================= 工作簿 =================
shuf = sel.sample(frac=1, random_state=SEED).reset_index(drop=True)
shuf['blind_id'] = ['SG-%03d' % (i + 1) for i in range(len(shuf))]

sheet1 = pd.DataFrame({
    'Blind ID': shuf.blind_id,
    'Screening year': shuf.year,
    'Left eye: original report text (Chinese)': shuf.Ltxt,
    'Right eye: original report text (Chinese)': shuf.Rtxt,
    'Fundus image file / archive ID (to be filled when retrieving)': '',
    'Grade LEFT (0-4)': '',
    'Grade RIGHT (0-4)': '',
    'Worst-eye grade for this record (0-4)': '',
    'Any DR for this record (1 = yes, 0 = no)': '',
    'Record gradable? (1 = yes, 0 = not gradable)': '',
    'Image/text quality (1 good / 2 fair / 3 poor)': '',
    'Comment (optional)': '',
})

sheet2 = pd.DataFrame({
    'Blind ID': shuf.blind_id,
    'Tier': shuf.tier,
    'Stratum': shuf.stratum,
    'Patient identifier': shuf.pid,
    'Screening year': shuf.year,
    'Age': shuf.age,
    'Sex': shuf.sex,
    'Gradeable (primary)': shuf.gradeable,
    'Primary worst-eye grade': shuf.severity,
    'Primary any-DR': shuf.any_dr,
    'Primary grade label': shuf.grade_label,
    'Positive lesion flag (primary)': shuf.les_flag,
    'Grader worst-eye grade': '',
    'Grader any-DR': '',
    'Grader gradable': '',
    'Grade agrees? (enter formula)': '',
    'any-DR agrees? (enter formula)': '',
    'Adjudicated grade (third reader, disagreements only)': '',
}).sort_values('Blind ID').reset_index(drop=True)

sheet3 = pd.DataFrame({
    'Item': [
        'Tier 1 - what it answers', 'Tier 2 - what it adds', 'Tier 3 - what it adds',
        'Stratum A (150)', 'Stratum B1 (17)', 'Stratum B2 (27)', 'Stratum D3 (6)',
        'Stratum C (40)', 'Stratum D1 / D2 (20)', 'Why the sample is enriched',
        'Why prevalence needs stratum C', 'Blinding', 'Reading order', 'Grading standard',
        'Photographs are not on this workstation', 'Reproducibility', 'Data protection',
    ],
    'Design': [
        'Agreement between the original reader and the second reader: quadratic-weighted kappa for the five-level ICDR grade and unweighted kappa for any-DR, on a sample that contains every rare-grade record.',
        'A random sample of the analysis population, to show that the 24.05% prevalence and the layer-2 discrimination are reproducible rather than artefacts of the enriched design; plus the six records that produced the pre-adjudication disagreement in the independent re-extraction.',
        'Records the primary pipeline excluded, to test whether "not gradable" is a property of the images or of the source text.',
        'All 24 proliferative, all 15 severe NPDR and all 49 mild NPDR records (the rare grades), plus 31 randomly chosen moderate NPDR and 31 randomly chosen no-DR records.',
        'Every record whose report states "suspected diabetic retinopathy" (可疑) without an ICDR grade - 17 records. These are the records where the interpretation of a verbal hedge decides case status.',
        'Every record graded no-DR that nevertheless carries a positive structured lesion descriptor - 27 records. This is the only internal inconsistency found in the outcome derivation (Supplementary Table S13).',
        'The six records that the independent re-extraction initially graded differently before the pre-specified adjudication rule was applied (Supplementary Table S14b).',
        'Simple random sample of the analysis population, used for prevalence and discrimination checks.',
        'D1: ten excluded records whose text mentions diabetic retinopathy without a grade. D2: ten excluded records that do not mention it. Together they test whether the exclusion rule removed records that a reader could have graded.',
        'Rare grades would contribute almost nothing to a simple random sample: proliferative DR is 1.2% of the cohort (24/1,971), so a random sample of 150 would contain about two such records and no usable estimate of agreement at the severe end. Enrichment buys precision where the clinical decisions are made.',
        'An enriched sample cannot estimate prevalence without weights. Stratum C is drawn at random so that the weighted prevalence can be compared with the 24.05% reported, and so that a systematic reclassification in one direction would be visible.',
        'The grader receives sheet 1 only: no primary grade, no stratum, no lipid values, no outcome. Records appear in random order and carry only a blind ID. The grader must not see sheets 2-4, the code repository or the manuscript.',
        'One sitting, random order, each record judged independently, no returning to earlier records.',
        'ICDR severity scale: 0 no apparent retinopathy, 1 mild NPDR, 2 moderate NPDR, 3 severe NPDR, 4 proliferative DR. The record-level grade is the worst eye.',
        'The fundus photographs are not stored on the analysis workstation. Whoever retrieves them from the screening system fills the image column of sheet 1 (or replaces the text columns with image links) before the sheet is sent out.',
        'Randomisation seed 20260929. Running phase9b_second_grader_sampling.py again reproduces this exact list.',
        'The blind ID is the only identifier in sheet 1. The key that links it to the patient identifier, the stratum and the primary grade is sheet 2 and stays with the study team.',
    ],
})

sheet4 = pd.DataFrame({
    'Analysis (pre-specified)': [
        'Five-level ICDR agreement', 'Binary any-DR agreement', 'Percent agreement and PABAK',
        'Not-gradable agreement', 'Prevalence reproducibility', 'Disagreement direction',
        'Effect on layer-2 discrimination', 'Worst-case sensitivity', 'Adjudication', 'Reporting',
    ],
    'Statistic and rule': [
        'Quadratic-weighted Cohen kappa with a 2,000-replicate bootstrap 95% CI, over records both readers call gradable.',
        'Unweighted Cohen kappa with bootstrap 95% CI over all sampled records, with "not gradable" carried as a separate category.',
        'Report raw agreement and the prevalence-adjusted bias-adjusted kappa alongside kappa. In an enriched sample kappa falls as prevalence departs from 50%, so kappa alone understates agreement and a reviewer will ask.',
        'Kappa for the binary gradable / not-gradable decision, plus the number of records the second reader could grade that the pipeline excluded.',
        'Weight by stratum so the sample reproduces the population composition; report the weighted prevalence with its CI against 24.05%.',
        'Five-level disagreement matrix; report how many disagreements are one grade apart and how many are two or more.',
        'Recompute the layer-2 complete-case AUC using the adjudicated grade in place of the primary grade. Compute the AUC at which the conclusion would change: with 367 events, 0.579 would have to move through roughly 0.70.',
        'Reclassify every disagreement against our reconstruction and recompute prevalence and AUC.',
        'A third reader (retina specialist, not one of the two above) resolves disagreements. Report both the pre-adjudication and the post-adjudication statistics.',
        'Report in Methods (design, blinding, grading standard, sample size rationale) and as one supplementary table; state explicitly that the sample was stratified and that kappa therefore refers to the enriched sample.',
    ]})

with pd.ExcelWriter(OUT, engine='openpyxl') as xw:
    sheet1.to_excel(xw, sheet_name='1_grader_worksheet_blinded', index=False)
    sheet2.to_excel(xw, sheet_name='2_master_key_internal', index=False)
    sheet3.to_excel(xw, sheet_name='3_design_notes', index=False)
    sheet4.to_excel(xw, sheet_name='4_analysis_template', index=False)

wb = openpyxl.load_workbook(OUT)
for ws in wb.worksheets:
    ws.freeze_panes = 'A2'
    for c in ws[1]:
        c.font = Font(bold=True)
        c.fill = PatternFill('solid', fgColor='DDEBF7')
        c.alignment = Alignment(vertical='center', wrap_text=True)
    for j in range(1, ws.max_column + 1):
        mx = max((len(str(c.value or '')) for c in ws[get_column_letter(j)]), default=10)
        ws.column_dimensions[get_column_letter(j)].width = min(max(mx + 2, 10), 46)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            if isinstance(c.value, str) and len(c.value) > 60:
                c.alignment = Alignment(wrap_text=True, vertical='top')
wb.save(OUT)

sel.to_csv(os.path.join(REPO_OT, 'layer2_second_grader_sampling.csv'), index=False, encoding='utf-8-sig')
summary = {
    'total': int(len(sel)),
    'by_tier': {str(t): int((sel.tier == t).sum()) for t in [1, 2, 3]},
    'cumulative_tier_1_2': int((sel.tier <= 2).sum()),
    'by_stratum': {k: int(v) for k, v in sel.stratum.value_counts().items()},
    'gradable': int(sel.gradeable.sum()),
    'not_gradable': int(sel.gradeable.eq(0).sum()),
    'enriched_DR_prevalence_pct': round(float(100 * sel.loc[sel.gradeable.eq(1), 'any_dr'].mean()), 2),
    'seed': SEED,
}
json.dump(summary, open(os.path.join(REPO_OT, 'layer2_second_grader_sampling_summary.json'), 'w'),
          ensure_ascii=False, indent=2)
print("\n已写出:", OUT)
print("已写出:", os.path.join(REPO_OT, 'layer2_second_grader_sampling.csv'))
