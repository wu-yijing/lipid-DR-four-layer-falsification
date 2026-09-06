#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
阶段 0.3 - 合并 NHANES 2005-2006 + 2007-2008 为"验证集 B"并计算衍生变量 (任务3)
输出:
  <项目>/外部验证_NHANES/验证集B/NHANES_验证集B.csv
  <项目>/外部验证_NHANES/验证集B/data_dictionary.md
  <项目>/外部验证_NHANES/验证集B/qc_report.md
合并键: SEQN (周期内固定); 跨周期用 cycle+SEQN 唯一标识。
DR 结局(基于 OPXRET worse-eye 分级 OPDURL4):
  any_dr  = OPDURL4>=2  (任意 DR: 轻度NPDR及以上)
  vtpr_dr = OPDURL4>=3  (视力威胁型: 中/重度NPDR + PDR)
  pdr     = OPDURL4==4  (增殖期 DR)
eGFR: CKD-EPI 2009 (race-based) 与 2021 (race-free) 双算。
"""
import os, numpy as np, pandas as pd, pyreadstat

RAW = r"E:/workbuddy/2026-09-04-18-19-16/nhanes_raw_fixed"
OUTDIR = r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/外部验证_NHANES/验证集B"
os.makedirs(OUTDIR, exist_ok=True)

# 组件 -> (周期后缀列表, 需提取的列)
SPEC = {
    "DEMO":   (["D", "E"], ["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "WTMEC2YR"]),
    "TCHOL":  (["D", "E"], ["SEQN", "LBXTC"]),
    "HDL":    (["D", "E"], ["SEQN", "LBDHDD"]),
    "TRIGLY": (["D", "E"], ["SEQN", "LBXTR", "LBDLDL", "LBXAPB"]),   # LBXAPB=ApoB(两周期均有)
    "GHB":    (["D", "E"], ["SEQN", "LBXGH"]),
    "GLU":    (["D", "E"], ["SEQN", "LBXGLU"]),
    "DIQ":    (["D", "E"], ["SEQN", "DIQ010", "DID040", "DIQ050", "DIQ080", "DID070"]),
    "BIOPRO": (["D", "E"], ["SEQN", "LBXSCR"]),
    "BMX":    (["D", "E"], ["SEQN", "BMXBMI"]),
    "BPX":    (["D", "E"], ["SEQN", "BPXSY1", "BPXDI1"]),
    "SMQ":    (["D", "E"], ["SEQN", "SMQ020"]),
    "OPXRET": (["D", "E"], ["SEQN", "OPDURET", "OPDURL4", "OPDUARMA", "OPDUMEC"]),
}

def load(comp, suffix):
    path = os.path.join(RAW, f"{comp}_{suffix}.xpt")
    cols = SPEC[comp][1]
    df, _ = pyreadstat.read_xport(path, usecols=cols)
    return df

# 每个周期内: 以 OPXRET(有眼底分级者)为锚, 左连其余组件; 再纵向拼接两周期
CYCLES = {"D": "2005-2006", "E": "2007-2008"}
frames = []
for suff, cyc in CYCLES.items():
    anchor = load("OPXRET", suff).rename(columns={"SEQN": "seqn"})
    cyc_df = anchor.copy()           # 锚: 仅保留有眼底分级的受试者
    print(f"[load] {cyc} OPXRET(锚): {anchor.shape}")
    for comp, (suffs, cols) in SPEC.items():
        if comp == "OPXRET" or suff not in suffs:
            continue
        d = load(comp, suff).rename(columns={"SEQN": "seqn"})
        cyc_df = cyc_df.merge(d, on="seqn", how="left")
        print(f"[merge] {cyc} {comp}: -> {cyc_df.shape}")
    cyc_df["_cycle"] = cyc
    frames.append(cyc_df)
merged = pd.concat(frames, ignore_index=True)
print(f"\n[concat] 两周期合并 n={len(merged)} cols={merged.shape[1]}")

# ---------- 衍生变量 ----------
m = merged
# 人口学
m["age"] = m["RIDAGEYR"]
m["male"] = (m["RIAGENDR"] == 1).astype("float")
m["race"] = m["RIDRETH1"]  # 1=MexAm 2=OthHis 3=NHW 4=NHB 5=Other
m["black"] = (m["race"] == 4).astype("float")

# 血脂 (mg/dL)
m["tc"] = m["LBXTC"]
m["hdl"] = m["LBDHDD"]
m["tg"] = m["LBXTR"]
m["ldl"] = m["LBDLDL"]
m["apob"] = m["LBXAPB"]  # 两周期均有(2005-2006来自TRIGLY, 2007-2008来自APOB)
# 非HDL-C
m["nonhdl"] = m["tc"] - m["hdl"]

# 血糖/糖尿病
m["hba1c_pct"] = m["LBXGH"]          # %
m["glucose"] = m["LBXGLU"]           # mg/dL
m["diabetes"] = (m["DIQ010"] == 1).astype("float")   # 医生告知患糖尿病
# 糖尿病病程(年) = 当前年龄 - 确诊年龄
m["dm_duration"] = np.where(
    (m["diabetes"] == 1) & (m["DID040"] >= 1) & (m["DID040"] < 90) & (m["DID040"] <= m["age"]),
    (m["age"] - m["DID040"]).clip(lower=0), np.nan)
m["insulin_use"] = (m["DIQ050"] == 1).astype("float")
m["oral_med"] = (m["DIQ080"] == 1).astype("float")

# 肾: eGFR (CKD-EPI 2009 race-based; 2021 race-free)
scr = m["LBXSCR"]
age = m["age"]
female = (m["male"] == 0)
def ckdepi2009(scr, age, female, black):
    k = np.where(female, 0.7, 0.9)
    a = np.where(female, -0.329, -0.411)
    minScr = np.minimum(scr, k)
    maxScr = np.maximum(scr, k)
    eGFR = 141 * (minScr / k) ** a * (maxScr / k) ** -1.209 * 0.993 ** age
    eGFR = np.where(female, eGFR * 1.018, eGFR)
    eGFR = np.where(black == 1, eGFR * 1.159, eGFR)
    return eGFR
def ckdepi2021(scr, age, female):
    k = np.where(female, 0.7, 0.9)
    a = np.where(female, -0.241, -0.302)
    minScr = np.minimum(scr, k)
    maxScr = np.maximum(scr, k)
    eGFR = 142 * (minScr / k) ** a * (maxScr / k) ** -1.200 * 0.9938 ** age
    eGFR = np.where(female, eGFR * 1.012, eGFR)
    return eGFR
m["egfr_ckdepi2009"] = ckdepi2009(scr, age, female, m["black"])
m["egfr_ckdepi2021"] = ckdepi2021(scr, age, female)
m["ckd_stage3plus"] = (m["egfr_ckdepi2009"] < 60).astype("float")

# 体测/生活
m["bmi"] = m["BMXBMI"]
m["sys_bp"] = m["BPXSY1"]
m["dia_bp"] = m["BPXDI1"]
m["smoker_100"] = (m["SMQ020"] == 1).astype("float")   # 终生吸烟>=100支

# DR 结局 (worse eye)
m["any_dr"] = (m["OPDURL4"] >= 2).astype("float")
m["vtpr_dr"] = (m["OPDURL4"] >= 3).astype("float")     # 视力威胁型(中/重度NPDR+PDR)
m["pdr"] = (m["OPDURL4"] == 4).astype("float")         # 增殖期
m["dr_level4"] = m["OPDURL4"]                          # 1-4
m["dr_level_detail"] = m["OPDURET"]                    # 细分级
m["macular_edema_center"] = (m["OPDUMEC"] == 1).astype("float")  # 中心性黄斑水肿

# ---------- 导出 ----------
keep_cols = ["seqn", "_cycle",
             "age", "male", "race", "black",
             "tc", "hdl", "tg", "ldl", "apob", "nonhdl",
             "hba1c_pct", "glucose", "diabetes", "dm_duration", "insulin_use", "oral_med",
             "LBXSCR", "egfr_ckdepi2009", "egfr_ckdepi2021", "ckd_stage3plus",
             "bmi", "sys_bp", "dia_bp", "smoker_100",
             "any_dr", "vtpr_dr", "pdr", "dr_level4", "dr_level_detail", "macular_edema_center",
             "OPDUARMA", "WTMEC2YR"]
out = m[keep_cols].copy()
csv_path = os.path.join(OUTDIR, "NHANES_验证集B.csv")
out.to_csv(csv_path, index=False)
print(f"\n[save] {csv_path}  shape={out.shape}")

# ---------- QC 报告 ----------
def pct(x): return f"{100*x:.1f}%"
qc = []
qc.append(f"# NHANES 验证集 B —— QC 报告")
qc.append(f"生成时间: 自动  | 数据周期: 2005-2006 + 2007-2008")
qc.append(f"")
qc.append(f"## 样本量")
qc.append(f"- 含眼底分级(OPXRET)总人数: **{len(out):,}**")
for cyc, g in out.groupby("_cycle"):
    qc.append(f"  - {cyc}: {len(g):,} 人")
qc.append(f"")
qc.append(f"## DR 患病率 (基于 OPDURL4 worse-eye)")
for lab, col in [("任意DR (OPDURL4>=2)", "any_dr"), ("视力威胁型 (OPDURL4>=3)", "vtpr_dr"), ("增殖期PDR (==4)", "pdr")]:
    sub = out.dropna(subset=[col])
    prev = sub[col].mean() if len(sub) else float('nan')
    qc.append(f"- {lab}: {sub[col].sum():.0f}/{len(sub)} = {pct(prev) if len(sub) else 'NA'}")
qc.append(f"- 中心性黄斑水肿 (OPDUMEC==1): {out['macular_edema_center'].sum():.0f}")
qc.append(f"")
qc.append(f"## 糖尿病状态")
dm = out.dropna(subset=["diabetes"])
qc.append(f"- 自报糖尿病 (DIQ010==1): {out['diabetes'].sum():.0f}/{len(out.dropna(subset=['diabetes']))} = {pct(out['diabetes'].mean())}")
qc.append(f"- 其中胰岛素治疗: {out['insulin_use'].sum():.0f}")
qc.append(f"- 其中口服药: {out['oral_med'].sum():.0f}")
qc.append(f"- 平均病程(年, 有值): {out['dm_duration'].mean():.1f} (n={out['dm_duration'].notna().sum()})")
qc.append(f"")
qc.append(f"## 关键预测因子缺失率")
for col in ["tc","hdl","tg","ldl","apob","hba1c_pct","LBXSCR","bmi","sys_bp","smoker_100","dm_duration"]:
    nn = out[col].notna().sum()
    qc.append(f"- {col}: 有值 {nn:,} ({pct(nn/len(out))}), 缺失 {len(out)-nn:,}")
qc.append(f"")
qc.append(f"## 人口学分布")
qc.append(f"- 年龄: 中位 {out['age'].median():.0f}, 范围 {out['age'].min():.0f}-{out['age'].max():.0f}")
qc.append(f"- 男性比例: {pct(out['male'].mean())}")
qc.append(f"- 非西班牙裔黑人比例: {pct(out['black'].mean())}")
qc.append(f"- eGFR(CKD-EPI2009) 中位: {out['egfr_ckdepi2009'].median():.0f}")
qc.append(f"- CKD 3期+ (<60): {pct(out['ckd_stage3plus'].mean())}")
qc_path = os.path.join(OUTDIR, "qc_report.md")
open(qc_path, "w", encoding="utf-8").write("\n".join(qc))
print(f"[save] {qc_path}")

# 数据字典
dic = []
dic.append("# NHANES 验证集 B —— 数据字典")
dic.append("")
dic.append("| 变量 | 含义 | 来源组件.原始变量 | 单位/编码 |")
dic.append("|---|---|---|---|")
rows = [
 ("age","年龄","DEMO.RIDAGEYR","岁"),
 ("male","男性(1)/女性(0)","DEMO.RIAGENDR","1=男 2=女→male"),
 ("race","种族","DEMO.RIDRETH1","1=MexAm 2=OthHis 3=NHW 4=NHB 5=Other"),
 ("black","非西班牙裔黑人","派生","race==4"),
 ("tc","总胆固醇","TCHOL.LBXTC","mg/dL"),
 ("hdl","HDL胆固醇","HDL.LBDHDD","mg/dL"),
 ("tg","甘油三酯","TRIGLY.LBXTR","mg/dL"),
 ("ldl","LDL胆固醇(计算)","TRIGLY.LBDLDL","mg/dL"),
 ("apob","载脂蛋白B","TRIGLY/APOB.LBXAPB","mg/dL"),
 ("nonhdl","非HDL-C","派生","tc-hdl"),
 ("hba1c_pct","糖化血红蛋白","GHB.LBXGH","%"),
 ("glucose","血糖","GLU.LBXGLU","mg/dL"),
 ("diabetes","自报糖尿病","DIQ.DIQ010","1=是 2=否"),
 ("dm_duration","糖尿病病程","DIQ.DID040派生","年=age-确诊年龄"),
 ("insulin_use","胰岛素治疗","DIQ.DIQ050","1=是"),
 ("oral_med","口服降糖药","DIQ.DIQ080","1=是"),
 ("LBXSCR","血清肌酐","BIOPRO.LBXSCR","mg/dL"),
 ("egfr_ckdepi2009","eGFR(CKD-EPI2009)","派生","mL/min/1.73m2"),
 ("egfr_ckdepi2021","eGFR(CKD-EPI2021)","派生","mL/min/1.73m2"),
 ("ckd_stage3plus","CKD 3期+","派生","egfr<60"),
 ("bmi","体质指数","BMX.BMXBMI","kg/m2"),
 ("sys_bp","收缩压","BPX.BPXSY1","mmHg"),
 ("dia_bp","舒张压","BPX.BPXDI1","mmHg"),
 ("smoker_100","终生吸烟≥100支","SMQ.SMQ020","1=是"),
 ("any_dr","任意DR(结局1)","OPXRET.OPDURL4>=2","0/1"),
 ("vtpr_dr","视力威胁型DR(结局2)","OPXRET.OPDURL4>=3","0/1"),
 ("pdr","增殖期DR(结局3)","OPXRET.OPDURL4==4","0/1"),
 ("dr_level4","DR 4级严重度","OPXRET.OPDURL4","1-4"),
 ("dr_level_detail","DR 细分级","OPXRET.OPDURET","数值编码"),
 ("macular_edema_center","中心性黄斑水肿","OPXRET.OPDUMEC","0/1"),
]
for r in rows:
    dic.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} |")
dic_path = os.path.join(OUTDIR, "data_dictionary.md")
open(dic_path, "w", encoding="utf-8").write("\n".join(dic))
print(f"[save] {dic_path}")
print("\n=== DONE 任务3 ===")
