"""Phase 2.1: export Chinese-xlsx dev cohort to ASCII-path CSV for R.
Dev cohort = '血脂-糖网关联队列' in cleaned_dataset.xlsx (n=2482).
R cannot read GBK Chinese paths -> Python reads xlsx, writes ASCII csv.
"""
import os, pandas as pd, numpy as np

BASE = "C:/Users/Administrator/Desktop/血脂，糖网等数据清洗"
SRC = os.path.join(BASE, "清洗结果", "cleaned_dataset.xlsx")
NHANES_B = os.path.join(BASE, "外部验证_NHANES", "验证集B", "NHANES_验证集B.csv")
OUTDIR = "E:/workbuddy/2026-09-04-18-19-16/phase2_dev"
EXTDIR = "E:/workbuddy/2026-09-04-18-19-16/phase2_ext"
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(EXTDIR, exist_ok=True)

coh = pd.read_excel(SRC, sheet_name="血脂-糖网关联队列")
print("raw dev cohort:", coh.shape)

out = pd.DataFrame()
out["id"] = coh["患者编号"]
# age: clip impossible values (raw max 1069 is data entry error)
age = pd.to_numeric(coh["年龄"], errors="coerce")
out["age"] = age.where(age.between(14, 120))
out["male"] = (coh["性别"].astype(str).str.strip() == "男").astype("float")

# 6 lipids (raw; tg also stored pre-log for log2 transform in R)
lip_map = [("TG","tg"),("TC","tc"),("HDL_C","hdl_c"),("LDL_C","ldl_c"),
           ("APO_A1","apo_a1"),("APO_B","apo_b")]
for src, dst in lip_map:
    out[dst] = pd.to_numeric(coh[src], errors="coerce")

# outcomes
out["dr_raw"] = coh["患者DR"].astype("float")       # raw bool (incl ungradable marked True)
sev = pd.to_numeric(coh["患者最重眼严重度"], errors="coerce")
out["severity"] = sev
out["gradeable"] = (~sev.isin([4, 9])).astype("float")  # exclude ungradable (4,9)
grad = out["gradeable"].astype(bool)
# PRIMARY any-DR = worst-eye grade >=1 among gradeable (基于最重眼分级, 不含不可分级误标)
out["any_dr"] = ((sev >= 1) & grad).astype("float")
out["pdr"] = ((sev == 3) & grad).astype("float")        # proliferative DR
out["vtpr_dr"] = ((sev >= 2) & grad).astype("float")    # moderate+ / threat-to-vision

# quantitative lesions (stage 3 1b备用)
les_map = [("左眼出血数","L_hem_num"),("左眼出血面积","L_hem_area"),("左眼渗出数","L_exu_num"),
           ("左眼渗出面积","L_exu_area"),("右眼出血数","R_hem_num"),("右眼出血面积","R_hem_area"),
           ("右眼渗出数","R_exu_num"),("右眼渗出面积","R_exu_area"),
           ("左眼微血管瘤","L_micro"),("左眼硬性渗出","L_hardex"),("左眼黄斑水肿","L_macular"),
           ("左眼新生血管","L_neov"),("右眼微血管瘤","R_micro"),("右眼硬性渗出","R_hardex"),
           ("右眼黄斑水肿","R_macular"),("右眼新生血管","R_neov")]
for src, dst in les_map:
    out[dst] = coh[src].astype("float")

out.to_csv(os.path.join(OUTDIR, "dev_set.csv"), index=False, encoding="utf-8")
print("dev_set.csv written:", out.shape)

# ---- QC summary ----
print("\n=== QC ===")
print("age clipped NA:", int(out["age"].isna().sum()))
g = out[grad]
print("gradeable n:", len(g), "| any_dr(>=grade1) rate: %.1f%%" % (100*g["any_dr"].mean()))
print("dr_raw rate in gradeable: %.1f%% (vs %.1f%% all-sample)" % (100*g["dr_raw"].mean(), 100*out["dr_raw"].mean()))
print("pdr n (gradeable):", int(g["pdr"].sum()), "| vtpr_dr n:", int(g["vtpr_dr"].sum()))
lip_cols = [d for _, d in lip_map]
for c in lip_cols:
    print("  %s missing %.1f%%" % (c, 100*out[c].isna().mean()))
all6 = out[lip_cols].notna().all(axis=1)
print("6-lipid complete n: %d (%.1f%%)" % (all6.sum(), 100*all6.mean()))
# modeling subset: gradeable + no NA in age/sex/any_dr/6 lipids
model_ok = grad & out[["age","male","any_dr"]+lip_cols].notna().all(axis=1)
print("full model listwise-complete n: %d (any_dr=%d)" % (model_ok.sum(), int(out.loc[model_ok,"any_dr"].sum())))

# ---- also copy NHANES B to ASCII working dir for stage 4 ----
import shutil
shutil.copy(NHANES_B, os.path.join(EXTDIR, "NHANES_B.csv"))
print("\nNHANES_B.csv copied to", EXTDIR)
