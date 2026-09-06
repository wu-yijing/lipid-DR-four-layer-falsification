# Assemble paper-ready Tables 1-5 (markdown) from analysis CSVs
import pandas as pd, numpy as np

W = "E:/workbuddy/2026-09-04-18-19-16"
DST = "C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/清洗结果/方案A_论文表格_Table1-5.md"

def fp(p):
    if pd.isna(p): return "—"
    return "<0.001" if p < 0.001 else f"{p:.3f}"

L = []
A = L.append

A("# 方案 A 论文表格（Table 1–5）")
A("")
A("> 由分析产物 CSV 自动汇编，2026-09-05。所有数值与各阶段报告一致，可直接用于投稿正文。")
A("")

# ---------- Table 1 ----------
A("---")
A("")
A("## Table 1. Four-layer falsification of lipid-based bedside DR risk scores")
A("")
A("| Layer | Design | Data | Key result | Verdict |")
A("|---|---|---|---|---|")
A("| 1 Causal | Two-sample + multivariable MR | 6 lipid GWAS → FinnGen R13 DR (15,353/88,091) | TG univariable OR 0.86 (p=0.009) → MVMR OR 0.94 (p=0.24); no lipid independent | ✗ No causal support |")
A("| 2 Discrimination | Chinese screening cohort, logistic + LASSO, bootstrap 1,000 | n=1,207 (121 events) | AUC 0.669 (optimism-corrected); IDI 0.032; ApoB/LDL-C collinear (coefs +5.03 / −1.25) | △ Weak and uninterpretable |")
A("| 3 Transportability | Chinese coefficients → NHANES diabetes subgroup | n=1,041 / 494 / 190 | AUC 0.513 / 0.530 / 0.477; calibration slope 0.31 / 0.38 / −0.12 | ✗ Fails to generalize |")
A("| 4 Benchmark | Full model, survey-weighted, 10-fold oof (+ temporal validation) | n=888 (267 events) | AUC 0.735 (temporal 0.748/0.722); lipid-only 0.505; ΔAUC 0.230 | ✓ HbA1c, duration, insulin robust |")
A("")

# ---------- Table 2 ----------
A("---")
A("")
A("## Table 2. Mendelian randomization: univariable (IVW / MR-Egger / weighted median) and multivariable estimates")
A("")
m = pd.read_csv(f"{W}/phase1_out/mr_main.csv")
A("| Lipid | Method | n SNP | OR (95% CI) | p |")
A("|---|---|---|---|---|")
order = ["TG", "LDL", "HDL", "TC", "ApoA1", "ApoB"]
meth_short = {"Inverse variance weighted": "IVW", "MR Egger": "MR-Egger", "Weighted median": "Weighted median"}
for lip in order:
    sub = m[m.exposure == lip]
    for _, r in sub.iterrows():
        ms = meth_short.get(r["method"], r["method"])
        ci = f"{r['or']:.2f} ({r.or_lci95:.2f}–{r.or_uci95:.2f})"
        star = " **" if r.pval < 0.05 else ""
        end = "**" if r.pval < 0.05 else ""
        A(f"| {lip if ms=='IVW' else ''} | {ms} | {r.nsnp} | {star}{ci}{end} | {fp(r.pval)} |")
mv = pd.read_csv(f"{W}/phase1_out/mvmr_core_TG_LDL_HDL.csv")
A("| **Multivariable MR (TG+LDL-C+HDL-C mutually adjusted)** | | | | |")
for _, r in mv.iterrows():
    short = {"Triglycerides": "TG", "LDL cholesterol": "LDL-C", "HDL cholesterol": "HDL-C"}
    key = [v for k, v in short.items() if k in r.exposure][0]
    A(f"| {key} | MVMR | {r.nsnp} | {r['or']:.2f} ({r.or_lci95:.2f}–{r.or_uci95:.2f}) | {fp(r.pval)} |")
ps = pd.read_csv(f"{W}/phase1_out/presso10k_summary.csv")
A("| **MR-PRESSO (Nb = 10,000, outlier-corrected)** | | | | |")
for _, r in ps.iterrows():
    star = " **" if r.corrected_p < 0.05 else ""
    end = "**" if r.corrected_p < 0.05 else ""
    A(f"| {r.exposure} | PRESSO | {r.n_snp} | {star}{r.corrected_or:.2f}{end} | {fp(r.corrected_p)} |")
A("")
A("IVs: p<5e-8, LD clumped r²<0.001 (10,000 kb); palindromic SNPs harmonized by allele frequency. Outcome: FinnGen R13 DM_RETINOPATHY. Lipids per 1-SD genetically predicted level (GLGC) or per allele (UKB ApoA1/ApoB). All Global Tests P < 1e-4 (residual heterogeneity); corrected estimates agree with IVW — the univariable TG signal is not outlier-driven but is not independent of other lipids (see MVMR).")
A("")

# ---------- Table 3 ----------
A("---")
A("")
A("## Table 3. Chinese cohort model (n=1,207, 121 events): coefficients and collinearity diagnosis")
A("")
fc = pd.read_csv(f"{W}/phase2_out/phase2_model_fit.csv")
ls = pd.read_csv(f"{W}/phase2_out/phase2_lasso.csv").set_index("var")["coef"]
A("| Variable | Full-model coef | OR | OR per 1-SD | p | LASSO coef |")
A("|---|---|---|---|---|---|")
lab = {"age10": "Age (per 10 y)", "male": "Male sex", "tg_log": "log2 TG", "tc": "TC (mmol/L)",
       "hdl_c": "HDL-C (mmol/L)", "ldl_c": "LDL-C (mmol/L)", "apo_a1": "ApoA-I (g/L)", "apo_b": "ApoB (g/L)"}
for _, r in fc.iterrows():
    v = r.variable
    if v == "Intercept":
        A(f"| Intercept | {r.coef:.4f} | {r.OR:.3f} | — | {fp(r.p)} | {ls.get('(Intercept)', np.nan):.4f} |")
        continue
    psd = f"{r.OR_perSD:.3f}" if not pd.isna(r.OR_perSD) and v not in ("male",) else "—"
    lc = ls.get(v, np.nan)
    A(f"| {lab.get(v, v)} | {r.coef:.4f} | {r.OR:.3f} | {psd} | {fp(r.p)} | {lc:.4f} |")
A("")
A("Performance: apparent AUC 0.692; optimism-corrected AUC 0.669; Brier 0.0876; calibration slope 1.000 / intercept −0.000; lipid IDI 0.032.")
A("")

# ---------- Table 4 ----------
A("---")
A("")
A("## Table 4. External validation on NHANES diabetes subgroup (any DR): Chinese coefficients applied")
A("")
e = pd.read_csv(f"{W}/phase4_out/phase4_external_auc.csv")
A("| Model version | Variables beyond age/sex | n | Events | AUC (95% CI) | Brier | Calibration slope | Calibration intercept | Chinese-cohort AUC |")
A("|---|---|---|---|---|---|---|---|---|")
vd = {"core": "TC + HDL-C", "std": "+ TG + LDL-C", "lasso": "+ TG + LDL-C + ApoB (phase-2 LASSO)"}
for _, r in e.iterrows():
    cn = "—" if pd.isna(r.auc_chinese_dev) else f"{r.auc_chinese_dev:.3f}"
    A(f"| {r.model} | {vd[r.model]} | {r.n} | {r.events} | {r.auc:.3f} ({r.auc_lo:.3f}–{r.auc_hi:.3f}) | {r.brier:.3f} | {r.calib_slope:.3f} | {r.calib_intercept:.3f} | {cn} |")
A("")
A("NHANES lipids converted mg/dL → mmol/L (TC/HDL/LDL ÷38.67; TG ÷88.57) and ApoB → g/L (÷100). AUC CIs by bootstrap (1,000).")
A("")

# ---------- Table 5 ----------
A("---")
A("")
A("## Table 5. Full reference model, NHANES diabetes subgroup (n=888, 267 events; survey-weighted logistic)")
A("")
rf = pd.read_csv(f"{W}/phase4_out/phase4_reference_model.csv")
lab5 = {"age10": "Age (per 10 y)", "male": "Male sex", "hba1c_pct": "HbA1c (per 1 %)",
        "dm_duration": "Diabetes duration (per 1 y)", "egfr_ckdepi2021": "eGFR (CKD-EPI 2021, per mL/min/1.73m²)",
        "bmi": "BMI (per kg/m²)", "smoker_100": "Smoked ≥100 cigarettes", "insulin_use": "Insulin use",
        "oral_med": "Oral glucose-lowering medication", "sys_bp": "Systolic BP (per mmHg)"}
A("| Variable | OR | p |")
A("|---|---|---|")
for _, r in rf.iterrows():
    if r["var"] == "(Intercept)":
        continue
    star = " **" if r.p < 0.05 else ""
    end = "**" if r.p < 0.05 else ""
    A(f"| {lab5.get(r['var'], r['var'])} | {star}{r.OR:.2f}{end} | {fp(r.p)} |")
A("")
A("10-fold out-of-fold AUC = 0.735. Temporal external validation: train 2005-06 → test 2007-08 AUC 0.748; train 2007-08 → test 2005-06 AUC 0.722. Lipid-only comparator AUC 0.505 (ΔAUC = 0.230; NRI 0.410; IDI 0.151).")
A("")

open(DST, "w", encoding="utf-8").write("\n".join(L))
print("written:", DST, f"({len(L)} lines)")
