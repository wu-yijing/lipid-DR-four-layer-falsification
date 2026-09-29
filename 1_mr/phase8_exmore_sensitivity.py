# -*- coding: utf-8 -*-
"""phase8 (MR): control-exclusion sensitivity analysis against FinnGen R13 DM_RETINOPATHY_EXMORE
   Rationale: the primary endpoint compares DR cases with controls not restricted to people without
   diabetes, so the estimate could capture an effect of lipids on type 2 diabetes (index-event bias).
   This script re-runs the six-lipid and the six cis drug-target analyses on the EXMORE endpoint.
   Inputs : phase1_out/harmonised_data.csv, phase1_out/iv_druggable.csv, <EXMORE>.gz
   Outputs: output_tables/MR_exmore_sensitivity.csv
   The same harmonisation rule as Table S2 is applied (palindromic variants with unavailable EAF dropped).
"""
import sys, os, csv, gzip, numpy as np, pandas as pd
from scipy import stats

P = r"E:/workbuddy/2026-09-04-18-19-16/phase1_out"
EXMORE = r"E:/workbuddy/BMC Genomics投稿资料/DR，DN，DM芬兰原始数据/finngen_R13_DM_RETINOPATHY_EXMORE.gz"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output_tables")
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
comp = lambda a: COMP.get(str(a).upper(), str(a).upper())

def ivw(b, se):
    w = 1 / se ** 2
    est = np.sum(w * b) / np.sum(w)
    se0 = np.sqrt(1 / np.sum(w))
    Q = np.sum(w * (b - est) ** 2); df = len(b) - 1
    return est, se0 * np.sqrt(max(1.0, Q / df)) if df > 0 else se0

def egger(b, se, bx):
    X = np.column_stack([np.ones(len(bx)), bx]); W = np.diag(1 / se ** 2)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ b)
    resid = b - X @ beta
    cov = np.linalg.inv(X.T @ W @ X) * ((resid.T @ W @ resid) / (len(b) - 2))
    return beta[1], np.sqrt(cov[1, 1]), beta[0], np.sqrt(cov[0, 0])

def wmed(b, se):
    w = 1 / se ** 2; o = np.argsort(b); b, w = b[o], w[o]
    p = np.cumsum(w) / np.sum(w) - 0.5 * w / np.sum(w)
    k = min(max(np.searchsorted(p, 0.5), 1), len(b) - 1)
    return b[k - 1] + (b[k] - b[k - 1]) * (0.5 - p[k - 1]) / (p[k] - p[k - 1])

def align(ea, nea, ref, alt):
    ea, nea, ref, alt = [str(x).upper() for x in (ea, nea, ref, alt)]
    if ea == alt and nea == ref: return +1
    if ea == ref and nea == alt: return -1
    if comp(ea) == alt and comp(nea) == ref: return +1
    if comp(ea) == ref and comp(nea) == alt: return -1
    return None

# ---- 1) SNP list + extraction from the EXMORE release -------------------------------
h = pd.read_csv(f"{P}/harmonised_data.csv", low_memory=False)
h = h[h['mr_keep'] == True]
iv = pd.read_csv(f"{P}/iv_druggable.csv")
iv.columns = ['target', 'trait', 'SNP', 'chr', 'pos', 'ea', 'nea', 'b_exp', 'se_exp', 'p_exp', 'eaf']
snps = set(h['SNP']) | set(iv['SNP'])
print("SNPs to extract:", len(snps), flush=True)
rows = []
CACHE1, CACHE2 = f"{P}/finngen_exmore_subset.tsv", f"{P}/finngen_exmore_cis_subset.tsv"
if os.path.exists(CACHE1) and os.path.exists(CACHE2):
    frames = []
    for path in (CACHE1, CACHE2):
        d = pd.read_csv(path, sep="\t").rename(columns={"rsid": "SNP"})
        frames.append(d[["SNP", "ref", "alt", "beta", "sebeta"]])
    exm = pd.concat(frames, ignore_index=True)
    exm[["beta", "sebeta"]] = exm[["beta", "sebeta"]].astype(float)
    exm = exm.drop_duplicates("SNP").set_index("SNP")
else:
    with gzip.open(EXMORE, 'rt', encoding='utf-8', errors='replace') as f:
        rd = csv.reader(f, delimiter='\t'); hdr = next(rd); i = {c: k for k, c in enumerate(hdr)}
        for r in rd:
            hit = [s for s in r[i['rsids']].split(',') if s in snps]
            if hit:
                rows.append((hit[0], r[i['ref']], r[i['alt']], r[i['beta']], r[i['sebeta']], r[i['pval']]))
                snps.discard(hit[0])
            if not snps: break
    exm = pd.DataFrame(rows, columns=['SNP', 'ref', 'alt', 'beta', 'sebeta', 'pval'])
    exm[['beta', 'sebeta']] = exm[['beta', 'sebeta']].astype(float)
    exm = exm.set_index('SNP')
print("extracted:", len(exm), "not found in the EXMORE release:", len(snps), flush=True)

# ---- 2) six lipids ------------------------------------------------------------------
out = []
for lbl in ['TG', 'LDL', 'HDL', 'TC', 'ApoA1', 'ApoB']:
    s = h[h['exposure'] == lbl].merge(exm, left_on='SNP', right_index=True, how='inner')
    b, se, bx = [], [], []
    for _, r in s.iterrows():
        sg = align(r['effect_allele.exposure'], r['other_allele.exposure'], r['ref'], r['alt'])
        if sg is None or r['beta.exposure'] == 0: continue
        b.append(sg * r['beta'] / r['beta.exposure']); se.append(r['sebeta'] / abs(r['beta.exposure']))
        bx.append(r['beta.exposure'])
    b, se, bx = np.array(b), np.array(se), np.array(bx)
    e, s_re = ivw(b, se)
    eg, egse, ic, icse = egger(b, se, bx)
    out.append(dict(analysis='six_lipids', term=lbl, n=len(b),
                    OR=np.exp(e), lo=np.exp(e - 1.96 * s_re), hi=np.exp(e + 1.96 * s_re),
                    p=2 * (1 - stats.norm.cdf(abs(e / s_re))),
                    egger_OR=np.exp(eg), egger_intercept_p=2 * (1 - stats.norm.cdf(abs(ic / icse))),
                    wmedian_OR=np.exp(wmed(b, se))))

# ---- 3) six cis drug targets (same harmonisation as Table S2) -----------------------
for tgt, gsub in iv.groupby('target'):
    m = gsub[gsub.eaf.notna()].merge(exm, left_on='SNP', right_index=True, how='inner', suffixes=('', '_o'))
    b, se = [], []
    for _, r in m.iterrows():
        sg = align(r['ea'], r['nea'], r['ref'], r['alt'])
        if sg is None or r['b_exp'] == 0: continue
        b.append(sg * r['beta'] / r['b_exp']); se.append(r['sebeta'] / abs(r['b_exp']))
    b, se = np.array(b), np.array(se)
    if len(b) == 1: e, s_re = b[0], se[0]
    else: e, s_re = ivw(b, se)
    out.append(dict(analysis='cis_drug_target', term=tgt, n=len(b),
                    OR=np.exp(e), lo=np.exp(e - 1.96 * s_re), hi=np.exp(e + 1.96 * s_re),
                    p=2 * (1 - stats.norm.cdf(abs(e / s_re))), egger_OR=np.nan,
                    egger_intercept_p=np.nan, wmedian_OR=np.nan))

res = pd.DataFrame(out)
os.makedirs(OUT, exist_ok=True)
res.to_csv(f"{OUT}/MR_exmore_sensitivity.csv", index=False)
print(res.to_string(index=False, float_format=lambda x: "%.4g" % x))
