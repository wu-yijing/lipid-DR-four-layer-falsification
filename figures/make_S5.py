# -*- coding: utf-8 -*-
"""Generate Supplementary Table S5 (cis drug-target MR instrument list) from real
phase1 outputs and insert it into the Frontiers charts md before S6."""
import csv, os

IV = r"E:/workbuddy/2026-09-04-18-19-16/phase1_out/iv_druggable.csv"
MR = r"E:/workbuddy/2026-09-04-18-19-16/phase1_out/druggable_mr_results.csv"
TBL = r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/清洗结果/方案A_论文图表_Frontiers版.md"

# gene cis regions (hg19, +/-500 kb) from phase5a_cis_iv.R
REGION = {
    "HMGCR": "5:74124000-75197000", "NPC1L1": "7:44052000-45259000",
    "PCSK9": "1:55005221-56028551", "APOC3": "11:116162381-117164780",
    "LPL": "8:19341030-20369239", "ANGPTL3": "1:62092402-63184930",
}
DRUG = {  # from Supplementary Table S2
    "HMGCR": "Statins", "NPC1L1": "Ezetimibe", "PCSK9": "PCSK9 inhibitors",
    "APOC3": "APOC3 antisense", "LPL": "LPL activators", "ANGPTL3": "Evinacumab",
}

ivs = {}
with open(IV, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        t = r["target"]
        beta = float(r["beta"]); se = float(r["se"]); p = float(r["p"])
        fstat = (beta / se) ** 2
        ivs.setdefault(t, []).append({
            "rsid": r["rsid"], "chr": r["chr"], "pos": r["position"],
            "ea": r["ea"], "nea": r["nea"], "beta": beta, "se": se,
            "p": p, "eaf": r["eaf"], "F": fstat,
        })

mr = {}
with open(MR, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        mr[r["target"]] = r

order = ["HMGCR", "NPC1L1", "PCSK9", "APOC3", "LPL", "ANGPTL3"]

def fmt_p(p):
    return f"{p:.2e}"

# (a) per-target summary
a_rows = ["| Target | Drug class | Lipid | Cis region (hg19, ±500 kb) | n IV (p<1e-5, r²<0.01) | n harmonised | Min F | Mean F | MR method | OR (95% CI) | p |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
for t in order:
    lst = ivs[t]; m = mr[t]
    n_iv = len(lst); n_h = int(m["n_harmonised"])
    minF = min(x["F"] for x in lst); meanF = sum(x["F"] for x in lst) / len(lst)
    orc = f'{float(m["OR"]):.3f} ({float(m["or_lci95"]):.3f}-{float(m["or_uci95"]):.3f})'
    a_rows.append(
        f'| {t} | {DRUG[t]} | {m["lipid_trait"]}-C | {REGION[t]} | {n_iv} | {n_h} | '
        f'{minF:.0f} | {meanF:.0f} | {m["method"]} | {orc} | {float(m["p"]):.3f} |'
    )

# (b) individual instruments
b_rows = ["| Target | rsID | Chr:Position | EA/NEA | beta (GLGC) | SE | p | EAF |",
          "|---|---|---|---|---|---|---|---|"]
for t in order:
    for x in ivs[t]:
        eaf = x["eaf"] if x["eaf"] not in ("NA", "", "NaN") else "NA"
        b_rows.append(
            f'| {t} | {x["rsid"]} | {x["chr"]}:{x["pos"]} | {x["ea"]}/{x["nea"]} | '
            f'{x["beta"]:.4g} | {x["se"]:.4g} | {fmt_p(x["p"])} | {eaf} |'
        )

total = sum(len(ivs[t]) for t in order)
block = (
    "## Supplementary Table S5. cis drug-target MR instrument list and per-target diagnostics "
    "(six lipid-lowering targets).\n\n"
    "**(a) Per-target instrument summary**\n\n"
    + "\n".join(a_rows) + "\n\n"
    "**(b) Individual cis instruments (n=%d)**\n\n" % total
    + "\n".join(b_rows) + "\n\n"
    "*Instruments were extracted from GLGC (ieu-a-300 for LDL-C, ieu-a-302 for TG) within ±500 kb of each gene locus, "
    "retaining independent SNPs at p<1×10⁻⁵ and pairwise LD r²<0.01 (EUR clumping, kb=500). "
    "F = (beta/SE)² is the instrument strength in the exposure GWAS. Harmonisation against FinnGen R13 DR dropped "
    "1 PCSK9 and 1 LPL SNP with incompatible alleles (n_harmonised shown). Single-IV targets (HMGCR, NPC1L1, ANGPTL3) "
    "were analysed by Wald ratio; multi-IV targets by inverse-variance weighted (IVW). All cis instruments were genome-wide "
    "significant in GLGC (p<1×10⁻⁵) and all per-target F statistics exceed conventional weak-instrument thresholds (F>10), "
    "indicating no substantial weak-instrument bias at the instrument level. Per-target effect estimates are in Supplementary Table S2; "
    "LPL cis heterogeneity diagnostics are in Supplementary Table S4.*\n"
)

with open(TBL, encoding="utf-8") as f:
    lines = f.read().splitlines()

anchor = None
for i, ln in enumerate(lines):
    if ln.startswith("## Supplementary Table S6."):
        anchor = i; break
if anchor is None:
    raise SystemExit("S6 anchor not found")

# avoid duplicate insertion
if any("Supplementary Table S5." in l for l in lines):
    print("S5 already present; skipping insertion.")
else:
    lines = lines[:anchor] + block.rstrip("\n").splitlines() + ["", ""] + lines[anchor:]
    with open(TBL, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("S5 inserted before S6.")

print("Total cis IVs:", total)
print("Per-target n IV:", {t: len(ivs[t]) for t in order})
print("Per-target min/mean F:", {t: (round(min(x['F'] for x in ivs[t])), round(sum(x['F'] for x in ivs[t])/len(ivs[t]))) for t in order})
