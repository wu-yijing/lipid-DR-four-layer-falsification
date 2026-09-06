# -*- coding: utf-8 -*-
"""Build 方案A_中文手稿_v3_FrontEndocrinol修改版.docx (v2 of builder)."""
import re, os
from docx import Document
from docx.shared import Pt, Mm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

MD = r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/清洗结果/方案A_中文手稿_v3_FrontEndocrinol修改版.md"
TBL = r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/清洗结果/方案A_论文图表_Frontiers版.md"
FIGDIR = r"E:/workbuddy/2026-09-04-18-19-16/figs_frontiers"
OUT = r"C:/Users/Administrator/Desktop/血脂，糖网等数据清洗/清洗结果/方案A_中文手稿_v3_FrontEndocrinol修改版.docx"

FIGS = {
    1: (os.path.join(FIGDIR, "Figure1_study_design.png"), 160),
    2: (os.path.join(FIGDIR, "Figure2_MR_forest.png"), 160),
    3: (os.path.join(FIGDIR, "Figure3_calibration.png"), 160),
    4: (os.path.join(FIGDIR, "Figure4_ROC.png"), 85),
}

def read_lines(p):
    with open(p, encoding="utf-8") as f:
        return f.read().splitlines()

def split_bold(s):
    parts = re.split(r"\*\*(.+?)\*\*", s)
    out = []
    for i, seg in enumerate(parts):
        seg = seg.replace("`", "")
        if seg:
            out.append((seg, i % 2 == 1))
    return out or [(" ", False)]

def parse_md_table_block(lines, i):
    rows = []
    while i < len(lines) and lines[i].strip().startswith("|"):
        cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
        if not all(re.fullmatch(r":?-{3,}:?", c) for c in cells):
            rows.append(cells)
        i += 1
    return rows, i

def parse_table_sections(tbl_md):
    lines = read_lines(tbl_md)
    sections, cur, i = {}, None, 0
    while i < len(lines):
        st = lines[i].strip()
        if st.startswith("## ") and "Table" in st:
            m = re.match(r"##\s+(?:Supplementary\s+)?Table\s+(S?\d)\.", st)
            if m:
                cur = {"title": st[3:].strip(), "elems": []}
                sections[m.group(1)] = cur
                i += 1; continue
        if cur is not None:
            if st.startswith("|"):
                rows, i2 = parse_md_table_block(lines, i)
                if rows:
                    cur["elems"].append(("table", rows))
                i = i2; continue
            if re.match(r"^\*\*\([abc]\)", st):
                cur["elems"].append(("sub", st.strip("*")))
                i += 1; continue
            if st.startswith("*") and st.endswith("*") and not st.startswith("**"):
                cur["elems"].append(("note", st.strip("*").strip()))
                i += 1; continue
            if st.startswith("## ") or st == "---":
                cur = None
        i += 1
    for k, sec in sections.items():
        merged = []
        for e in sec["elems"]:
            if e[0] == "table" and merged and merged[-1][0] == "table":
                merged[-1] = ("table", merged[-1][1] + e[1])
            else:
                merged.append(e)
        sec["elems"] = merged
    return sections

TSEC = parse_table_sections(TBL)
print("sections:", sorted(TSEC.keys()))
for k in sorted(TSEC):
    kinds = [e[0] for e in TSEC[k]["elems"]]
    nr = [len(e[1]) for e in TSEC[k]["elems"] if e[0] == "table"]
    print(" ", k, kinds, nr)

doc = Document()
stn = doc.styles["Normal"]
stn.font.name = "Times New Roman"; stn.font.size = Pt(10.5)
stn.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
stn.paragraph_format.space_after = Pt(4)

def para(text="", size=10.5, bold=False, align=None, color=None,
         space_after=4, east="SimSun"):
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(space_after)
    if align is not None:
        p.alignment = align
    for seg, b in split_bold(text):
        r = p.add_run(seg)
        r.font.size = Pt(size); r.font.bold = bold or b
        r.font.name = "Times New Roman"
        r.element.rPr.rFonts.set(qn("w:eastAsia"), east)
        if color:
            r.font.color.rgb = RGBColor(*color)
    return p

def heading(text, level):
    h = doc.add_heading("", level=level)
    h.paragraph_format.space_before = Pt(10 if level <= 2 else 6)
    h.paragraph_format.space_after = Pt(4)
    clean = re.sub(r"\s*\[(?:Table|Figure|附表)[^\]]*\]", "", text)
    for seg, b in split_bold(clean):
        r = h.add_run(seg); r.font.bold = True
        r.font.color.rgb = RGBColor(0, 0, 0)
        r.font.name = "Times New Roman"
        r.font.size = Pt({1: 15, 2: 13, 3: 11.5}.get(level, 11))
        r.element.rPr.rFonts.set(qn("w:eastAsia"), "SimHei")

def set_three_line(table):
    tblPr = table._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge, sz, val in (("top", 12, "single"), ("bottom", 12, "single"),
                          ("left", 0, "none"), ("right", 0, "none"),
                          ("insideH", 0, "none"), ("insideV", 0, "none")):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), val); el.set(qn("w:sz"), str(sz))
        el.set(qn("w:color"), "000000")
        borders.append(el)
    tblPr.append(borders)
    for cell in table.rows[0].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tb = OxmlElement("w:tcBorders")
        bot = OxmlElement("w:bottom")
        bot.set(qn("w:val"), "single"); bot.set(qn("w:sz"), "6")
        bot.set(qn("w:color"), "000000")
        tb.append(bot); tcPr.append(tb)

def add_md_table(rows, fontsize=8.5):
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for ri, row in enumerate(rows):
        for ci in range(ncol):
            cell = t.cell(ri, ci)
            pcell = cell.paragraphs[0]
            pcell.paragraph_format.space_after = Pt(1)
            text = row[ci] if ci < len(row) else ""
            for seg, b in split_bold(text):
                r = pcell.add_run(seg)
                r.font.size = Pt(fontsize)
                r.font.bold = (ri == 0) or b
                r.font.name = "Times New Roman"
                r.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    set_three_line(t)
    return t

def add_caption(text, size=9, align=WD_ALIGN_PARAGRAPH.JUSTIFY):
    p = doc.add_paragraph()
    p.alignment = align
    p.paragraph_format.space_after = Pt(8)
    for seg, b in split_bold(text):
        r = p.add_run(seg)
        r.font.size = Pt(size); r.font.bold = b
        r.font.name = "Times New Roman"
        r.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    return p

def add_table_section(key, fontsize=8.5):
    sec = TSEC[key]
    add_caption("**" + sec["title"] + "**", size=9.5)
    for kind, val in sec["elems"]:
        if kind == "sub":
            para("**" + val + "**", size=8.5, space_after=2)
        elif kind == "table":
            add_md_table(val, fontsize)
            para("", size=4, space_after=2)
        elif kind == "note":
            add_caption("*" + val + "*", size=8.5)

FIG_LEGENDS = {
1: "**Figure 1. Four-layer falsification framework for lipid-based diabetic retinopathy risk stratification.** The framework tests three premises required for a lipid-based bedside risk tool: (a) causal support, (b) stable and interpretable discrimination, and (c) transferability. Layer 1: two-sample Mendelian randomization of six lipids (GLGC/UK Biobank genome-wide instruments) against FinnGen R13 diabetic retinopathy (DR; 15,353 cases/88,091 controls), with multivariable MR (MVMR), MR-PRESSO and cis drug-target MR. Layer 2: logistic and LASSO modelling in a Chinese EMR-based screening cohort (manually graded fundus photographs; gradeable n=1,966, 275 events), with bootstrap internal validation and multiple imputation (MI) as the primary analysis. Layer 3: external validation of Chinese coefficients in the NHANES 2005-2008 diabetes subgroup (n=1,119), unweighted and MEC-weighted. Layer 4: survey-weighted full reference model with 10-fold out-of-fold prediction and bi-directional temporal validation. Two generalisable failure modes were quantified: (i) complete-case analysis inflated discrimination relative to MI (dAUC 0.03-0.12), with missingness itself strongly DR-associated (OR 2.31); and (ii) MR causal signals were instrument-dependent (genome-wide TG OR 0.86 vs LPL-cis OR 1.267 - discordant and neither robust). Positive replacement factors were glycaemic control, diabetes duration and treatment intensity.",
2: "**Figure 2. Mendelian randomization estimates for genetically predicted lipids on diabetic retinopathy.** (a) Univariable inverse-variance weighted (IVW) estimates for six lipids (dark blue) and mutually adjusted multivariable MR estimates for TG, LDL-C and HDL-C (light blue). (b) cis drug-target MR estimates for six lipid-lowering targets (HMGCR, NPC1L1, PCSK9 for LDL-C; APOC3, LPL, ANGPTL3 for TG). Squares indicate ORs, horizontal lines the 95% CIs (log scale); the dashed line marks the null (OR=1). The univariable TG signal (OR 0.86, p=0.009) is shown alongside its attenuation after MVMR (0.94) and the opposite-direction (but non-significant and internally homogeneous, I2=0%) LPL-cis estimate (1.267). *p<0.05.",
3: "**Figure 3. External calibration of the transferred Chinese lipid scores in the NHANES diabetes subgroup.** Decile-wise mean predicted risk versus observed any-DR rate for the (a) core (age, sex, TC, HDL-C), (b) std (+TG, LDL-C) and (c) LASSO (+ApoB) versions. The dashed diagonal denotes perfect calibration. All three versions cluster far from the diagonal with unweighted calibration slopes of 0.31, 0.38 and -0.12 respectively; MEC-weighted calibration slopes were 0.79, 0.94 and -0.17 (Table 4), indicating the apparent severe miscalibration is partly a survey-design artefact while discrimination remains at chance under both weightings. The LASSO version (n=190) is exploratory.",
4: "**Figure 4. Discrimination of the lipid-only score versus the full reference model (NHANES diabetes subgroup, n=888).** Receiver operating characteristic curves for the survey-weighted 10-fold out-of-fold full reference model (AUC 0.735) and the Chinese-coefficient lipid-only score (AUC 0.505). The lipid-only score performs at chance level, whereas the full model - driven by HbA1c, diabetes duration and insulin use - shows robust discrimination (weighted OOF AUC 0.749; bi-directional temporal validation 0.748/0.722).",
}

def add_figure(n):
    path, w = FIGS[n]
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(4)
    p.add_run().add_picture(path, width=Mm(w))
    add_caption(FIG_LEGENDS[n], size=9)

lines = read_lines(MD)
i = 0
in_figlist = False
cur_section = ""
pending_table3 = False
pending_fig2 = False
done = set()
while i < len(lines):
    ln = lines[i].rstrip()
    st = ln.strip()
    if st == "---" or st == "":
        i += 1; continue
    if st.startswith("# ") and not st.startswith("## "):
        heading(st[2:], 1); i += 1; continue
    if st.startswith("## "):
        title = st[3:]
        if title.startswith("图表清单"):
            in_figlist = True; i += 1; continue
        in_figlist = False
        cur_section = title
        if title.startswith("参考文献"):
            heading(title, 2); i += 1
            while i < len(lines) and lines[i].strip():
                para(lines[i].strip(), size=9.5, space_after=2)
                i += 1
            heading("补充材料（Supplementary Material）", 2)
            para("以下补充表在投稿时作为 Separate files 上传；编号与正文引用一致。", size=9)
            for k in ["S2", "S3", "S4", "S5", "S6"]:
                add_table_section(k)
            supp_done = True
            continue
        heading(title, 2); i += 1; continue
    if st.startswith("### "):
        cur_section = st[4:]
        pending_table3 = cur_section.startswith("3.2")
        pending_fig2 = cur_section.startswith("3.1")
        heading(cur_section, 3); i += 1; continue
    if in_figlist:
        i += 1; continue
    if st.startswith(">"):
        para(st.lstrip("> ").strip(), size=9, color=(90, 90, 90), space_after=2)
        i += 1; continue
    if st.startswith("|"):
        rows, i2 = parse_md_table_block(lines, i)
        if "3.3" in cur_section and "4" not in done:
            add_table_section("4"); add_figure(3); done.add("4")
        elif "3.4" in cur_section and "5" not in done:
            add_table_section("5"); add_figure(4); done.add("5")
        elif "3.5" in cur_section and "S1" not in done:
            add_table_section("S1"); done.add("S1")
        i = i2; continue

    para(st)
    if pending_table3:
        add_table_section("3")
        pending_table3 = False
    if pending_fig2:
        add_figure(2)
        add_table_section("2")
        pending_fig2 = False
    if "[Figure 1]" in st:
        add_figure(1)
        add_table_section("1")
    i += 1

doc.save(OUT)
print("saved:", OUT, os.path.getsize(OUT), "bytes")
