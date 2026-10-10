#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Disclosure gate: run this before every commit that leaves the author's machine.

The repository is intended to become public on acceptance, so a leak is permanent
even if it happens while the repository is still private. This gate answers one
question: **does anything about to be committed contain individual-level or
patient-identifying material?**

It checks three surfaces.

1. **Paths.** Files that must never be committed: record-level extracts, raw data
   directories, Office documents (which carry the manuscript), and statistical
   binaries. Several of these are also listed in ``.gitignore``; the gate
   re-checks, because an ignore rule added after a file was tracked does nothing.
2. **Content.** Long digit runs (a bare 9+ digit integer in an aggregate table is
   almost always an identifier), clinical Chinese vocabulary (source-report text),
   and personal-data field names.
3. **History** (``--history``). Whether such material was committed earlier. A
   file deleted in a later commit is still readable in the history, and becomes
   readable to everyone the moment the repository is made public.

Chinese *file paths* are expected in this repository (the author's local analysis
directory has Chinese in its name), so CJK is reported as a warning unless it
matches clinical vocabulary, which is fatal.

Usage
-----
    python tools/disclosure_gate.py                 # scan what a commit would include
    python tools/disclosure_gate.py --staged        # scan only the index
    python tools/disclosure_gate.py --history       # also scan the whole commit history
    python tools/disclosure_gate.py --history --strict    # history findings are fatal

Exit status is non-zero if any fatal finding is present.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# --------------------------------------------------------------------------
# never commit these
# --------------------------------------------------------------------------
FORBIDDEN_NAMES = [
    # record-level and raw extracts (also in .gitignore; listed because a file that
    # was tracked once stays tracked)
    "layer2_independent_reextraction.csv",
    "layer2_independent_reextraction_discordances.csv",
    "layer2_independent_reextraction_discordances_preadjudication.csv",
    "layer2_second_grader_sampling.csv",
    "layer2_text_mode_practice_key.csv",
    "layer2_text_mode_discordances.csv",
    "layer2_text_mode_5level_matrix__selftest_returnedxlsx.csv",
    "layer2_text_mode_agreement__selftest_returnedxlsx.json",
    "dev_set.csv",
    "dev_set_v2.csv",
    "dev_set_v2_dedup.csv",
    "severity_rebuild_audit.csv",
    "NHANES_B.csv",
]
FORBIDDEN_SUFFIXES = [".docx", ".doc", ".pdf", ".pptx", ".ppt", ".xlsx", ".xls",
                      ".sav", ".xpt", ".rds", ".dta", ".zip", ".7z", ".png", ".tif",
                      ".tiff", ".jpg", ".jpeg"]
FORBIDDEN_PATH_PARTS = ["nhanes_raw", "raw/", "/raw/", "data/", "ethics", "伦理",
                        "waiver", "_out/", "phase1_out/", "phase2_dev/", "phase2_ext/",
                        "phase2_out", "phase4_out", "phase6_out", "phase7_out"]

# --------------------------------------------------------------------------
# content patterns
# --------------------------------------------------------------------------
# a bare integer token of nine or more digits; decimals and exponents are excluded
#
# dbSNP rsIDs are excluded too: "rs112945592" is a variant identifier, not a record
# identifier, and its nine-digit numeric part would otherwise read as one. The two
# controls below (--self-test) keep that distinction falsifiable in both directions.
DIGITS = re.compile(r"(?<![\d.])(?<![rR][sS])[0-9]{9,}(?![\d.])")

#: Controls for DIGITS. Built from expressions rather than written as literals so that
#: this file -- which other scanners read as ordinary source -- does not itself contain
#: a bare nine-digit run.
DIGIT_CONTROLS = [
    ("screening record " + str(2 ** 34), True),    # a record number: must be caught
    ("rs" + str(2 ** 29), False),                  # a dbSNP rsID: must be ignored
]
CJK = re.compile(r"[\u4e00-\u9fff]+")
# clinical / personal vocabulary that has no business in an aggregate table or a script
CLINICAL = re.compile(
    r"患者|病历|住院|门诊|身份证|姓名|电话|地址|出生|随访|主诉|诊断|检查|报告|"
    r"糖尿病视网膜病变|增殖|非增殖|视网膜|眼底|糖化血红蛋白"
)
PERSONAL_FIELDS = re.compile(r"\b(patient_id|patid|pid|mrn|id_card|national_id|phone|address)\b",
                             re.I)

#: allowlist for CJK that is expected: the author's local analysis directory name,
#: which appears in the documented default paths.
CJK_ALLOW = re.compile(r"投稿资料|芬兰原始数据|糖尿病.*原始数据")

#: files that are scanned by content rules but must not be matched by them: this gate
#: states the forbidden vocabulary, and .gitignore explains in a comment which files
#: are protected. Neither is data.
SKIP_CONTENT = {"tools/disclosure_gate.py"}
COMMENT_ONLY = {".gitignore", ".gitattributes", ".gitmodules"}

#: genomic coordinates look like 11:116162381-117164780 and are not identifiers
GENOMIC = re.compile(r"[0-9XY]+:[0-9]{5,}")


def git(*args: str, check: bool = True) -> str:
    r = subprocess.run(["git", *args], cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    if check and r.returncode != 0:
        raise SystemExit("git %s failed: %s" % (" ".join(args), r.stderr.strip()))
    return r.stdout


def candidate_files(staged_only: bool) -> list[str]:
    out = git("diff", "--cached", "--name-only") if staged_only else \
        git("status", "--porcelain", "-uall")
    if staged_only:
        return [l for l in out.splitlines() if l.strip()]
    files = []
    for line in out.splitlines():
        if len(line) < 4:
            continue
        code, path = line[:2], line[3:].strip()
        if code == "!!":      # ignored and deliberately not scanned
            continue
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        files.append(path.strip('"'))
    return files


def read(path: str) -> str | None:
    p = ROOT / path
    if not p.is_file():
        return None
    try:
        return p.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return None


def scan() -> tuple[list[tuple], list[tuple]]:
    files = candidate_files(args.staged)
    fatal: list[tuple] = []
    warn: list[tuple] = []
    print("scanning %d path(s) about to be committed\n" % len(files))
    for f in files:
        base = Path(f).name
        if base in FORBIDDEN_NAMES:
            fatal.append(("path", f, "record-level extract must never be committed"))
            continue
        if Path(f).suffix.lower() in FORBIDDEN_SUFFIXES:
            fatal.append(("path", f, "binary or Office document"))
            continue
        if any(part in f for part in FORBIDDEN_PATH_PARTS):
            fatal.append(("path", f, "path matches a raw/intermediate-data pattern"))
            continue
        if f in SKIP_CONTENT:
            continue
        body = read(f)
        if body is None:
            continue
        if Path(f).name in COMMENT_ONLY:
            # the explanatory comments say which files are protected; that is not data
            body = "\n".join(l for l in body.splitlines() if not l.strip().startswith("#"))
        hits = [h for h in DIGITS.findall(body)
                if not GENOMIC.search(body[max(0, body.find(h) - 40):body.find(h) + 40])]
        if hits:
            fatal.append(("content", f, "%d bare digit run(s) of 9+ digits, e.g. %s"
                          % (len(hits), hits[:3])))
        for m in CLINICAL.finditer(body):
            ctx = body[max(0, m.start() - 30):m.end() + 30].replace("\n", " ")
            if CJK_ALLOW.search(ctx):
                continue
            fatal.append(("content", f, "clinical Chinese vocabulary: ...%s..." % ctx))
            break
        pf = PERSONAL_FIELDS.findall(body)
        if pf:
            warn.append(("content", f, "personal-data field name(s): %s" % sorted(set(pf))[:4]))
        cjk = [c for c in CJK.findall(body) if not CJK_ALLOW.search(c)]
        cjk = [c for c in cjk if not CLINICAL.search(c)]
        if cjk:
            warn.append(("content", f, "%d non-clinical CJK run(s), e.g. %s"
                         % (len(cjk), cjk[:3])))
    return fatal, warn


#: paths known to have been committed earlier and since removed. Listed so the gate
#: reports them as informational rather than alarming; all four are aggregate figures.
HISTORY_KNOWN = {
    "figures/Figure1_study_design.png": "aggregate study-design figure, removed from HEAD",
    "figures/Figure2_MR_forest.png": "aggregate MR forest plot, removed from HEAD",
    "figures/Figure3_calibration.png": "aggregate calibration plot, removed from HEAD",
    "figures/Figure4_ROC.png": "aggregate ROC plot, removed from HEAD",
}


def scan_committed() -> list[tuple]:
    """Binaries and Office documents that are already tracked. Not a breach by itself,
    but they are regenerable and a reviewer comparing the repository to the manuscript
    should know they are there."""
    out = []
    for f in git("ls-tree", "-r", "--name-only", "HEAD", check=False).splitlines():
        f = f.strip()
        if f and Path(f).suffix.lower() in FORBIDDEN_SUFFIXES:
            out.append(("committed", f, "tracked binary/Office file (pre-existing)"))
    return out


def scan_history(strict: bool) -> list[tuple]:
    """Report (a) any record-level filename ever committed and (b) files that once
    existed, are no longer in HEAD, and contain identifiers or clinical vocabulary.
    Chinese in scripts that are still present is not reported: it is reviewable."""
    findings: list[tuple] = []

    for f in FORBIDDEN_NAMES:
        out = git("log", "--all", "--pretty=format:%h %ad", "--date=short", "--", "**/" + f,
                  check=False)
        if out.strip():
            findings.append(("history", f, "record-level file committed at some point: %s"
                             % out.strip().splitlines()[0]))

    commits = [c for c in git("log", "--all", "--pretty=format:%h").split() if c]
    latest: dict[str, str] = {}
    for c in commits:
        for line in git("ls-tree", "-r", "--name-only", c).splitlines():
            if line.strip():
                latest[line.strip()] = c
    try:
        in_head = set(git("ls-tree", "-r", "--name-only", "HEAD").split())
    except Exception:
        in_head = set()

    print("\nhistory: %d commit(s), %d distinct path(s) ever present" % (len(commits), len(latest)))
    removed = sorted(p for p in latest if p not in in_head)
    print("history: %d path(s) no longer in HEAD" % len(removed))

    for p in removed:
        t = git("show", "%s:%s" % (latest[p], p))
        if not t:
            continue
        dd = [d for d in DIGITS.findall(t)
              if not GENOMIC.search(t[max(0, t.find(d) - 40):t.find(d) + 40])]
        cj = CLINICAL.search(t)
        note = HISTORY_KNOWN.get(p)
        if note:
            findings.append(("history-known", p, note))
            continue
        if dd:
            findings.append(("history", p, "removed file contains %d bare digit run(s), e.g. %s"
                             % (len(dd), dd[:3])))
        if cj and not CJK_ALLOW.search(t):
            findings.append(("history", p, "removed file contains clinical Chinese vocabulary: %s"
                             % cj.group(0)))
    return findings


def self_test() -> int:
    """Check that the digit-run detector still discriminates, in both directions.

    A gate that has stopped matching is worse than no gate: it reports clean while
    letting identifiers through. Run this before trusting a PASS.
    """
    print("digit-run detector controls")
    bad = 0
    for sample, must_match in DIGIT_CONTROLS:
        got = bool(DIGITS.search(sample))
        ok = got == must_match
        bad += 0 if ok else 1
        print("  %-4s %-30s expected=%-7s got=%s"
              % ("ok" if ok else "FAIL", sample,
                 "caught" if must_match else "ignored", "caught" if got else "ignored"))
    print("\nSELF-TEST: %s" % ("PASS" if bad == 0 else "FAIL (%d control(s))" % bad))
    return 1 if bad else 0


def main(argv=None) -> int:
    global args
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--staged", action="store_true", help="scan only the index")
    ap.add_argument("--history", action="store_true", help="also scan the full commit history")
    ap.add_argument("--strict", action="store_true", help="treat history findings as fatal")
    ap.add_argument("--self-test", action="store_true",
                    help="check the digit-run detector, then exit")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()

    fatal, warn = scan()

    print("FATAL (%d)" % len(fatal))
    for kind, f, why in fatal:
        print("  [%s] %s\n        %s" % (kind, f, why))
    if not fatal:
        print("  none")

    print("\nWARN (%d)" % len(warn))
    for kind, f, why in warn:
        print("  [%s] %s\n        %s" % (kind, f, why))
    if not warn:
        print("  none")

    pre = scan_committed()
    print("\nALREADY COMMITTED (%d)" % len(pre))
    for kind, f, why in pre:
        print("  [%s] %s\n        %s" % (kind, f, why))
    if not pre:
        print("  none")

    hist: list[tuple] = []
    if args.history:
        hist = scan_history(args.strict)
        print("\nHISTORY (%d)%s" % (len(hist), "  [fatal in --strict]" if args.strict else "  [advisory]"))
        for kind, f, why in hist:
            print("  [%s] %s\n        %s" % (kind, f, why))
        if not hist:
            print("  none")

    bad = len(fatal) + (len(hist) if args.strict else 0)
    print("\n" + "=" * 72)
    print("DISCLOSURE GATE: %s" % ("PASS" if bad == 0 else "BLOCKED (%d fatal)" % bad))
    print("=" * 72)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
