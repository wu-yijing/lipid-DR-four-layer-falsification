#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
阶段 0.3 - 重新下载正确的 NHANES XPT 数据
原 raw/ 下 25 个文件全是 CDC "Page Not Found" 404 错误页(假 http=200)。
正确端点已定位: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/<YEAR>/DataFiles/<COMP>.xpt
注意: 眼底照相/DR 分级组件官方名为 OPXRET (非 RETINO)。
"""
import os, sys, urllib.request, time, hashlib

BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
OUT_DIR = r"E:/workbuddy/2026-09-04-18-19-16/nhanes_raw_fixed"
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"

# (year_dir, component, out_name)
# 2005-2006 -> year_dir=2005, suffix _D ; 2007-2008 -> year_dir=2007, suffix _E
FILES = [
    # 2005-2006 cycle (suffix _D)
    ("2005", "DEMO_D"), ("2005", "TCHOL_D"), ("2005", "HDL_D"), ("2005", "TRIGLY_D"),
    ("2005", "GHB_D"), ("2005", "GLU_D"), ("2005", "DIQ_D"), ("2005", "BIOPRO_D"),
    ("2005", "BMX_D"), ("2005", "BPX_D"), ("2005", "SMQ_D"), ("2005", "OPXRET_D"),
    # 2007-2008 cycle (suffix _E)
    ("2007", "DEMO_E"), ("2007", "TCHOL_E"), ("2007", "HDL_E"), ("2007", "TRIGLY_E"),
    ("2007", "GHB_E"), ("2007", "GLU_E"), ("2007", "DIQ_E"), ("2007", "BIOPRO_E"),
    ("2007", "BMX_E"), ("2007", "BPX_E"), ("2007", "SMQ_E"), ("2007", "OPXRET_E"),
    ("2007", "APOB_E"),
]

def xpt_valid(path):
    """真实 SAS XPT 以 'HEADER R' (或 'SAS'/'XPT') 起始; 错误页以 <!DOCTYP 起始。"""
    try:
        with open(path, "rb") as f:
            head = f.read(8)
        return head.startswith(b"HEADER R") or head.startswith(b"SAS     ") or head.startswith(b"XPORT")
    except Exception:
        return False

def download(url, dest, retries=3):
    last_err = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            with open(dest, "wb") as f:
                f.write(data)
            if xpt_valid(dest) and os.path.getsize(dest) > 100000:
                return True, os.path.getsize(dest), None
            # 下载到但内容不是有效 XPT (可能是错误页)
            last_err = f"invalid XPT magic (size={os.path.getsize(dest)})"
        except Exception as e:
            last_err = str(e)
        time.sleep(2)
    return False, 0, last_err

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    manifest = []
    ok = 0
    fail = 0
    for year, comp in FILES:
        url = f"{BASE}/{year}/DataFiles/{comp}.xpt"
        dest = os.path.join(OUT_DIR, f"{comp}.xpt")
        print(f"[fetch] {comp:10s} <- {url}", flush=True)
        success, size, err = download(url, dest)
        if success:
            ok += 1
            print(f"        OK  size={size:,} bytes", flush=True)
            manifest.append((comp, "OK", size, url))
        else:
            fail += 1
            print(f"        FAIL err={err}", flush=True)
            manifest.append((comp, "FAIL", 0, url))
    # write manifest
    with open(os.path.join(OUT_DIR, "_manifest.tsv"), "w", encoding="utf-8") as f:
        f.write("component\tstatus\tsize_bytes\turl\n")
        for comp, st, sz, url in manifest:
            f.write(f"{comp}\t{st}\t{sz}\t{url}\n")
    print(f"\n=== DONE: {ok} OK, {fail} FAIL (total {len(FILES)}) ===")
    print(f"manifest -> {os.path.join(OUT_DIR, '_manifest.tsv')}")
    if fail:
        sys.exit(1)

if __name__ == "__main__":
    main()
