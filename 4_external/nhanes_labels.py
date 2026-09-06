#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""导出 OPXRET 与 DIQ 的关键变量标签 + 取值标签, 定位 DR 分级变量与糖尿病病程变量。"""
import os, pyreadstat

RAW = r"E:/workbuddy/2026-09-04-18-19-16/nhanes_raw_fixed"

def dump_labels(comp, patterns):
    path = os.path.join(RAW, f"{comp}.xpt")
    df, meta = pyreadstat.read_xport(path)
    vlabels = meta.variable_value_labels  # dict var -> {code: label}
    var_labels = meta.column_names_to_labels
    print(f"\n===== {comp} =====")
    for col in meta.column_names:
        full = var_labels.get(col, "")
        if any(p in col.upper() or p in full.upper() for p in patterns):
            print(f"\n[{col}] {full}")
            vl = vlabels.get(col)
            if vl:
                for code, lab in list(vl.items())[:25]:
                    print(f"    {code!r}: {lab}")

# DR 分级相关 / 视网膜病变 / 视力威胁
dump_labels("OPXRET_D", ["RET", "DR", "GRD", "GRADE", "VCDR", "MACUL", "EDEMA", "RPE", "EXU", "HEME", "MICRO"])
# 糖尿病确诊/病程/用药
dump_labels("DIQ_D", ["DIAB", "TOLD", "AGE", "INSUL", "PILL", "DURATION", "040", "050", "070", "080", "010"])
