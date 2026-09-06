#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""读取关键变量的 取值->标签 映射 (value_labels='raw' 保留数值编码)。"""
import os, pyreadstat

RAW = r"E:/workbuddy/2026-09-04-18-19-16/nhanes_raw_fixed"

def show(comp, vars_):
    path = os.path.join(RAW, f"{comp}.xpt")
    df, meta = pyreadstat.read_xport(path, value_labels='raw')
    vl = meta.variable_value_labels
    print(f"\n===== {comp} =====")
    for v in vars_:
        if v in vl:
            print(f"\n[{v}]")
            for code, lab in vl[v].items():
                print(f"    {code!r}: {lab}")
        else:
            print(f"\n[{v}] (无取值标签/变量缺失)")

show("OPXRET_D", ["OPDURET", "OPDURL4", "OPDUARMA", "OPDDRET", "OPDSRET", "OPDDRL4", "OPDSRL4"])
show("DIQ_D", ["DIQ010", "DID040", "DIQ050", "DIQ080", "DID070", "DIQ220"])
show("DEMO_D", ["RIAGENDR", "RIDRETH1"])
show("BIOPRO_D", ["LBXSCR"])
