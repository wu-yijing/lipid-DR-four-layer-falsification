#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
阶段 0.3 - 解析各 XPT, 探查关键变量名与维度 (任务2)
读取 nhanes_raw_fixed/ 下真实 XPT, 打印每个组件的关键变量, 供后续合并/衍生参考。
"""
import os, pyreadstat, pandas as pd

RAW = r"E:/workbuddy/2026-09-04-18-19-16/nhanes_raw_fixed"

# 想重点确认的关键变量(部分名为 NHANES 标准命名, 实际以文件为准)
KEYWORDS = {
    "DEMO": ["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "WTINT2YR", "WTMEC2YR"],
    "TCHOL": ["SEQN", "LBXTC"],
    "HDL": ["SEQN", "LBDHDD"],
    "TRIGLY": ["SEQN", "LBXTR"],
    "APOB": ["SEQN", "LBXAPB"],
    "GHB": ["SEQN", "LBXGH"],
    "GLU": ["SEQN", "LBXGLU"],
    "DIQ": ["SEQN", "DIQ010", "DIQ040", "DIQ050", "DIQ070", "DIQ080"],
    "BIOPRO": ["SEQN", "LBXSCR", "LBXSAL", "URXUMA", "URXUCR"],
    "BMX": ["SEQN", "BMXBMI", "BMXWT", "BMXHT"],
    "BPX": ["SEQN", "BPXSY1", "BPXDI1"],
    "SMQ": ["SEQN", "SMQ020"],
    "OPXRET": ["SEQN", "VODRGRD", "VOSDRGRD", "VODRSOP", "VOSRSOP", "VODRDTY", "VOSRDTY"],
}

def main():
    comps = ["DEMO_D", "TCHOL_D", "HDL_D", "TRIGLY_D", "APOB_E", "GHB_D", "GLU_D",
             "DIQ_D", "BIOPRO_D", "BMX_D", "BPX_D", "SMQ_D", "OPXRET_D"]
    for comp in comps:
        path = os.path.join(RAW, f"{comp}.xpt")
        if not os.path.exists(path):
            print(f"\n##### {comp}: 文件缺失 #####")
            continue
        try:
            df, meta = pyreadstat.read_xport(path)
        except Exception as e:
            print(f"\n##### {comp}: 读取失败 {e} #####")
            continue
        print(f"\n##### {comp} #####")
        print(f"  维度: {df.shape[0]} 行 x {df.shape[1]} 列")
        # 显示关键词变量是否存在
        cols = list(df.columns)
        kws = KEYWORDS.get(comp.split('_')[0], [])
        found = [c for c in kws if c in cols]
        missing = [c for c in kws if c not in cols]
        if kws:
            print(f"  关键变量 命中: {found}")
            if missing:
                print(f"  关键变量 缺失: {missing}")
        # 列出所有列名(便于人工核对)
        print(f"  全部列({len(cols)}): {cols}")
        # 展示前 2 行关键变量
        show = [c for c in (kws if kws else cols[:8]) if c in cols]
        if show:
            print("  前2行样例:")
            print(df[show].head(2).to_string(index=False))

if __name__ == "__main__":
    main()
