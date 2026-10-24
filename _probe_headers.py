# -*- coding: utf-8 -*-
"""探查三个数据文件的真实列名与结构"""
import gzip, os, sys

RAW = "data/raw"
FILES = {
    "eQTLGen 主文件": "eqtlgen_cis_fdr005.txt.gz",
    "MAF 文件":       "eqtlgen_maf.txt.gz",
    "AD GWAS":        "AD_Bellenguez2022.txt.gz",
}

for name, fn in FILES.items():
    path = os.path.join(RAW, fn)
    print("=" * 90)
    print(f"### {name}   ({fn})")
    if not os.path.exists(path):
        print("  文件不存在 (可能仍在移动中)")
        continue
    sz = os.path.getsize(path) / 1024**2
    print(f"  大小: {sz:.1f} MB")
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            header = f.readline().rstrip("\n")
            l1 = f.readline().rstrip("\n")
            l2 = f.readline().rstrip("\n")
        cols = header.split("\t") if "\t" in header else header.split(",")
        print(f"  分隔符: {'TAB' if chr(9) in header else ('逗号' if ',' in header else '未知')}")
        print(f"  列数: {len(cols)}")
        print(f"  列名: {cols}")
        v1 = l1.split("\t") if "\t" in l1 else l1.split(",")
        print(f"  首行: {v1[:12]}")
    except Exception as e:
        print(f"  读取失败: {e}")
    print()
