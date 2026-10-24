#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用 eQTLGen 的等位基因频率表 (欧洲人群, 31,684) 给缺失 eaf 的 IV 文件补频率。
用途: IMSGC 2019 (MS) 的 harmonised 汇总统计没有 effect_allele_frequency 列,
      导致 Steiger 定向检验无法执行。这里用同源欧洲人群的频率替代。

用法: python _fill_maf.py --iv data/processed/iv_MS.tsv
"""
import gzip, csv, sys, argparse

MAF_FILE = "data/raw/eqtlgen_maf.txt.gz"   # SNP, hg19_chr, hg19_pos, AlleleA, AlleleB,
                                           # allA_total, allAB_total, allB_total, AlleleB_all

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--iv', required=True)
    ap.add_argument('--maf', default=MAF_FILE)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.iv), delimiter='\t'))
    hdr = list(rows[0].keys())
    if 'eaf' not in hdr:
        sys.exit("IV 文件无 eaf 列")
    need = {r['SNP'] for r in rows if not r.get('eaf')}
    if not need:
        print("无需补充, eaf 已完整"); return
    print("待补频率 SNP: %d / %d" % (len(need), len(rows)))

    fill = {}
    with gzip.open(a.maf, 'rt', errors='replace') as fh:
        rdr = csv.reader(fh, delimiter='\t')
        h = next(rdr)
        i = {c: k for k, c in enumerate(h)}
        for row in rdr:
            s = row[i['SNP']]
            if s in need:
                try:
                    fill[s] = (row[i['AlleleB']].upper(), float(row[i['AlleleB_all']]))
                except (ValueError, IndexError):
                    pass
    print("在 eQTLGen MAF 中命中: %d" % len(fill))

    n_ok = 0
    for r in rows:
        if r.get('eaf'): continue
        if r['SNP'] in fill:
            b, f = fill[r['SNP']]
            ea = r['effect_allele'].strip().upper()
            # AlleleB_all 是 AlleleB 的频率; 与 effect_allele 对齐
            r['eaf'] = '%.6g' % (f if b == ea else 1.0 - f)
            n_ok += 1
    with open(a.iv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter='\t'); w.writeheader(); w.writerows(rows)
    print("已补齐 %d 个 -> %s" % (n_ok, a.iv))
    miss = [r['SNP'] for r in rows if not r.get('eaf')]
    if miss: print("仍缺失 %d 个: %s" % (len(miss), ', '.join(miss[:10])))

if __name__ == '__main__':
    main()
