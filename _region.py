#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 GWAS 汇总统计中提取指定区域的全部 SNP (供共定位用)。

用法:
  python _region.py --gwas data/raw/DED_MVP_GCST90475895.h.tsv.gz \
                    --chr 2 --start 126100000 --end 128100000 \
                    --out data/processed/region_DED_chr2_127Mb.tsv
"""
import gzip, csv, math, argparse, sys, os
from _mr_prep import iter_rows, get_beta_se, fnum


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gwas', required=True)
    ap.add_argument('--chr', required=True)
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--end', type=int, required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    ncase = nctrl = None
    out = []
    for ci, rdr in iter_rows(a.gwas):
        need = [ci['chr'], ci['pos'], ci['p']]
        if any(k not in ci for k in ('chr', 'pos')):
            sys.exit("该 GWAS 无 chr/pos 列, 无法按区域提取")
        for row in rdr:
            if len(row) <= max(ci.values()):
                continue
            if str(row[ci['chr']]).strip() != str(a.chr):
                continue
            try:
                pos = int(float(row[ci['pos']]))
            except Exception:
                continue
            if pos < a.start or pos > a.end:
                continue
            beta, se = get_beta_se(row, ci)
            if beta is None or se is None or se <= 0:
                continue
            ea = (row[ci['ea']] if 'ea' in ci else '').strip().upper()
            oa = (row[ci['oa']] if 'oa' in ci else '').strip().upper()
            if {ea, oa} in ({'A', 'T'}, {'G', 'C'}):
                continue          # 回文 SNP 无法判定方向
            if 'ncase' in ci and ncase is None:
                ncase = fnum(row[ci['ncase']])
            if 'nctrl' in ci and nctrl is None:
                nctrl = fnum(row[ci['nctrl']])
            eaf = fnum(row[ci['eaf']]) if 'eaf' in ci else None
            out.append([row[ci['snp']].strip(), pos, ea, oa, beta, se,
                        fnum(row[ci['p']]), eaf])

    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['SNP', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'eaf'])
        w.writerows(out)
    print("区域 chr%s:%d-%d -> %d 个 SNP (cases=%s controls=%s) -> %s"
          % (a.chr, a.start, a.end, len(out), ncase, nctrl, a.out))


if __name__ == '__main__':
    main()
