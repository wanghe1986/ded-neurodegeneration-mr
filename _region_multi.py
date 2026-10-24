#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
一次扫描 GWAS, 同时提取多个区域 (供共定位用)。大文件只解压一遍。
用法:
  python _region_multi.py --gwas data/raw/DED_xxx.gz --prefix data/processed/reg_DED \
      --regions "chr10:21035798-22035798,chr2:59302564-61302564"
输出: <prefix>_chr10_21035798.tsv 等
"""
import gzip, csv, argparse, sys, os
from _mr_prep import iter_rows, get_beta_se, fnum


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gwas', required=True)
    ap.add_argument('--regions', required=True,
                    help='"chr:start-end,chr:start-end" 窗口需事先算好')
    ap.add_argument('--prefix', required=True)
    a = ap.parse_args()

    regs = []
    for r in a.regions.split(','):
        c, se = r.split(':')
        s, e = se.split('-')
        regs.append((c.strip(), int(s), int(e), '%s_%s_%s' % (c.strip(), s, e)))
    print("待提取区域: %d 个" % len(regs))

    out = {r[3]: [] for r in regs}
    rmap = {}
    for r in regs:
        rmap.setdefault(r[0], []).append(r)
    ncase = nctrl = None

    for ci, rdr in iter_rows(a.gwas):
        if 'chr' not in ci or 'pos' not in ci:
            sys.exit("该 GWAS 无 chr/pos 列")
        for row in rdr:
            if len(row) <= max(ci.values()):
                continue
            ch = str(row[ci['chr']]).strip().lstrip('#')
            if ch not in rmap:
                continue
            try:
                pos = int(float(row[ci['pos']]))
            except Exception:
                continue
            for (rc, rs_, re_, tag) in rmap[ch]:
                if pos < rs_ or pos > re_:
                    continue
                beta, se = get_beta_se(row, ci)
                if beta is None or se is None or se <= 0:
                    continue
                ea = (row[ci['ea']] if 'ea' in ci else '').strip().upper()
                oa = (row[ci['oa']] if 'oa' in ci else '').strip().upper()
                if {ea, oa} in ({'A', 'T'}, {'G', 'C'}):
                    continue
                if 'ncase' in ci and ncase is None:
                    ncase = fnum(row[ci['ncase']])
                if 'nctrl' in ci and nctrl is None:
                    nctrl = fnum(row[ci['nctrl']])
                eaf = fnum(row[ci['eaf']]) if 'eaf' in ci else None
                out[tag].append([row[ci['snp']].strip(), pos, ea, oa, beta, se,
                                 fnum(row[ci['p']]), eaf])
                break

    os.makedirs(os.path.dirname(a.prefix) or '.', exist_ok=True)
    for tag, rows in out.items():
        fn = '%s_%s.tsv' % (a.prefix, tag)
        with open(fn, 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['SNP', 'pos', 'effect_allele', 'other_allele',
                        'beta', 'se', 'pval', 'eaf'])
            w.writerows(rows)
        print("  %s: %d 个 SNP -> %s" % (tag, len(rows), fn))
    print("样本量: cases=%s controls=%s" % (ncase, nctrl))


if __name__ == '__main__':
    main()
