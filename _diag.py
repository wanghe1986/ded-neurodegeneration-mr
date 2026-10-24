#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断: 统计 meta GWAS 的 P 值分布与可用统计量比例"""
import gzip, csv, sys, collections

PATH = sys.argv[1]
cnt = collections.Counter()
buckets = collections.Counter()
top = []

with gzip.open(PATH, 'rt', encoding='utf-8', errors='replace') as fh:
    rdr = csv.reader(fh, delimiter='\t')
    hdr = next(rdr)
    h = [c.strip().lstrip('#').lower() for c in hdr]
    i = {c: k for k, c in enumerate(h)}
    print("表头:", hdr)
    ip = i['p_value']; ise = i['standard_error']; ior = i['odds_ratio']
    ilo = i.get('ci_lower'); ihi = i.get('ci_upper')
    isnp = i['rsid']; ichr = i['chromosome']; ipos = i['base_pair_location']
    iea = i['effect_allele']; ioa = i['other_allele']

    def f(x):
        x = x.strip()
        if x in ('', 'NA', 'NaN', '.', 'None'): return None
        try: return float(x)
        except ValueError: return None

    for row in rdr:
        cnt['total'] += 1
        p = f(row[ip])
        if p is None:
            cnt['p_na'] += 1; continue
        # 可用统计量?
        se = f(row[ise]); lo = f(row[ilo]) if ilo is not None else None
        hi = f(row[ihi]) if ihi is not None else None
        usable = (se is not None and se > 0) or (lo and hi and lo > 0 and hi > 0)
        for thr, nm in ((5e-8, 'p5e8'), (5e-7, 'p5e7'), (5e-6, 'p5e6'), (1e-5, 'p1e5')):
            if p < thr:
                cnt[nm] += 1
                if usable: cnt[nm + '_usable'] += 1
        if p < 5e-6:
            ea = row[iea].strip().upper(); oa = row[ioa].strip().upper()
            pal = {ea, oa} in ({'A', 'T'}, {'G', 'C'})
            cnt['p5e6_palin'] += 1 if pal else 0
            if usable and not pal:
                buckets[(row[ichr], int(float(row[ipos])) // 10000000)] += 1
                top.append((p, row[isnp], row[ichr], row[ipos], ea, oa, f(row[ior])))

print()
for k in ('total', 'p_na', 'p5e8', 'p5e8_usable', 'p5e7', 'p5e7_usable',
          'p5e6', 'p5e6_usable', 'p5e6_palin', 'p1e5', 'p1e5_usable'):
    print("%-14s %d" % (k, cnt[k]))
print()
print("10 Mb 分箱 (P<5e-6 且统计量可用、非回文), 前 25:")
for (c, b), n in sorted(buckets.items(), key=lambda x: -x[1])[:25]:
    print("   chr%s bin %d-%d Mb : %d" % (c, b * 10, b * 10 + 10, n))
print()
print("Top 20 位点:")
for t in sorted(top)[:20]:
    print("   %-14s chr%s:%s  %s/%s  OR=%s  p=%.2e" % (t[1], t[2], t[3], t[4], t[5], t[6], t[0]))
