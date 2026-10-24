#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
在 GWAS Catalog 里找某个性状的所有研究, 并判断哪些有可下载的 harmonised 汇总统计
------------------------------------------------
用法: python _find_gwas.py <assoc.zip> <harmonised_list.txt> <性状关键字>
输出: 按 (病例数+对照数) 排序的研究清单, 并标注是否有可下载汇总统计。
"""
import csv, sys, zipfile, io, os, re


def open_maybe_zip(path):
    if path.endswith('.zip'):
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist() if n.endswith(('.tsv', '.txt'))]
        return io.TextIOWrapper(z.open(names[0]), encoding='utf-8', errors='replace')
    return open(path, encoding='utf-8', errors='replace')


def main():
    assoc, hlist, kw = sys.argv[1], sys.argv[2], sys.argv[3].lower()

    have = set()
    with open(hlist, encoding='utf-8', errors='replace') as f:
        for line in f:
            m = re.search(r'(GCST\d+)', line)
            if m:
                have.add(m.group(1))
    print("可下载汇总统计的 GCST 共 %d 个\n" % len(have))

    # STUDY ACCESSION -> {trait, n, author, pubmed}
    stud = {}
    with open_maybe_zip(assoc) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        cols = rdr.fieldnames
        ctrait = next(c for c in cols if c.strip().lower() == 'disease/trait')
        cmap = next((c for c in cols if 'mapped_trait' in c.lower()), None)
        cacc = next((c for c in cols if 'study accession' in c.lower()), None)
        csamp = next((c for c in cols if 'initial sample' in c.lower()), None)
        cauth = next((c for c in cols if 'first author' in c.lower()), None)
        cpm = next((c for c in cols if 'pubmed' in c.lower()), None)
        cdate = next((c for c in cols if 'date added' in c.lower()), None)
        for r in rdr:
            tr = (r.get(ctrait) or '')
            mp = (r.get(cmap) or '') if cmap else ''
            if kw not in (tr + ' ' + mp).lower():
                continue
            acc = (r.get(cacc) or '').strip()
            if not acc or acc in stud:
                continue
            stud[acc] = {
                'trait': tr.strip(), 'mapped': mp.strip(),
                'sample': (r.get(csamp) or '').strip(),
                'author': (r.get(cauth) or '').strip(),
                'pmid': (r.get(cpm) or '').strip(),
                'date': (r.get(cdate) or '').strip(),
            }

    def ntotal(s):
        """从 sample 描述里解析总数, 形如 '14,802 cases, 26,703 controls'"""
        nums = [int(x.replace(',', '')) for x in re.findall(r'(\d[\d,]*)', s)]
        return max(nums) if nums else 0

    rows = sorted(stud.items(), key=lambda kv: -ntotal(kv[1]['sample']))
    print("性状包含 '%s' 的研究: %d 个\n" % (kw, len(rows)))
    print("%-14s %-9s %-11s %-46s %-30s %s" %
          ("GCST", "可下载", "PMID", "性状", "样本量", "第一作者"))
    print("-" * 140)
    for acc, v in rows[:30]:
        print("%-14s %-9s %-11s %-46s %-30s %s" %
              (acc, "是" if acc in have else "-", v['pmid'][:11],
               v['trait'][:46], v['sample'][:30], v['author']))
    print()
    # 只看可下载的
    dl = [(a, v) for a, v in rows if a in have]
    print("其中**有可下载汇总统计**的 %d 个 (按样本量降序):" % len(dl))
    for acc, v in dl[:12]:
        print("   %-14s %-46s %s  (PMID %s, %s)" %
              (acc, v['trait'][:46], v['sample'][:40], v['pmid'], v['author']))


if __name__ == '__main__':
    main()
