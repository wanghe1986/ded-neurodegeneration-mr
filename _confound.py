#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
工具变量混杂排查 (PhenoScanner 的离线替代)
------------------------------------------------
思路: GWAS Catalog 收录了所有已发表的 GWAS 显著位点及其注释性状。
      把 25 个干眼工具变量拿到里面比对, 看它们是否同时是其他性状的已知位点。
      重点关注会混淆「干眼 -> 帕金森病」的通路: 吸烟、饮酒、BMI、教育、2型糖尿病等。

用法:
  python _confound.py <assoc.tsv> <snp列表文件> [--focus 关键字文件]
输出: 每个 SNP 命中的性状清单, 并按关注通路高亮。
"""
import csv, sys, os, zipfile, io, collections, argparse

# 关注通路: 这些若命中, 说明工具变量可能通过这些通路影响结局 (水平多效性)
FOCUS = {
    'smok': '吸烟', 'tobacco': '吸烟', 'cigarette': '吸烟', 'nicotine': '吸烟',
    'alcohol': '饮酒', 'drink': '饮酒',
    'body mass': 'BMI/肥胖', 'bmi': 'BMI/肥胖', 'obesity': 'BMI/肥胖',
    'educational': '教育', 'education': '教育', 'cognitive': '教育',
    'type 2 diabetes': '2型糖尿病', 'diabetes': '2型糖尿病',
    'glycemic': '2型糖尿病', 'glucose': '2型糖尿病', 'insulin': '2型糖尿病',
    'depression': '抑郁', 'neuroticism': '抑郁', 'schizophrenia': '精神分裂',
    'bipolar': '双相', 'anxiety': '焦虑',
    'blood pressure': '血压', 'hypertension': '血压',
    'coronary': '冠心病', 'stroke': '卒中',
    'c-reactive': '炎症', 'c reactive': '炎症',
    'rheumatoid': '自身免疫', 'lupus': '自身免疫', "sjogren": '自身免疫/干燥',
    'autoimmune': '自身免疫', 'inflammatory bowel': '自身免疫',
    'asthma': '过敏/哮喘', 'allerg': '过敏/哮喘', 'eczema': '过敏/哮喘',
    'hay fever': '过敏/哮喘', 'atopic': '过敏/哮喘',
    'parkinson': '帕金森病', 'alzheimer': '阿尔茨海默病',
    'multiple sclerosis': '多发性硬化',
    'loneliness': '社会因素', 'income': '社会因素',
    'physical activity': '体力活动', 'sleep': '睡眠',
    'urate': '代谢', 'lipid': '血脂', 'cholesterol': '血脂',
    'height': '身高',
}

P_CUT = 5e-8   # 只认全基因组显著的既往报道


def load_snps(path):
    out = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if s and s.lower() != 'snp':
                out.append(s)
    return out


def open_assoc(path):
    """支持 .tsv / .zip(内含单个 tsv)"""
    if path.endswith('.zip'):
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist() if n.endswith(('.tsv', '.txt'))]
        if not names:
            sys.exit("zip 内无 tsv")
        print("zip 内文件:", names)
        return io.TextIOWrapper(z.open(names[0]), encoding='utf-8', errors='replace')
    return open(path, encoding='utf-8', errors='replace')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('assoc')
    ap.add_argument('snps')
    ap.add_argument('--all', action='store_true',
                    help='输出全部命中而非仅关注通路')
    a = ap.parse_args()

    want = set(load_snps(a.snps))
    print("待查 SNP: %d 个" % len(want))

    hits = collections.defaultdict(list)
    nrow = 0
    with open_assoc(a.assoc) as fh:
        rdr = csv.DictReader(fh, delimiter='\t')
        cols = rdr.fieldnames
        print("列名(%d): %s" % (len(cols), ', '.join(cols[:8])) + ' ...')
        # !! 注意 !!  SNP_ID_CURRENT 列存的是裸数字 (13107325), 不带 rs 前缀。
        #            必须用 SNPS 列 (rs13107325)。顺序不能反。
        csnp = next((c for c in cols if c.strip().lower() == 'snps'), None) or \
               next((c for c in cols if c.strip().lower() ==
                     'snp_id_current'), None)
        ctrait = next((c for c in cols if c.strip().lower() ==
                       'disease/trait'), None)
        cmapped = next((c for c in cols if 'mapped_trait' in c.lower()), None)
        cp = next((c for c in cols if c.strip().lower().startswith('p-value')
                   or c.strip().lower() == 'p_value'), None)
        cgene = next((c for c in cols if 'mapped_gene' in c.lower()), None)
        cstudy = next((c for c in cols if 'study accession' in c.lower()), None)
        cauth = next((c for c in cols if 'first author' in c.lower()), None)
        if not csnp or not ctrait:
            sys.exit("找不到 SNP/性状 列")
        print("使用列: SNP=%s, 性状=%s, P=%s, 基因=%s" % (csnp, ctrait, cp, cgene))

        for r in rdr:
            nrow += 1
            s = (r.get(csnp) or '').strip()
            if not s or s not in want:
                continue
            p = r.get(cp) or ''
            try:
                pf = float(p.replace('E', 'e'))
            except Exception:
                pf = None
            if pf is not None and pf > P_CUT:
                continue
            tr = (r.get(ctrait) or '').strip()
            mp = (r.get(cmapped) or '').strip() if cmapped else ''
            ge = (r.get(cgene) or '').strip() if cgene else ''
            au = (r.get(cauth) or '').strip() if cauth else ''
            hits[s].append((pf, tr, mp, ge, au))

    print("扫描 %d 行, 命中 %d / %d 个 SNP\n" % (nrow, len(hits), len(want)))

    flagged = collections.defaultdict(set)
    print("=" * 100)
    for i, s in enumerate(sorted(want, key=lambda x: (x not in hits, x))):
        if s not in hits:
            print("%-16s  —  GWAS Catalog 无其他性状的既往报道" % s)
            continue
        rows = sorted(hits[s], key=lambda x: (x[0] is None, x[0] or 9))
        foc = []
        for pf, tr, mp, ge, au in rows:
            low = (tr + ' ' + mp).lower()
            for k, v in FOCUS.items():
                if k in low:
                    foc.append((v, tr, pf, ge, au))
                    flagged[s].add(v)
        if not foc and not a.all:
            print("%-16s  —  仅命中非关注性状 (%d 条: %s)" %
                  (s, len(rows), '; '.join(sorted(set(r[1] for r in rows))[:3])))
        else:
            print("%-16s  ⚠ 关注通路命中 %d 条" % (s, len(foc)))
            seen = set()
            for v, tr, pf, ge, au in foc:
                key = (v, tr)
                if key in seen:
                    continue
                seen.add(key)
                print("      [%s] %-52s p=%-10.1e %s (%s)" % (v, tr[:52], pf, ge[:22], au))

    print("=" * 100)
    if flagged:
        print("\n⚠ 需要警惕的工具变量:")
        for s, vs in sorted(flagged.items(), key=lambda x: -len(x[1])):
            print("   %-16s %s" % (s, '、'.join(sorted(vs))))
    else:
        print("\n未发现工具变量命中关注的混杂通路。")


if __name__ == '__main__':
    main()
