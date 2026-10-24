#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1/B 方案 · 通用 GWAS 汇总统计预处理
------------------------------------------------
两个模式:
  iv      从暴露 GWAS 构建工具变量 (P 阈值 + 距离剪枝去 LD + F 统计)
  extract 给定 SNP 列表, 从结局 GWAS 提取对应统计量

用法:
  python _mr_prep.py iv      --gwas data/raw/DED_MVP_GCST90475895.h.tsv.gz \
                             --out  data/processed/iv_DED.tsv --p 5e-8 --window 10000000
  python _mr_prep.py extract --gwas data/raw/MS_IMSGC2019_GCST005531.h.tsv.gz \
                             --snps data/processed/iv_DED.tsv \
                             --out  data/processed/out_MS_for_DED.tsv

说明: 本机无 PLINK + 1000G 参考面板, 故用「按位点距离剪枝」代替 LD clumping。
      窗口默认 10 Mb (严于 PLINK 常用的 r2<0.001/1 Mb), 保守但绝无残留连锁。
"""
import gzip, sys, csv, math, argparse, os

# 说明: 候选按顺序匹配, 新增的格式名一律放末尾, 以免影响已解析过的文件。
#       FinnGen R12: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...
COLMAP = {
    'snp':   ['rsid', 'hm_rsid', 'variant_id', 'SNP', 'rs_id', 'rsids'],
    'chr':   ['chromosome', 'hm_chrom', 'chr', 'chrom'],
    'pos':   ['base_pair_location', 'hm_pos', 'pos', 'bp'],
    'ea':    ['effect_allele', 'hm_effect_allele', 'A1', 'alt'],
    'oa':    ['other_allele', 'hm_other_allele', 'A2', 'ref'],
    'beta':  ['beta', 'hm_beta'],
    'se':    ['standard_error', 'se', 'SE', 'sebeta'],
    'or':    ['odds_ratio', 'hm_odds_ratio', 'OR'],
    'p':     ['p_value', 'pval', 'p_value_gc', 'P'],
    'eaf':   ['effect_allele_frequency', 'hm_effect_allele_frequency', 'freq', 'maf'],
    'cilo':  ['ci_lower', 'hm_ci_lower'],
    'cihi':  ['ci_upper', 'hm_ci_upper'],
    'ncase': ['num_cases', 'n_cases', 'Ncase', 'ncase'],
    'nctrl': ['num_controls', 'n_controls', 'Ncontrol', 'ncontrol'],
    'n':     ['n', 'N', 'n_total', 'Neff'],
}

def open_maybe_gz(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'rt', encoding='utf-8', errors='replace')

def resolve_cols(header):
    """把实际表头映射到标准列名; 返回 {std: idx}"""
    h = [c.strip().lstrip('#') for c in header]
    idx = {}
    for c in h:
        idx.setdefault(c.lower(), None)
    low = [c.lower() for c in h]
    out = {}
    for std, cands in COLMAP.items():
        for c in cands:
            if c.lower() in low:
                out[std] = low.index(c.lower())
                break
    return out, h

def fnum(x):
    if x is None: return None
    x = x.strip()
    if x in ('', 'NA', 'NaN', 'None', '.', '-'): return None
    try: return float(x)
    except ValueError: return None

def get_beta_se(row, ci):
    """按优先级取 (beta, se): beta/se -> log(OR) + CI 反推 SE -> log(OR) + P 反推"""
    beta = fnum(row[ci['beta']]) if 'beta' in ci else None
    se   = fnum(row[ci['se']])   if 'se'   in ci else None
    orv  = fnum(row[ci['or']])   if 'or'   in ci else None
    lo   = fnum(row[ci['cilo']]) if 'cilo' in ci else None
    hi   = fnum(row[ci['cihi']]) if 'cihi' in ci else None
    p    = fnum(row[ci['p']])    if 'p'    in ci else None

    if beta is None and orv is not None and orv > 0:
        beta = math.log(orv)
    if beta is None:
        return None, None
    if se is not None:
        return beta, se
    # CI 反推 (95% CI)
    if lo is not None and hi is not None and lo > 0 and hi > 0:
        se = (math.log(hi) - math.log(lo)) / (2 * 1.959964)
        return beta, se
    # P 反推
    if p is not None and p > 0 and abs(beta) > 0:
        z = abs(_qnorm(p / 2.0))
        if z > 0:
            return beta, abs(beta) / z
    return beta, None

def _qnorm(q):
    """标准正态分位数 (Acklam 近似), q in (0,1)"""
    if q <= 0: return -40.0
    if q >= 1: return 40.0
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if q < pl:
        t = math.sqrt(-2 * math.log(q))
        return (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5]) / ((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1)
    if q > ph:
        t = math.sqrt(-2 * math.log(1 - q))
        return -(((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5]) / ((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1)
    t = q - 0.5; r = t * t
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*t / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def iter_rows(path):
    with open_maybe_gz(path) as fh:
        rdr = csv.reader(fh, delimiter='\t')
        try:
            header = next(rdr)
        except StopIteration:
            return
        ci, raw_h = resolve_cols(header)
        for need in ('snp', 'p'):
            if need not in ci:
                sys.exit("错误: %s 缺少可识别的 %s 列。表头: %s" % (path, need, raw_h))
        yield ci, rdr


def mode_iv(args):
    ci = None
    rows = []
    ncase = nctrl = None
    for ci, rdr in iter_rows(args.gwas):
        for row in rdr:
            if len(row) <= max(ci.values()): continue
            p = fnum(row[ci['p']])
            if p is None or p > args.p: continue
            beta, se = get_beta_se(row, ci)
            if beta is None or se is None or se <= 0: continue
            ea = (row[ci['ea']] if 'ea' in ci else '').strip().upper()
            oa = (row[ci['oa']] if 'oa' in ci else '').strip().upper()
            # 剔除回文 SNP (A/T, G/C) —— 方向不可判定
            if {ea, oa} in ({'A', 'T'}, {'G', 'C'}): continue
            try:
                ch = row[ci['chr']].strip(); pos = int(float(row[ci['pos']]))
            except Exception:
                continue
            if 'ncase' in ci and ncase is None:
                ncase = fnum(row[ci['ncase']])
            if 'nctrl' in ci and nctrl is None:
                nctrl = fnum(row[ci['nctrl']])
            eaf = fnum(row[ci['eaf']]) if 'eaf' in ci else None
            rows.append([row[ci['snp']].strip(), ch, pos, ea, oa, beta, se, p, eaf])
    if not rows:
        sys.exit("没有 SNP 通过 P < %g 阈值" % args.p)

    rows.sort(key=lambda r: (str(r[1]), r[2]))
    kept = []
    for r in sorted(rows, key=lambda x: x[7]):      # 按 P 升序
        ok = True
        for k in kept:
            if k[1] == r[1] and abs(k[2] - r[2]) < args.window:
                ok = False; break
        if ok: kept.append(r)

    kept.sort(key=lambda r: (str(r[1]), r[2]))
    with open(args.out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        cols = ['SNP', 'chr', 'pos', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'eaf']
        if ncase is not None and nctrl is not None:
            for r in kept: r.append(int(ncase)); r.append(int(nctrl))
            cols += ['n_cases', 'n_controls']
        for r in kept:
            r.append(round(r[5] * r[5] / (r[6] * r[6]), 2));  # F
        cols.append('F')
        w.writerow(cols); w.writerows(kept)
    print("IV 构建: 全基因组显著 %d 个 -> 距离剪枝后 %d 个独立位点" % (len(rows), len(kept)))
    print("样本量: cases=%s controls=%s" % (ncase, nctrl))
    print("输出: %s" % args.out)
    lo = min(r[6] for r in kept) if False else None
    print("F 范围: %.1f - %.1f" % (min(r[-1] for r in kept), max(r[-1] for r in kept)))
    for r in kept[:30]:
        print("   %-16s %s:%-10d %s/%s  beta=%+.4f  p=%.2e  F=%.0f" %
              (r[0], r[1], r[2], r[3], r[4], r[5], r[7], r[-1]))


def mode_extract(args):
    want = {}
    with open(args.snps) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            want[row['SNP'].strip()] = True
    print("待提取 SNP: %d" % len(want))
    got = {}
    ncase = nctrl = None
    for ci, rdr in iter_rows(args.gwas):
        for row in rdr:
            if len(row) <= max(ci.values()): continue
            s = row[ci['snp']].strip()
            if s not in want or s in got: continue
            beta, se = get_beta_se(row, ci)
            if beta is None or se is None or se <= 0: continue
            if 'ncase' in ci and ncase is None: ncase = fnum(row[ci['ncase']])
            if 'nctrl' in ci and nctrl is None: nctrl = fnum(row[ci['nctrl']])
            ea = (row[ci['ea']] if 'ea' in ci else '').strip().upper()
            oa = (row[ci['oa']] if 'oa' in ci else '').strip().upper()
            eaf = fnum(row[ci['eaf']]) if 'eaf' in ci else None
            got[s] = [s, ea, oa, beta, se, fnum(row[ci['p']]), eaf]
    with open(args.out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        cols = ['SNP', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'eaf']
        if ncase is not None and nctrl is not None:
            for v in got.values(): v += [int(ncase), int(nctrl)]
            cols += ['n_cases', 'n_controls']
        w.writerow(cols)
        for v in got.values(): w.writerow(v)
    print("命中: %d / %d  -> %s" % (len(got), len(want), args.out))
    if len(got) < len(want):
        miss = [s for s in want if s not in got]
        print("未命中(结局 GWAS 未覆盖): %d 个, 例: %s" % (len(miss), ', '.join(miss[:8])))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['iv', 'extract'])
    ap.add_argument('--gwas', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--p', type=float, default=5e-8)
    ap.add_argument('--window', type=int, default=10000000)
    ap.add_argument('--snps', help='extract 模式: 含 SNP 列的 TSV')
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    (mode_iv if a.mode == 'iv' else mode_extract)(a)
