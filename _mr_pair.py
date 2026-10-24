# -*- coding: utf-8 -*-
"""通用配对 MR（Python 复刻 mr_core.R::mr_estimate），供 R 不可用时重算。

用法:
  python _mr_pair.py <iv.tsv> <out.tsv> <标签> [ncase_x nctrl_x ncase_y nctrl_y] [--csv 文件]

与 mr_core.R 的差异仅在于等位基因协调同时支持「链翻转」：
  direct           暴露/结局同链同序                     -> 不翻转
  swapped          暴露/结局同链反序                     -> 两侧 beta 同时取反
  flipped          结局记在互补链上、同序                 -> 不翻转
  flipped_swapped  结局记在互补链上、反序                 -> 两侧 beta 同时取反
（暴露侧回文 SNP 构建工具集时已剔除，故翻转类无方向歧义）

复刻细节（与 mr_core.R 一致）：
  IVW   : w = 1/sy^2；b = S(w*bx*by)/S(w*bx^2)；se = sqrt(1/S(w*bx^2))
  Q     : S(w*(by - b*bx)^2)，df = n-1，p = chi2 精确上尾
  RE    : se *= sqrt(max(1, Q/df))
  Egger : 加权 lm(by ~ bx)，系数 p 用 t 分布（df = n-2）
  WM    : 权重 w*bx^2，bootstrap 重抽样 SNP（B=1000，seed 42），se = sd
"""
import csv
import math
import random
import sys

from scipy.stats import chi2, norm, t as tdist

COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def rd(path):
    with open(path, encoding='utf-8') as f:
        return {r['SNP']: r for r in csv.DictReader(f, delimiter='\t')}


def harmonize(iv, out):
    """返回 (bx, sx, by, sy, names, 状态计数)"""
    bx, sx, by, sy, names = [], [], [], [], []
    cnt = {'direct': 0, 'swapped': 0, 'flipped': 0, 'flipped_swapped': 0,
           'ambiguous': 0, 'mismatch': 0, 'missing': 0}
    for snp, a in iv.items():
        if snp not in out:
            cnt['missing'] += 1
            continue
        o = out[snp]
        e1, o1 = a['effect_allele'].upper(), a['other_allele'].upper()
        e2, o2 = o['effect_allele'].upper(), o['other_allele'].upper()
        b_x, b_y = float(a['beta']), float(o['beta'])
        c1, c2 = COMP.get(e1), COMP.get(o1)
        if (e1, o1) == (e2, o2):
            st = 'direct'
        elif (e1, o1) == (o2, e2):
            st = 'swapped'
            b_x, b_y = -b_x, -b_y
        elif (c1, c2) == (e2, o2):
            st = 'flipped'
        elif (c1, c2) == (o2, e2):
            st = 'flipped_swapped'
            b_x, b_y = -b_x, -b_y
        elif c1 == o1:          # 回文，方向不可判定
            st = 'ambiguous'
        else:
            st = 'mismatch'
        cnt[st] += 1
        if st not in ('direct', 'swapped', 'flipped', 'flipped_swapped'):
            continue
        bx.append(b_x); sx.append(float(a['se']))
        by.append(b_y); sy.append(float(o['se']))
        names.append(snp)
    return bx, sx, by, sy, names, cnt


def estimate(bx, sx, by, sy):
    n = len(bx)
    w = [1 / s ** 2 for s in sy]
    den = sum(wi * a * a for wi, a in zip(w, bx))
    b_ivw = sum(wi * a * b for wi, a, b in zip(w, bx, by)) / den
    se_ivw = den ** -0.5
    Q = sum(wi * (b - b_ivw * a) ** 2 for wi, a, b in zip(w, bx, by))
    df = n - 1
    Q_p = chi2.sf(Q, df)
    se_re = se_ivw * max(1.0, Q / df) ** 0.5

    S = sum(w)
    Sx = sum(wi * a for wi, a in zip(w, bx))
    Sy = sum(wi * b for wi, b in zip(w, by))
    Sxx = den
    Sxy = sum(wi * a * b for wi, a, b in zip(w, bx, by))
    D = S * Sxx - Sx ** 2
    b_eg = (S * Sxy - Sx * Sy) / D
    b0 = (Sy - b_eg * Sx) / S
    Syy = sum(wi * b * b for wi, b in zip(w, by))
    s2 = (Syy - b0 * Sy - b_eg * Sxy) / (n - 2)
    se_eg = (s2 * S / D) ** 0.5
    se_b0 = (s2 * Sxx / D) ** 0.5
    p_eg = 2 * tdist.sf(abs(b_eg / se_eg), n - 2)
    p_b0 = 2 * tdist.sf(abs(b0 / se_b0), n - 2)

    ratio = [b / a for a, b in zip(bx, by)]
    wm_w = [wi * a * a for wi, a in zip(w, bx)]

    def wm(rr, ww):
        order = sorted(range(len(rr)), key=lambda i: rr[i])
        tot = sum(ww)
        acc = 0.0
        for i in order:
            acc += ww[i]
            if acc >= tot / 2:
                return rr[i]
        return rr[order[-1]]

    b_wm = wm(ratio, wm_w)
    random.seed(42)
    bs = []
    for _ in range(1000):
        s = [random.randrange(n) for _ in range(n)]
        bs.append(wm([ratio[i] for i in s], [wm_w[i] for i in s]))
    mean = sum(bs) / len(bs)
    se_wm = (sum((x - mean) ** 2 for x in bs) / (len(bs) - 1)) ** 0.5

    rows = [('ivw_fe', b_ivw, se_ivw, 2 * norm.sf(abs(b_ivw / se_ivw))),
            ('ivw_re', b_ivw, se_re, 2 * norm.sf(abs(b_ivw / se_re))),
            ('mr_egger', b_eg, se_eg, p_eg),
            ('weighted_median', b_wm, se_wm, 2 * norm.sf(abs(b_wm / se_wm)))]
    out = []
    for name, b, se, p in rows:
        out.append(dict(method=name, n=n, b=b, se=se, pval=p,
                        OR=math.exp(b), LCI=math.exp(b - 1.959964 * se),
                        UCI=math.exp(b + 1.959964 * se)))
    return out, Q, df, Q_p, b0, se_b0, p_b0


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    ivf, outf, label = args[0], args[1], args[2]
    nums = [float(x) for x in args[3:7]]
    iv, out = rd(ivf), rd(outf)
    bx, sx, by, sy, names, cnt = harmonize(iv, out)
    print('标签: %s' % label)
    print('IV %d 个 | 结局命中 %d | 协调后可用 %d' % (len(iv), len(iv) - cnt['missing'], len(bx)))
    print('协调状态:', {k: v for k, v in cnt.items() if v})
    rows, Q, df, Q_p, b0, se_b0, p_b0 = estimate(bx, sx, by, sy)
    print('\n%-16s %4s %8s %8s %8s %10s' % ('method', 'n', 'OR', 'LCI', 'UCI', 'P'))
    for r in rows:
        print('%-16s %4d %8.4f %8.4f %8.4f %10.4g' %
              (r['method'], r['n'], r['OR'], r['LCI'], r['UCI'], r['pval']))
    print('\nQ = %.4f (df %d), Q_P = %.4f' % (Q, df, Q_p))
    print('Egger intercept = %.6g (SE %.6g), P = %.4f' % (b0, se_b0, p_b0))

    if len(nums) == 4 and 'n_cases' in next(iter(iv.values())):
        pass
    # Steiger：需要研究层面样本量
    if len(nums) == 4:
        ncase_x, nctrl_x, ncase_y, nctrl_y = nums
        nx = 4 * ncase_x * nctrl_x / (ncase_x + nctrl_x)
        ny = 4 * ncase_y * nctrl_y / (ncase_y + nctrl_y)
        r2x = [a ** 2 / (a ** 2 + s ** 2 * nx) for a, s in zip(bx, sx)]
        r2y = [b ** 2 / (b ** 2 + s ** 2 * ny) for b, s in zip(by, sy)]
        ok = sum(1 for a, b in zip(r2x, r2y) if a > b)
        print('Steiger: 方向正确 %d/%d, r2暴露均值=%.3g, r2结局均值=%.3g'
              % (ok, len(bx), sum(r2x) / len(r2x), sum(r2y) / len(r2y)))

    if '--csv' in sys.argv:
        fn = sys.argv[sys.argv.index('--csv') + 1]
        new = not __import__('os').path.exists(fn)
        with open(fn, 'a', newline='', encoding='utf-8') as f:
            wr = csv.writer(f)
            if new:
                wr.writerow(['pair', 'nsnp_iv', 'nsnp_used', 'method', 'b', 'se', 'pval',
                             'nsnp', 'Q', 'Q_pval', 'egger_intercept', 'egger_pval',
                             'OR', 'LCI', 'UCI'])
            for r in rows:
                wr.writerow([label, len(iv), len(bx), r['method'], r['b'], r['se'], r['pval'],
                             r['n'], Q, Q_p, b0, p_b0, r['OR'], r['LCI'], r['UCI']])
        print('\n已追加 %s' % fn)


if __name__ == '__main__':
    main()
