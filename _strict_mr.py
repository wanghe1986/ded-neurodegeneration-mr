# -*- coding: utf-8 -*-
"""补算严格工具集（P<5e-8, 9 loci）DED -> PD 的 MR 估计。

严格复刻 mr_core.R::mr_estimate()：
  IVW   : w = 1/sy^2；b = S(w*bx*by)/S(w*bx^2)；se = sqrt(1/S(w*bx^2))
  Q     : S(w*(by - b*bx)^2)，df = n-1，p = pchisq(upper)
  RE    : se *= sqrt(max(1, Q/df))
  Egger : lm(by ~ bx, weights = w)，系数 p 用 t 分布（df = n-2）
  WM    : 权重 w*bx^2（= 1/se_i^2），bootstrap 重抽样 SNP（B=1000），se = sd
"""
import csv
import random

from scipy.stats import chi2, norm, t as tdist

IV = 'data/processed/iv_DED.tsv'
OUT = 'data/processed/out_PD_for_DED_p5e8.tsv'


def rd(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    iv = {r['SNP']: r for r in rd(IV)}
    out = {r['SNP']: r for r in rd(OUT)}

    bx, sx, by, sy, names = [], [], [], [], []
    for snp in iv:
        if snp not in out:
            continue
        a, b = iv[snp], out[snp]
        b_x, b_y = float(a['beta']), float(b['beta'])
        if (a['effect_allele'], a['other_allele']) == (b['other_allele'], b['effect_allele']):
            b_x, b_y = -b_x, -b_y
        elif (a['effect_allele'], a['other_allele']) != (b['effect_allele'], b['other_allele']):
            raise SystemExit('等位基因不兼容: %s' % snp)
        bx.append(b_x); sx.append(float(a['se']))
        by.append(b_y); sy.append(float(b['se']))
        names.append(snp)

    n = len(bx)
    print('匹配 SNP: %d / IV %d  ->  %s' % (n, len(iv), ', '.join(names)))

    w = [1 / s ** 2 for s in sy]
    num = sum(wi * a * b for wi, a, b in zip(w, bx, by))
    den = sum(wi * a * a for wi, a in zip(w, bx))
    b_ivw = num / den
    se_ivw = den ** -0.5
    Q = sum(wi * (b - b_ivw * a) ** 2 for wi, a, b in zip(w, bx, by))
    df = n - 1
    Q_p = chi2.sf(Q, df)
    se_re = se_ivw * max(1.0, Q / df) ** 0.5

    # MR-Egger：加权线性回归，t 分布 df = n-2
    S = sum(w)
    Sx = sum(wi * a for wi, a in zip(w, bx))
    Sy = sum(wi * b for wi, b in zip(w, by))
    Sxx = sum(wi * a * a for wi, a in zip(w, bx))
    Sxy = sum(wi * a * b for wi, a, b in zip(w, bx, by))
    D = S * Sxx - Sx ** 2
    b_eg = (S * Sxy - Sx * Sy) / D
    b0 = (Sy - b_eg * Sx) / S
    # lm() 的残差方差：SSE/(n-2)，SSE = Syy - b0*Sy - b_eg*Sxy（加权）
    Syy = sum(wi * b * b for wi, b in zip(w, by))
    SSE = Syy - b0 * Sy - b_eg * Sxy
    s2 = SSE / (n - 2)
    se_eg = (s2 * S / D) ** 0.5
    se_b0 = (s2 * Sxx / D) ** 0.5
    p_eg = 2 * tdist.sf(abs(b_eg / se_eg), n - 2)
    p_b0 = 2 * tdist.sf(abs(b0 / se_b0), n - 2)

    # Weighted median：权重 w*bx^2，bootstrap 重抽样
    ratio = [b / a for a, b in zip(bx, by)]
    wm_w = [wi * a * a for wi, a in zip(w, bx)]

    def wm(rr, ww):
        o = sorted(range(len(rr)), key=lambda i: rr[i])
        tot = sum(ww)
        acc = 0.0
        for i in o:
            acc += ww[i]
            if acc >= tot / 2:
                return rr[i]
        return rr[o[-1]]

    b_wm = wm(ratio, wm_w)
    random.seed(42)
    bs = []
    for _ in range(1000):
        s = [random.randrange(n) for _ in range(n)]
        bs.append(wm([ratio[i] for i in s], [wm_w[i] for i in s]))
    mean = sum(bs) / len(bs)
    se_wm = (sum((x - mean) ** 2 for x in bs) / (len(bs) - 1)) ** 0.5

    print('\n%-16s %4s %8s %8s %8s %10s' % ('method', 'n', 'OR', 'LCI', 'UCI', 'P'))
    res = []
    for name, b, se, p in [('ivw_fe', b_ivw, se_ivw, 2 * norm.sf(abs(b_ivw / se_ivw))),
                           ('ivw_re', b_ivw, se_re, 2 * norm.sf(abs(b_ivw / se_re))),
                           ('mr_egger', b_eg, se_eg, p_eg),
                           ('weighted_median', b_wm, se_wm, 2 * norm.sf(abs(b_wm / se_wm)))]:
        OR = pow(2.718281828459045, b)
        lo, hi = pow(2.718281828459045, b - 1.959964 * se), pow(2.718281828459045, b + 1.959964 * se)
        print('%-16s %4d %8.4f %8.4f %8.4f %10.3g' % (name, n, OR, lo, hi, p))
        res.append([name, b, se, p, OR, lo, hi])

    print('\nQ = %.4f (df %d), Q_P = %.4f' % (Q, df, Q_p))
    print('Egger intercept = %.6g (SE %.6g), P = %.4f' % (b0, se_b0, p_b0))

    with open('results/strict_p5e8_PD.csv', 'w', newline='', encoding='utf-8') as f:
        wr = csv.writer(f)
        wr.writerow(['pair', 'nsnp', 'method', 'b', 'se', 'pval', 'Q', 'Q_pval',
                     'egger_intercept', 'egger_pval', 'OR', 'LCI', 'UCI'])
        for name, b, se, p, OR, lo, hi in res:
            wr.writerow(['DED_EUR9 -> PD', n, name, b, se, p, Q, Q_p, b0, p_b0, OR, lo, hi])
    print('\n已写 results/strict_p5e8_PD.csv')


if __name__ == '__main__':
    main()
