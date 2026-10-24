#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
为 Figure 3 (单 SNP Wald ratio) 与 Figure 4 (留一法) 生成精确绘图数据。
IVW 公式与 mr_core.R 完全一致:  w = 1/sy^2; b = sum(w*bx*by)/sum(w*bx^2); se = sqrt(1/sum(w*bx^2))
输出 JSON 到 results/figdata_*.json, 并打印与 _loo.R / _wald.py 的核对结果。
"""
import csv, json, math, os, sys

IV = "data/processed/iv_DED_eur_p5e7.tsv"
OUT = "data/processed/out_PD_for_iv_DED_eur_p5e7.tsv"

def f(x):
    x = (x or '').strip()
    if x in ('', 'NA', 'NaN', '.', 'None'):
        return None
    try:
        return float(x)
    except ValueError:
        return None

def rd(p):
    with open(p, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

iv = rd(IV)
out = rd(OUT)
o = {r['SNP']: r for r in out}

rows = []
dropped = []
for r in iv:
    s = r['SNP']
    if s not in o:
        dropped.append((s, 'not in outcome'))
        continue
    q = o[s]
    be, se_e = f(r['beta']), f(r['se'])
    bo, se_o = f(q['beta']), f(q['se'])
    if None in (be, se_e, bo, se_o):
        dropped.append((s, 'missing stat'))
        continue
    ea, oa = r['effect_allele'].upper(), r['other_allele'].upper()
    qea, qoa = q['effect_allele'].upper(), q['other_allele'].upper()
    if ea == qea and oa == qoa:
        pass
    elif ea == qoa and oa == qea:
        bo = -bo
    else:
        dropped.append((s, 'allele mismatch'))
        continue
    if be == 0:
        dropped.append((s, 'zero exposure effect'))
        continue
    rows.append({'snp': s, 'bx': be, 'sx': se_e, 'by': bo, 'sy': se_o})

print("纳入 %d 个 SNP, 剔除 %d 个" % (len(rows), len(dropped)))
for s, why in dropped:
    print("   剔除 %-16s (%s)" % (s, why))

# ---------- IVW (固定效应), 与 mr_core.R 同式 ----------
def ivw(rs):
    w = [1.0 / r['sy'] ** 2 for r in rs]
    num = sum(wi * r['bx'] * r['by'] for wi, r in zip(w, rs))
    den = sum(wi * r['bx'] ** 2 for wi, r in zip(w, rs))
    b = num / den
    se = math.sqrt(1.0 / den)
    z = b / se
    p = math.erfc(abs(z) / math.sqrt(2))
    return b, se, p

def ci(b, se):
    return math.exp(b - 1.959963985 * se), math.exp(b + 1.959963985 * se)

b0, se0, p0 = ivw(rows)
print("\n全样本 IVW: b=%+.5f  OR=%.4f (%.4f-%.4f)  p=%.5g" % ((b0,) + (math.exp(b0),) + ci(b0, se0) + (p0,)))

# ---------- Figure 3: 单 SNP Wald ratio ----------
wald = []
for r in rows:
    bw = r['by'] / r['bx']
    sw = abs(r['sy'] / r['bx'])
    lo, hi = ci(bw, sw)
    wald.append({'snp': r['snp'], 'OR': math.exp(bw), 'lci': lo, 'uci': hi,
                 'p': math.erfc(abs(bw / sw) / math.sqrt(2))})
wald.sort(key=lambda d: d['OR'])
n_neg = sum(1 for d in wald if d['OR'] < 1)
print("\n单 SNP Wald: OR<1 共 %d / %d" % (n_neg, len(wald)))
print("  OR 范围 %.3f - %.3f" % (wald[0]['OR'], wald[-1]['OR']))

# ---------- Figure 4: 留一法 ----------
loo = []
for i, r in enumerate(rows):
    sub = rows[:i] + rows[i + 1:]
    b, se, p = ivw(sub)
    lo, hi = ci(b, se)
    loo.append({'snp': r['snp'], 'OR': math.exp(b), 'lci': lo, 'uci': hi, 'p': p})
loo_sorted = sorted(loo, key=lambda d: d['snp'])
print("\n留一法 (n=%d):" % len(loo))
print("  OR 范围 %.3f - %.3f" % (min(d['OR'] for d in loo), max(d['OR'] for d in loo)))
print("  p  范围 %.4g - %.4g" % (min(d['p'] for d in loo), max(d['p'] for d in loo)))
print("  方向全部一致: %s" % (all(d['OR'] < 1 for d in loo) if b0 < 0 else all(d['OR'] > 1 for d in loo)))
print("\n对照 _loo.R 输出: OR 0.596-0.655, p 0.0003885-0.003997")

os.makedirs("results", exist_ok=True)
json.dump({'wald': wald, 'loo': sorted(loo, key=lambda d: d['OR']),
           'overall': {'OR': math.exp(b0), 'lci': ci(b0, se0)[0], 'uci': ci(b0, se0)[1], 'p': p0}},
          open("results/figdata_wald_loo.json", "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
print("\n已写出 results/figdata_wald_loo.json")
