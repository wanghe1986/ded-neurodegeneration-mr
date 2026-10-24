#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
单 SNP Wald ratio + Steiger 逐位点诊断
------------------------------------------------
目的: 判断 MR 的合并效应是否由个别 SNP 驱动, 以及每个位点的方向是否一致。
用法: python _wald.py <暴露iv.tsv> <结局out.tsv> [标签]
      [--ncx N] [--nkx N] [--ncy N] [--nky N]

!! 重要 !!  GWAS Catalog harmonised 文件里的 n_cases / n_controls 常常是
           「逐 SNP 的分批样本量」(如 PD 文件首行 1169/968), 不是研究总数。
           直接用会让 Steiger 的 r2 被放大数十倍、方向判断完全失真。
           必须用论文报告的总样本量覆盖: --ncy/--nky。
"""
import csv, math, sys, argparse


def f(x):
    x = (x or '').strip()
    if x in ('', 'NA', 'NaN', '.', 'None'): return None
    try: return float(x)
    except ValueError: return None


def rd(p):
    with open(p, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


ap = argparse.ArgumentParser(add_help=False)
ap.add_argument('--ncx', type=float)
ap.add_argument('--nkx', type=float)
ap.add_argument('--ncy', type=float)
ap.add_argument('--nky', type=float)
known, rest = ap.parse_known_args()
pos = [a for a in sys.argv[1:] if not a.startswith('-') and not a.startswith('--')]

iv = rd(pos[0]); out = rd(pos[1])
lab = pos[2] if len(pos) > 2 else ''
o = {r['SNP']: r for r in out}

# 有效样本量 (二分类: 4*ncase*nctrl/(ncase+nctrl))
def neff(rs, c0, k0):
    if c0 and k0: return 4 * c0 * k0 / (c0 + k0)
    c = f(rs[0].get('n_cases')); k = f(rs[0].get('n_controls'))
    if c and k: return 4 * c * k / (c + k)
    return None

nx = neff(iv, known.ncx, known.nkx)
ny = neff(out, known.ncy, known.nky)

rows = []
drop = []
for r in iv:
    s = r['SNP']
    if s not in o:
        drop.append((s, '结局未覆盖')); continue
    q = o[s]
    be, se_e, pe = f(r['beta']), f(r['se']), f(r['pval'])
    bo, se_o, po = f(q['beta']), f(q['se']), f(q['pval'])
    if None in (be, se_e, bo, se_o): drop.append((s, '统计量缺失')); continue
    ea, oa = r['effect_allele'].upper(), r['other_allele'].upper()
    qea, qoa = q['effect_allele'].upper(), q['other_allele'].upper()
    if ea == qea and oa == qoa:
        pass
    elif ea == qoa and oa == qea:
        bo = -bo                      # 等位基因互换, 结局效应取反
    else:
        drop.append((s, '等位基因不兼容')); continue
    if be == 0: drop.append((s, '暴露效应为0')); continue
    w = bo / be
    wse = abs(se_o / be)              # 一阶近似
    r2x = be ** 2 / (be ** 2 + se_e ** 2 * nx) if nx else float('nan')
    r2y = bo ** 2 / (bo ** 2 + se_o ** 2 * ny) if ny else float('nan')
    rows.append([s, r['chr'], int(float(r['pos'])), ea, oa, be, pe, bo, po,
                 w, wse, r2x, r2y, r2x > r2y])

rows.sort(key=lambda x: x[9])
print("\n===== 单 SNP 诊断: %s =====" % lab)
print("暴露有效 N=%s, 结局有效 N=%s" % (int(nx) if nx else '-', int(ny) if ny else '-'))
print("纳入 %d 个 SNP, 剔除 %d 个" % (len(rows), len(drop)))
for s, why in drop: print("   剔除 %-16s (%s)" % (s, why))
print()
print("%-16s %-14s %8s %9s %9s %10s  %8s %8s %s" %
      ("SNP", "位置", "beta_exp", "beta_out", "Wald_OR", "Wald_p", "r2_exp", "r2_out", "Steiger"))
print("-" * 108)
for x in rows:
    z = x[9] / x[10] if x[10] else 0
    pw = math.erfc(abs(z) / math.sqrt(2))
    print("%-16s chr%s:%-10d %+8.4f %+9.4f %9.3f %10.3f  %8.2e %8.2e %s" %
          (x[0], x[1], x[2], x[5], x[7], math.exp(x[9]), pw, x[11], x[12],
           "OK" if x[13] else "反向"))
print("-" * 108)
neg = sum(1 for x in rows if x[9] < 0)
print("方向: Wald<0 的 SNP 共 %d / %d" % (neg, len(rows)))
st = sum(1 for x in rows if x[13])
print("Steiger 方向正确: %d / %d (%.2f)" % (st, len(rows), st / len(rows) if rows else 0))
