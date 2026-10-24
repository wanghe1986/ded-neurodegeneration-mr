#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MR 最小可检出效应 (MDE) 计算
------------------------------------------------
双向 MR 全阴性时, 审稿人第一个问题就是"是不是功效不足"。
本脚本按观测到的 SE 反推: 在 alpha=0.05 (双侧)、power=80% 下能检出的最小真实效应。

    |b_MDE| = (z_{0.975} + z_{0.80}) * SE = 2.8016 * SE
    OR_MDE = exp(b_MDE)  /  exp(-b_MDE)

注意: 该值以"估计量无弱工具衰减"为前提, 是乐观下界 (即真实可检出效应只会更大),
     稿件中应如此声明。
"""
import csv, math, sys

Z = 1.959963985 + 0.841621234   # 2.801585


def mde(se):
    b = Z * se
    return b, math.exp(b), math.exp(-b)


def load(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f))


rows = load(sys.argv[1] if len(sys.argv) > 1 else 'results/mr_pairs.csv')

# 每个 pair 只取 ivw_re 作为保守口径 (异质性存在时更宽)
seen = {}
for r in rows:
    key = r['pair']
    m = r['method']
    if m == 'ivw_re':
        seen[key] = r
    elif m == 'ivw_fe' and key not in seen:
        seen[key] = r

print()
print("=" * 104)
print("最小可检出效应 (MDE) — alpha=0.05 双侧, power=80%, 基于观测 SE (取 IVW 随机效应, 更保守)")
print("=" * 104)
print("%-38s %5s %8s %10s  %-22s" % ("分析", "nSNP", "观测OR", "p", "80%功效可检出 OR 区间外"))
print("-" * 104)
out = []
for k, r in seen.items():
    se = float(r['se']); orv = float(r['OR']); p = float(r['pval']); n = r['nsnp_used']
    b, hi, lo = mde(se)
    print("%-38s %5s %8.3f %10.3f  OR >= %.3f  或  OR <= %.3f" % (k, n, orv, p, hi, lo))
    out.append((k, n, orv, float(r['LCI']), float(r['UCI']), p, hi, lo))
print("-" * 104)
print("解读: 若真实效应落在 [下限, 上限] 区间内, 本研究有 80% 以上把握漏检 (假阴性)。")
print("      区间越窄 -> 排除能力越强 -> 阴性结论越可靠。")
print()

with open('results/mde.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['pair', 'nsnp_used', 'OR', 'LCI', 'UCI', 'pval',
                'MDE_OR_upper', 'MDE_OR_lower'])
    for r in out:
        w.writerow(r)
print("已写出 results/mde.csv")
