#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
修正 GWAS 汇总统计里的样本量列
------------------------------------------------
问题: GWAS Catalog harmonised 文件的 n_cases / n_controls 常是「逐 SNP 分批样本量」
      (PD 文件首行只有 1169/968), 远小于研究总数 (33674/449056)。
      Steiger 定向检验用 4*ncase*nctrl/(ncase+nctrl) 作有效样本量,
      样本量被低估 38 倍会让 r2 虚高, 产生"方向错误"的假警报。

解决: 按文件来源把这两列统一改成论文报告的总样本量。
"""
import csv, os, sys

# 论文报告的总样本量
TOTAL = {
    'DED_MVP_EUR': (84262, 314374),     # GCST90475895
    'DED_MVP_META': (124106, 438620),   # GCST90481119
    'AD_Bellenguez2022': (85934, 401577),
    'PD_Nalls2019': (33674, 449056),
    'MS_FinnGen_R12': (2926, 495931),
    # !! GCST005531 不是 IMSGC 2019 !!
    # FTP 文件名: imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv.gz
    # = IMSGC 2013 (PMID 24076602, Nat Genet) 的 ImmunoChip 子研究,
    #   discovery 14,498 cases / 24,091 controls, 仅 161,311 个常染色体变异 (免疫区定制芯片)。
    # 之前误标为 IMSGC 2019 (14,802/26,703) —— 那其实是该文的 replication 阶段数字。
    'MS_IMSGC2013': (14498, 24091),
}

# 文件名关键字 -> 数据集
KEYMAP = [
    ('out_AD', 'AD_Bellenguez2022'),
    ('out_PD', 'PD_Nalls2019'),
    ('out_MSFinn', 'MS_FinnGen_R12'),
    ('out_MSIMSGC', 'MS_IMSGC2013'),
    ('iv_MS', 'MS_IMSGC2013'),
    ('iv_AD', 'AD_Bellenguez2022'),
    ('iv_PD', 'PD_Nalls2019'),
    ('out_DED_for', 'DED_MVP_EUR'),
    ('iv_DED_meta', 'DED_MVP_META'),
    ('iv_DED_eur', 'DED_MVP_EUR'),
    ('iv_DED_', 'DED_MVP_EUR'),
    ('out_MS_for_DED', 'MS_IMSGC2013'),
]

D = sys.argv[1] if len(sys.argv) > 1 else 'data/processed'


def pick(name):
    for k, v in KEYMAP:
        if k in name:
            return v
    return None


n = 0
for fn in sorted(os.listdir(D)):
    if not fn.endswith('.tsv'):
        continue
    ds = pick(fn)
    if ds is None:
        continue
    path = os.path.join(D, fn)
    with open(path, encoding='utf-8') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    if not rows:
        continue
    hdr = rows[0]
    try:
        ic, ik = hdr.index('n_cases'), hdr.index('n_controls')
    except ValueError:
        # 缺列则追加
        hdr += ['n_cases', 'n_controls']
        ic, ik = len(hdr) - 2, len(hdr) - 1
        for r in rows[1:]:
            r += ['', '']
    c, k = TOTAL[ds]
    old = (rows[1][ic], rows[1][ik]) if len(rows) > 1 else ('-', '-')
    for r in rows[1:]:
        r[ic], r[ik] = str(c), str(k)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        csv.writer(f, delimiter='\t', lineterminator='\n').writerows(rows)
    print("%-46s -> %-18s %6d / %-7d (原 %s/%s)" % (fn, ds, c, k, old[0], old[1]))
    n += 1
print("\n已修正 %d 个文件" % n)
