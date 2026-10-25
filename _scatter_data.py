# -*- coding: utf-8 -*-
"""导出六对比较的 SNP 级协调数据（供 Figure 5 散点图 / Figure S1 漏斗图使用）。

协调规则与 mr_core.R 一致：链翻转只翻结局 beta。
"""
import csv
import json
import math
import os

import numpy as np

sys_path = os.path.dirname(os.path.abspath(__file__))
import _mr_pair as M

PROC = os.path.join(sys_path, 'data', 'processed')

# (key, 面板标题, iv文件, out文件, 方向 forward/reverse)
PAIRS = [
    ('DED_AD', 'Dry eye → Alzheimer disease',
     'iv_DED_eur_p5e7.tsv', 'out_AD_for_iv_DED_eur_p5e7.tsv', 'forward'),
    ('DED_PD', 'Dry eye → Parkinson disease',
     'iv_DED_eur_p5e7.tsv', 'out_PD_for_iv_DED_eur_p5e7.tsv', 'forward'),
    ('DED_MS', 'Dry eye → multiple sclerosis',
     'iv_DED_eur_p5e7.tsv', 'out_MSFinn_for_iv_DED_eur_p5e7.tsv', 'forward'),
    ('AD_DED', 'Alzheimer disease → dry eye',
     'iv_AD.tsv', 'out_DED_for_AD.tsv', 'reverse'),
    ('PD_DED', 'Parkinson disease → dry eye',
     'iv_PD.tsv', 'out_DED_for_PD.tsv', 'reverse'),
    ('MS_DED', 'Multiple sclerosis → dry eye',
     'iv_MS.tsv', 'out_DED_for_MS.tsv', 'reverse'),
]


def main():
    out = {}
    for key, title, ivf, outf, direction in PAIRS:
        iv = M.rd(os.path.join(PROC, ivf))
        ot = M.rd(os.path.join(PROC, outf))
        bx, sx, by, sy, names, cnt = M.harmonize(iv, ot)
        rows, Q, df, Q_p, b0, se_b0, p_b0 = M.estimate(bx, sx, by, sy)
        est = {r['method']: r for r in rows}
        bx = np.array(bx); sx = np.array(sx)
        by = np.array(by); sy = np.array(sy)
        wald = by / bx
        out[key] = dict(
            title=title, direction=direction, iv_file=ivf, out_file=outf,
            nsnp=len(bx),
            snp=names,
            bx=bx.tolist(), sx=sx.tolist(),
            by=by.tolist(), sy=sy.tolist(),
            wald=wald.tolist(),
            ivw=est['ivw_fe'], ivw_re=est['ivw_re'],
            egger=est['mr_egger'], wm=est['weighted_median'],
            egger_intercept=b0, egger_intercept_se=se_b0, egger_intercept_p=p_b0,
            Q=Q, Q_df=df, Q_p=Q_p,
            harmonize={k: v for k, v in cnt.items() if v},
        )
        print('%-8s n=%2d  ivw OR %.4f (%.4f-%.4f) P=%.4g  egger OR %.4f P=%.4g  wm OR %.4f P=%.4g'
              % (key, len(bx), est['ivw_fe']['OR'], est['ivw_fe']['LCI'],
                 est['ivw_fe']['UCI'], est['ivw_fe']['pval'],
                 est['mr_egger']['OR'], est['mr_egger']['pval'],
                 est['weighted_median']['OR'], est['weighted_median']['pval']))
        print('         bx [%.4g, %.4g]  by [%.4g, %.4g]  wald [%.3g, %.3g]  1/sy [%.3g, %.3g]'
              % (bx.min(), bx.max(), by.min(), by.max(),
                 wald.min(), wald.max(), (1 / sy).min(), (1 / sy).max()))
        print('         协调:', {k: v for k, v in cnt.items() if v})

    dst = os.path.join(sys_path, 'results', 'figdata_scatter.json')
    with open(dst, 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print('\n已写出 %s' % dst)


if __name__ == '__main__':
    main()
