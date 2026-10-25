# -*- coding: utf-8 -*-
"""按生信期刊（Nature Genetics / Genome Medicine）风格重画 Figure 1–4。

风格要点
- 低饱和配色、浅填充 + 深边框、黑字（弃用高饱和蓝绿橙粉）。
- 森林图用 meta 分析惯例：方块（面积 ∝ 1/SE²）+ 细横线 + 右侧 OR (95% CI) 数值列，
  不用等大圆点；正向 / 反向分组。
- 三栏布局：左侧标签轴 + 中间图轴 + 右侧数值轴，行距固定，无大片空白。
- Arial，最小字号 8 pt，figsize 即最终尺寸（170 mm 宽，不用 bbox_inches 以免缩宽），600 dpi。
"""
import json
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

plt.rcParams.update({'font.family': 'Arial', 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7})

OUT = 'figures_MR'
MM = 1 / 25.4
DPI = 600
SRC = 'A1_眼表营养轴MR_分析管线/results/figdata_wald_loo.json'
SRC2 = 'A1_眼表营养轴MR_分析管线/results/figdata_scatter.json'

C_BLUE = '#4C72B0'
C_GREEN = '#55A868'
C_ORANGE = '#C98A5B'
C_RED = '#B85C5C'
C_TEXT = '#222222'
C_AXIS = '#444444'
C_REF = '#9A9A9A'


def save(fig, name):
    # ---- 门禁：逐元素量溢出 + 字号下限（不能只靠肉眼看图）----
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    Wp, Hp = fig.get_size_inches() * fig.dpi
    # 只查用户添加的文本（轴标签 / ax.text / fig.text）。
    # tick 标签由轴管理且 clip_on=True，log 轴的 minor tick 会落在轴外但不会显示，不能算溢出。
    targets = list(fig.texts)
    for a in fig.axes:
        targets.extend(a.texts)
        targets.append(a.xaxis.label)
        targets.append(a.yaxis.label)
    bad = []
    for t in targets:
        if not t.get_text().strip():
            continue
        fs = t.get_fontsize()
        if fs < 8.0 - 1e-6:
            bad.append('字号 %.1f < 8: %r' % (fs, t.get_text()[:36]))
        bb = t.get_window_extent(r)
        over = max(-bb.x0, -bb.y0, bb.x1 - Wp, bb.y1 - Hp)
        if over > 0.5:
            bad.append('溢出画布 %.1f px: %r' % (over, t.get_text()[:36]))
        # 轴内文本不得越出所属 axes 的包围盒（防"被相邻轴覆盖/裁切"的隐性截断）。
        # xlabel/ylabel 豁免——居中超轴宽属正常，只有超出画布才被裁。
        own = None
        for a in fig.axes:
            if t in a.texts:
                own = a
                break
        if own is not None:
            ab = own.get_window_extent(r)
            pad = 1.0
            if (bb.x0 < ab.x0 - pad or bb.x1 > ab.x1 + pad
                    or bb.y0 < ab.y0 - pad or bb.y1 > ab.y1 + pad):
                bad.append('越出所属轴 %.1f px: %r' % (
                    max(ab.x0 - bb.x0, bb.x1 - ab.x1, ab.y0 - bb.y0,
                        bb.y1 - ab.y1), t.get_text()[:36]))
    for b in bad:
        print('    ⚠ ' + b)
    if not bad:
        print('    门禁 OK（字号 ≥8、无溢出、无轴内截断）')
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, '%s.%s' % (name, ext)), dpi=DPI,
                    facecolor='white')
    plt.close(fig)
    print('  已生成 %s' % name)


# ============================================================ Figure 1 流程图
def fig1():
    W, H = 170, 98
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H)
    ax.axis('off')

    def box(x, y, w, h, fill, edge):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge,
                               linewidth=0.9, zorder=2))

    def txt(x, y, s, size=8, weight='normal', ha='center'):
        ax.text(x, y, s, fontsize=size, fontweight=weight, color=C_TEXT,
                ha=ha, va='center', zorder=3)

    def arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='-|>', color=C_AXIS, lw=0.9,
                                    shrinkA=1, shrinkB=1, mutation_scale=8),
                    zorder=1)

    c1x, c1w = 3, 58
    c2x, c2w = 72, 40
    c3x, c3w = 123, 44
    for cx, cw, lab in ((c1x, c1w, 'GWAS summary statistics'),
                        (c2x, c2w, 'Instrument sets'),
                        (c3x, c3w, 'Two-sample MR')):
        txt(cx + cw / 2, 92, lab, 9, 'bold')

    rows = [
        ('#D9E5F1', C_BLUE, 'Dry eye disease (exposure)',
         '84,262 / 314,374 (MVP, EUR)\nPheCode 375.1', '39 loci\nF 21.7–45.4', 68, 15),
        ('#F4E7D7', C_ORANGE, 'Alzheimer disease',
         '85,934 / 401,577\n(Bellenguez 2022)', '43 loci', 48.5, 12.5),
        ('#DDEBDD', C_GREEN, 'Parkinson disease',
         '33,674 / 449,056\n(Nalls 2019)', '25 loci', 32, 12.5),
        ('#F2DEDF', C_RED, 'Multiple sclerosis',
         '14,498 / 24,091 (IMSGC 2013)\n2,926 / 495,931 (FinnGen R12)',
         '35 loci\nImmunoChip', 12, 16),
    ]
    mid = []
    for fill, edge, name, sub, ivt, y, h in rows:
        box(c1x, y, c1w, h, fill, edge)
        txt(c1x + 2.2, y + h - 4.3, name, 8, 'bold', ha='left')
        txt(c1x + 2.2, y + h / 2 - 4.6, sub, 8, ha='left')
        box(c2x, y, c2w, h, '#F5F5F5', '#8A8A8A')
        txt(c2x + c2w / 2, y + h / 2, ivt, 8)
        arrow(c1x + c1w + 1.5, y + h / 2, c2x - 1.5, y + h / 2)
        mid.append(y + h / 2)

    for title, sub, fill, edge, y, h in [
            ('Forward MR', 'DED → AD / PD / MS', '#E8EEF5', C_BLUE, 62, 13),
            ('Reverse MR', 'AD / PD / MS → DED', '#E8F0E9', C_GREEN, 42, 13),
            ('Sensitivity and robustness',
             'leave-one-out · pleiotropy\ncolocalization · minimum\ndetectable effect',
             '#F2F2F2', '#8A8A8A', 14, 17)]:
        box(c3x, y, c3w, h, fill, edge)
        txt(c3x + c3w / 2, y + h - 4.2, title, 8, 'bold')
        txt(c3x + c3w / 2, y + h / 2 - 4.8, sub, 8)

    arrow(c2x + c2w + 1.2, mid[0] + 1.2, c3x - 1.2, 68.5)
    arrow(c2x + c2w + 1.2, mid[0] - 1.2, c3x - 1.2, 48.5)
    for m in mid[1:]:
        arrow(c2x + c2w + 1.2, m + 1.2, c3x - 1.2, 48.5)
    for m in mid:
        arrow(c2x + c2w + 1.2, m - 1.4, c3x - 1.2, 22)

    # 用 mathtext 渲染上标：Arial 无 U+207B（⁻）字形，直接用会出方框
    txt(3, 4, 'Thresholds: P < 5 × 10$^{-7}$ (primary), 5 × 10$^{-8}$ (strict), '
        '5 × 10$^{-6}$ (relaxed); 29-locus multi-ancestry set.', 8, ha='left')
    save(fig, 'fig1_diagram')


# ============================================================ 通用森林图
def draw_forest(rows, xlim, xticks, xlabel, name, refline=1.0, label_w=0.30):
    """rows: [{'type':'group','label':..} | {'type':'row','label','or','lo','hi','color'}]

    紧凑布局：H 按行数算，axes 占满除 xlabel 外的空间，底部不留大片空白。
    """
    ys = []
    y = 0.0
    for r in rows:
        y += 1.0
        ys.append(y)
        if r['type'] == 'group':
            y += 0.4
    ymax = y + 0.6
    nrow = sum(1 for r in rows if r['type'] == 'row')
    ngrp = sum(1 for r in rows if r['type'] == 'group')
    XM = 6.0                                   # xlabel 占用的毫米数
    H = 13 + XM + 4.3 * nrow + 7.2 * ngrp      # 顶留白 13 mm

    bot_mm = XM + 5.0                          # ax 底距画布底：xlabel 之下还要给刻度标签留 5 mm
    bot = bot_mm / H
    hh = 0.99 - bot
    fig = plt.figure(figsize=(170 * MM, H * MM))
    axl = fig.add_axes([0.005, bot, label_w, hh])
    ax = fig.add_axes([label_w + 0.005, bot, 0.74 - label_w, hh])
    axr = fig.add_axes([0.755, bot, 0.24, hh])
    for a in (axl, axr):
        a.axis('off'); a.set_xlim(0, 1); a.set_ylim(ymax, 0)

    ax.set_xscale('log'); ax.set_xlim(*xlim); ax.set_ylim(ymax, 0)
    # log 轴会自动补刻度（并生成 2×10⁻¹ 之类落到轴外的标签），用 Fixed 精确定位
    from matplotlib.ticker import FixedLocator, FixedFormatter
    ax.xaxis.set_major_locator(FixedLocator(xticks))
    ax.xaxis.set_major_formatter(FixedFormatter(['%g' % t for t in xticks]))
    ax.tick_params(labelsize=8)
    ax.set_yticks([])
    fig.text(0.5, 3.0 / H, xlabel, fontsize=9,
             ha='center', va='center', color=C_TEXT)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis='x', colors=C_AXIS, length=2.5, which='major')
    ax.minorticks_off()
    ax.axvline(refline, color=C_REF, lw=0.7, ls=(0, (4, 3)), zorder=1)

    wts = [1 / ((math.log(r['hi']) - math.log(r['or'])) / 1.96) ** 2
           for r in rows if r['type'] == 'row']
    wmin, wmax = min(wts), max(wts)
    area = lambda w: 24 + 120 * (w - wmin) / max(wmax - wmin, 1e-12)

    wi = 0
    for r, yy in zip(rows, ys):
        if r['type'] == 'group':
            axl.text(0.02, yy, r['label'], fontsize=8.5, fontweight='bold',
                     ha='left', va='center', color=C_TEXT)
            continue
        o, lo, hi, col = r['or'], r['lo'], r['hi'], r['color']
        ax.plot([lo, hi], [yy, yy], color=col, lw=0.9, solid_capstyle='butt', zorder=2)
        for v in (lo, hi):
            ax.plot([v, v], [yy - 0.13, yy + 0.13], color=col, lw=0.9, zorder=2)
        ax.scatter([o], [yy], marker='s', s=area(wts[wi]), color=col,
                   edgecolors='white', linewidths=0.4, zorder=3)
        wi += 1
        axl.text(0.02, yy, r['label'], fontsize=8, ha='left', va='center', color=C_TEXT)
        axr.text(0.02, yy, '%.2f (%.2f–%.2f)' % (o, lo, hi), fontsize=8,
                 ha='left', va='center', color=C_TEXT)
    axr.text(0.02, 0.5, 'OR (95% CI)', fontsize=8, fontweight='bold',
             ha='left', va='center', color=C_TEXT)
    save(fig, name)


def fig2():
    rows = [{'type': 'group', 'label': 'Forward: DED → neurodegenerative disease'}]
    for lab, o, lo, hi in [('DED → Parkinson disease', 0.6198, 0.4681, 0.8207),
                           ('DED → Alzheimer disease', 1.0348, 0.9417, 1.1372),
                           ('DED → multiple sclerosis', 1.2815, 0.8987, 1.8274)]:
        rows.append({'type': 'row', 'label': lab, 'or': o, 'lo': lo, 'hi': hi,
                     'color': C_BLUE})
    rows.append({'type': 'group', 'label': 'Reverse: neurodegenerative disease → DED'})
    for lab, o, lo, hi in [('Parkinson disease → DED', 1.0039, 0.9849, 1.0232),
                           ('Alzheimer disease → DED', 0.9754, 0.9482, 1.0034),
                           # 2026-10-25 修正：MS 反向曾用带符号 bug 的 Python 版算成
                           # 1.0146 (0.9998–1.0295)；按 mr_core.R 的「只翻结局 beta」重算为
                           # OR 1.0062 (0.9916–1.0210)，P = 0.41（见 results/mr_pairs.csv）。
                           ('Multiple sclerosis → DED', 1.0062, 0.9916, 1.0210)]:
        rows.append({'type': 'row', 'label': lab, 'or': o, 'lo': lo, 'hi': hi,
                     'color': C_GREEN})
    draw_forest(rows, (0.42, 2.4), [0.5, 0.7, 1.0, 1.5, 2.0],
                'Odds ratio (log scale)', 'fig2_forest', label_w=0.45)


def fig3():
    d = json.load(open(SRC, encoding='utf-8'))
    w = sorted(d['wald'], key=lambda x: x['OR'])
    rows = [{'type': 'row', 'label': r['snp'], 'or': r['OR'], 'lo': r['lci'],
             'hi': r['uci'], 'color': C_BLUE} for r in w]
    draw_forest(rows, (0.045, 145), [0.05, 0.2, 1, 5, 25, 120],
                'Wald ratio per dry-eye-increasing allele (log scale)',
                'fig3_wald', label_w=0.22)


def fig4():
    d = json.load(open(SRC, encoding='utf-8'))
    w = sorted(d['loo'], key=lambda x: x['OR'])
    rows = [{'type': 'row', 'label': 'Omitting ' + r['snp'], 'or': r['OR'],
             'lo': r['lci'], 'hi': r['uci'], 'color': C_BLUE} for r in w]
    overall = 0.6198
    draw_forest(rows, (0.40, 0.86), [0.4, 0.5, 0.6, 0.7, 0.8],
                'Inverse-variance weighted estimate omitting each instrument in turn',
                'fig4_loo', refline=overall, label_w=0.30)


# ============================================================ Figure 5 散点图
# (key, 面板字母, 标题, x 轴标签, y 轴标签, 点色)
SCATTER = [
    ('DED_AD', 'A', 'Dry eye → Alzheimer disease',
     'β on dry eye disease', 'β on Alzheimer disease', C_BLUE),
    ('DED_PD', 'B', 'Dry eye → Parkinson disease',
     'β on dry eye disease', 'β on Parkinson disease', C_BLUE),
    ('DED_MS', 'C', 'Dry eye → multiple sclerosis',
     'β on dry eye disease', 'β on multiple sclerosis', C_BLUE),
    ('AD_DED', 'D', 'Alzheimer disease → dry eye',
     'β on Alzheimer disease', 'β on dry eye disease', C_GREEN),
    ('PD_DED', 'E', 'Parkinson disease → dry eye',
     'β on Parkinson disease', 'β on dry eye disease', C_GREEN),
    ('MS_DED', 'F', 'Multiple sclerosis → dry eye',
     'β on multiple sclerosis', 'β on dry eye disease', C_GREEN),
]

# 2 × 3 栅格的公共几何（毫米，y 自底向上）
_W5, _H5 = 170.0, 142.0
_LEFT5, _RIGHT5, _GAPX5 = 16.0, 4.0, 12.0
_PW5 = (_W5 - _LEFT5 - _RIGHT5 - 2 * _GAPX5) / 3.0
_XS5 = [_LEFT5, _LEFT5 + _PW5 + _GAPX5, _LEFT5 + 2 * (_PW5 + _GAPX5)]
_PH5 = 40.0
_ROWS5 = [(84.0, 124.0), (20.0, 60.0)]        # (bottom, top)：第 0 行在上
_TITLE5 = [(137.0, 131.0), (73.0, 67.0)]      # (标题行, n= 行)
_PH = 40.0


def _grid_axes(fig, W, H, i):
    """第 i 个面板（0..5）的 axes：行优先，第 0 行在下方。"""
    r, c = divmod(i, 3)
    x0 = _XS5[c]
    bot, top = _ROWS5[r]
    return fig.add_axes([x0 / W, bot / H, _PW5 / W, (top - bot) / H])


def fig5():
    d = json.load(open(SRC2, encoding='utf-8'))
    W, H = _W5, _H5
    fig = plt.figure(figsize=(W * MM, H * MM))

    for i, (key, letter, title, xl, yl, col) in enumerate(SCATTER):
        v = d[key]
        bx = np.array(v['bx']); sx = np.array(v['sx'])
        by = np.array(v['by']); sy = np.array(v['sy'])
        b_ivw = v['ivw']['b']; b_eg = v['egger']['b']
        b0 = v['egger_intercept']; b_wm = v['wm']['b']

        ax = _grid_axes(fig, W, H, i)
        ax.errorbar(bx, by, xerr=1.959964 * sx, yerr=1.959964 * sy, fmt='none',
                    ecolor='#BFBFBF', elinewidth=0.5, capsize=0, zorder=1)
        ax.scatter(bx, by, s=9, facecolor=col, edgecolor='white', linewidths=0.3,
                   zorder=3)

        xr = [min(bx - 1.96 * sx), max(bx + 1.96 * sx)]
        yr = [min(by - 1.96 * sy), max(by + 1.96 * sy)]
        px = 0.10 * (xr[1] - xr[0]); py = 0.10 * (yr[1] - yr[0])
        xr = [xr[0] - px, xr[1] + px]; yr = [yr[0] - py, yr[1] + py]
        gx = np.linspace(xr[0], xr[1], 60)
        for slope, inter, ls, lw, c in (
                (b_ivw, 0.0, '-', 1.1, '#333333'),
                (b_eg, b0, (0, (4, 2)), 0.9, C_RED),
                (b_wm, 0.0, (0, (1.2, 1.2)), 1.0, C_ORANGE)):
            gy = inter + slope * gx
            yr = [min(yr[0], gy.min()), max(yr[1], gy.max())]
        pad = 0.06 * (yr[1] - yr[0])
        ax.set_xlim(*xr); ax.set_ylim(yr[0] - pad, yr[1] + pad)

        for slope, inter, ls, lw, c in (
                (b_ivw, 0.0, '-', 1.1, '#333333'),
                (b_eg, b0, (0, (4, 2)), 0.9, C_RED),
                (b_wm, 0.0, (0, (1.2, 1.2)), 1.0, C_ORANGE)):
            ax.plot(gx, inter + slope * gx, color=c, ls=ls, lw=lw, zorder=2)

        ax.axhline(0, color='#D0D0D0', lw=0.6, zorder=0)
        ax.axvline(0, color='#D0D0D0', lw=0.6, zorder=0)
        ax.tick_params(labelsize=8, colors=C_AXIS, length=2.5)
        ax.set_xlabel(xl, fontsize=8, color=C_TEXT, labelpad=2)
        ax.set_ylabel(yl, fontsize=8, color=C_TEXT, labelpad=2)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for sp in ('bottom', 'left'):
            ax.spines[sp].set_linewidth(0.7)
            ax.spines[sp].set_color(C_AXIS)

        r, c = divmod(i, 3)
        ty, ny = _TITLE5[r]
        fig.text(_XS5[c] / W, ty / H, '%s  %s' % (letter, title), fontsize=8.5,
                 fontweight='bold', color=C_TEXT, ha='left', va='center')
        fig.text(_XS5[c] / W, ny / H, 'n = %d instruments' % v['nsnp'],
                 fontsize=8, color=C_AXIS, ha='left', va='center')

    # 图例（独立 axes，让门禁能校验其内文本）
    axL = fig.add_axes([0.10, 0.013, 0.80, 0.048]); axL.axis('off')
    axL.set_xlim(0, 1); axL.set_ylim(0, 1)
    for x0, lab, c, ls in ((0.02, 'Inverse-variance weighted', '#333333', '-'),
                           (0.34, 'MR-Egger', C_RED, (0, (4, 2))),
                           (0.50, 'Weighted median', C_ORANGE, (0, (1.2, 1.2)))):
        axL.plot([x0, x0 + 0.05], [0.5, 0.5], color=c, ls=ls, lw=1.0,
                 clip_on=False)
        axL.text(x0 + 0.055, 0.5, lab, fontsize=8, color=C_TEXT, ha='left',
                 va='center')
    save(fig, 'fig5_scatter')


# ============================================================ Figure S1 漏斗图
def figS1():
    d = json.load(open(SRC2, encoding='utf-8'))
    W, H = _W5, _H5
    fig = plt.figure(figsize=(W * MM, H * MM))

    for i, (key, letter, title, xl, yl, col) in enumerate(SCATTER):
        v = d[key]
        bx = np.array(v['bx']); sy = np.array(v['sy'])
        wald = np.array(v['wald'])
        prec = np.abs(bx) / sy              # 1 / SE(Wald ratio)
        b = v['ivw']['b']
        ax = _grid_axes(fig, W, H, i)
        ax.scatter(wald, prec, s=9, facecolor=col, edgecolor='white',
                   linewidths=0.3, zorder=3)
        ymax = prec.max() * 1.08
        gy = np.linspace(prec.min() * 0.9, ymax, 80)
        half = 1.959964 / gy
        ax.plot(b - half, gy, color=C_RED, ls=(0, (4, 2)), lw=0.8, zorder=2)
        ax.plot(b + half, gy, color=C_RED, ls=(0, (4, 2)), lw=0.8, zorder=2)
        ax.axvline(b, color='#333333', lw=1.0, zorder=2)
        xr = [min(wald.min(), (b - half).min()), max(wald.max(), (b + half).max())]
        px = 0.06 * (xr[1] - xr[0])
        ax.set_xlim(xr[0] - px, xr[1] + px)
        ax.set_ylim(0, ymax)
        ax.tick_params(labelsize=8, colors=C_AXIS, length=2.5)
        ax.set_xlabel('Wald ratio (log odds)', fontsize=8, color=C_TEXT,
                      labelpad=2)
        ax.set_ylabel('1 / SE', fontsize=8, color=C_TEXT, labelpad=2)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for sp in ('bottom', 'left'):
            ax.spines[sp].set_linewidth(0.7)
            ax.spines[sp].set_color(C_AXIS)

        r, c = divmod(i, 3)
        ty, ny = _TITLE5[r]
        fig.text(_XS5[c] / W, ty / H, '%s  %s' % (letter, title), fontsize=8.5,
                 fontweight='bold', color=C_TEXT, ha='left', va='center')
        fig.text(_XS5[c] / W, ny / H, 'n = %d instruments' % v['nsnp'],
                 fontsize=8, color=C_AXIS, ha='left', va='center')

    axL = fig.add_axes([0.10, 0.013, 0.80, 0.048]); axL.axis('off')
    axL.set_xlim(0, 1); axL.set_ylim(0, 1)
    axL.plot([0.02, 0.07], [0.5, 0.5], color='#333333', lw=1.0, clip_on=False)
    axL.text(0.075, 0.5, 'Inverse-variance weighted estimate', fontsize=8,
             color=C_TEXT, ha='left', va='center')
    axL.plot([0.48, 0.53], [0.5, 0.5], color=C_RED, ls=(0, (4, 2)), lw=0.8,
             clip_on=False)
    axL.text(0.535, 0.5, '1.96 × SE of the Wald ratio', fontsize=8,
             color=C_TEXT, ha='left', va='center')
    save(fig, 'figS1_funnel')


# ============================================================ Figure 6 机制 DAG
def fig6():
    W, H = 170.0, 148.0
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H)
    ax.axis('off')

    def box(x, y, w, h, fill, edge):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge,
                               linewidth=0.8, zorder=2))

    def txt(x, y, s, size=8, weight='normal', ha='center', color=C_TEXT):
        ax.text(x, y, s, fontsize=size, fontweight=weight, color=color,
                ha=ha, va='center', zorder=3)

    def arrow(x1, y1, x2, y2, dashed=False):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='-|>', color=C_AXIS, lw=0.9,
                                    ls=(0, (3, 2)) if dashed else '-',
                                    shrinkA=1, shrinkB=1, mutation_scale=8),
                    zorder=1)

    # ---------------- Panel A：观察性关联 ----------------
    txt(6, 143, 'A  Observational association between dry eye and '
        'neurodegenerative disease', 9, 'bold', ha='left')
    txt(6, 137.5, 'Age, sex and comorbidity are common causes; medication and '
        'ascertainment follow the diagnosis', 8, ha='left', color=C_AXIS)

    # 共同原因（真混杂）
    box(58, 122, 54, 10, '#F7F7F7', '#9A9A9A')
    txt(85, 127, 'Age, sex and systemic comorbidity', 8)
    arrow(68, 122, 36, 115, dashed=True)
    arrow(102, 122, 134, 115, dashed=True)

    box(6, 99, 60, 16, '#D9E5F1', C_BLUE)
    txt(36, 109.6, 'Neurodegenerative disease', 8.5, 'bold')
    txt(36, 104.4, 'Alzheimer, Parkinson, multiple sclerosis', 8)
    box(104, 99, 60, 16, '#F4E7D7', C_ORANGE)
    txt(134, 109.6, 'Dry eye disease', 8.5, 'bold')
    txt(134, 104.4, 'as recorded in health records', 8)

    ax.annotate('', xy=(104, 103), xytext=(66, 103),
                arrowprops=dict(arrowstyle='<|-|>', color=C_AXIS, lw=0.9,
                                ls=(0, (3, 2)), shrinkA=0, shrinkB=0,
                                mutation_scale=8), zorder=1)
    txt(85, 110, 'observed association', 8, ha='center', color=C_AXIS)

    # 诊断的下游后果（介导路径，非混杂）
    box(24, 76, 48, 12, '#F7F7F7', '#9A9A9A')
    txt(48, 82, 'Anticholinergic medication', 8)
    box(98, 76, 48, 12, '#F7F7F7', '#9A9A9A')
    txt(122, 84.4, 'Healthcare contact and', 8)
    txt(122, 79.6, 'diagnostic ascertainment', 8)
    arrow(36, 99, 48, 88, dashed=True)
    arrow(50, 99, 116, 88, dashed=True)
    arrow(66, 88, 120, 99, dashed=True)
    arrow(140, 88, 140, 99, dashed=True)

    # ---------------- 分隔 ----------------
    ax.plot([4, 166], [70, 70], color='#D8D8D8', lw=0.7, zorder=0)

    # ---------------- Panel B：孟德尔随机化 ----------------
    txt(6, 64, 'B  Mendelian randomization: what the genetic instruments '
        'do capture', 9, 'bold', ha='left')

    chains = [
        (47, 'Dry eye instruments', '39 variants, F 21.7–45.4',
         'Dry eye liability', 'AD · PD · MS liability', 'forward MR'),
        (25, 'Disease instruments', '25–43 variants',
         'AD · PD · MS liability', 'Dry eye liability', 'reverse MR'),
    ]
    xs3 = [9, 65, 121]
    for y0, g1, g2, m1, m2, lab in chains:
        yc = y0 + 7
        box(xs3[0], y0, 40, 14, '#F2F2F2', '#8A8A8A')
        txt(xs3[0] + 20, y0 + 9.4, g1, 8, 'bold')
        txt(xs3[0] + 20, y0 + 4.6, g2, 8)
        box(xs3[1], y0, 40, 14, '#D9E5F1', C_BLUE)
        txt(xs3[1] + 20, yc, m1, 8)
        box(xs3[2], y0, 40, 14, '#F4E7D7', C_ORANGE)
        txt(xs3[2] + 20, yc, m2, 8)
        arrow(xs3[0] + 40, yc, xs3[1], yc)
        arrow(xs3[1] + 40, yc, xs3[2], yc)
        txt(113, yc + 4.5, lab, 8, ha='center', color=C_AXIS)

    txt(6, 19, 'Germline variants are fixed at conception and independent of '
        'medication use, healthcare contact and disease onset,', 8, ha='left',
        color=C_AXIS)
    txt(6, 14.5, 'and act on the outcome only through the exposure, so they '
        'are not subject to reverse causation.', 8, ha='left', color=C_AXIS)

    ax.plot([4, 12], [6.5, 6.5], color=C_AXIS, lw=0.9, zorder=2)
    txt(14, 6.5, 'causal path captured by Mendelian randomization', 8,
        ha='left', color=C_AXIS)
    ax.plot([84, 92], [6.5, 6.5], color=C_AXIS, lw=0.9, ls=(0, (3, 2)),
            zorder=2)
    txt(94, 6.5, 'non-genetic path, not instrumented', 8, ha='left',
        color=C_AXIS)
    save(fig, 'fig6_mechanism')


if __name__ == '__main__':
    print('Figure 1:'); fig1()
    print('Figure 2:'); fig2()
    print('Figure 3:'); fig3()
    print('Figure 4:'); fig4()
    print('Figure 5:'); fig5()
    print('Figure 6:'); fig6()
    print('Figure S1:'); figS1()
