# -*- coding: utf-8 -*-
"""生成 Eye and Vision 投稿稿件的 Table 1-6（Word，随后转三线表）。

数据来源：
  Table 1  稿件 Methods - Data sources（含 GCST accession）
  Table 2  results/mr_pairs.csv + results/strict_p5e8_PD.csv
  Table 3  results/mr_pairs.csv（反向三组）
  Table 4  results/confound_iv25.txt（25 工具 vs GWAS Catalog 1,203,840 条）
  Table 5  results/coloc_results.csv（九位点 PP.H0-H4）
  Table 6  results/mde.csv（80% power 最小可检出效应）
"""
import csv
import os

from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn

BASE = r'D:\workbuddy working'
RES = os.path.join(BASE, 'A1_眼表营养轴MR_分析管线', 'results')
OUT = os.path.join(BASE, 'Tables_1-6_干眼神经退行MR.docx')

SUP = str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻')


def sup(x):
    return str(x).translate(SUP)


def fp(p):
    """P 值格式化：>=0.001 三位小数，否则科学计数（Unicode 上标）。"""
    if p is None or p == '' or p == '—':
        return '—'
    p = float(p)
    if p < 0.001:
        n = 0
        m = p
        while m < 1:
            m *= 10
            n += 1
        return '%.1f × 10%s' % (m, sup(-n))
    return '%.3f' % p


def fci(OR, lo, hi):
    return '%.2f (%.2f–%.2f)' % (float(OR), float(lo), float(hi))


# ---------------------------------------------------------------- 读数据
def load_pairs():
    d = {}
    with open(os.path.join(RES, 'mr_pairs.csv'), encoding='utf-8') as f:
        for r in csv.DictReader(f):
            d.setdefault(r['pair'], {}).setdefault(r['method'], r)
    return d


def load_strict():
    d = {}
    with open(os.path.join(RES, 'strict_p5e8_PD.csv'), encoding='utf-8') as f:
        for r in csv.DictReader(f):
            d[r['method']] = r
    return d


PAIRS = load_pairs()
STRICT = load_strict()

METHODS = [('ivw_fe', 'IVW (fixed effects)'),
           ('ivw_re', 'IVW (random effects)'),
           ('mr_egger', 'MR-Egger'),
           ('weighted_median', 'Weighted median')]

# ---------------------------------------------------------------- 文档
doc = Document()
st = doc.styles['Normal']
st.font.name = 'Arial'
st.font.size = Pt(9)
st.element.rPr.rFonts.set(qn('w:eastAsia'), 'Arial')
sec = doc.sections[0]
sec.left_margin = sec.right_margin = Inches(0.8)
USABLE = 6.9


def caption(text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(14)
    p.paragraph_format.space_after = Pt(4)
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(9)
    r.font.name = 'Arial'
    return p


def footnote(text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(3)
    p.paragraph_format.space_after = Pt(10)
    r = p.add_run(text)
    r.font.size = Pt(8)
    r.font.name = 'Arial'
    r.italic = False
    return p


def add_table(header, rows, widths, align=None):
    """widths: 相对比例列表。"""
    t = doc.add_table(rows=1, cols=len(header))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    total = sum(widths)
    for j, h in enumerate(header):
        c = t.rows[0].cells[j]
        c.text = ''
        p = c.paragraphs[0]
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.space_before = Pt(2)
        r = p.add_run(h)
        r.bold = True
        r.font.size = Pt(8.5)
        r.font.name = 'Arial'
        if align and align[j] == 'l':
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
        else:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for row in rows:
        cells = t.add_row().cells
        for j, v in enumerate(row):
            cells[j].text = ''
            p = cells[j].paragraphs[0]
            p.paragraph_format.space_after = Pt(1)
            p.paragraph_format.space_before = Pt(1)
            r = p.add_run(str(v))
            r.font.size = Pt(8.5)
            r.font.name = 'Arial'
            if align and align[j] == 'l':
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            else:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for row in t.rows:
        for j, c in enumerate(row.cells):
            c.width = Inches(USABLE * widths[j] / total)
    return t


# ================================================================ Table 1
caption('Table 1. Data sources and sample sizes')
T1 = [
    ['Dry eye disease', 'Exposure (primary)',
     'Million Veteran Program phenome-wide scan, PheCode 375.1, GCST90475895 (ref. 12)',
     '84,262', '314,374', 'European'],
    ['Dry eye disease', 'Exposure (sensitivity)',
     'Million Veteran Program phenome-wide scan, multi-ancestry meta-analysis, GCST90481119 (ref. 12)',
     '124,106', '438,620', 'Multi-ancestry'],
    ['Alzheimer disease', 'Outcome and reverse exposure', 'Bellenguez et al. 2022, GCST90027158 (ref. 13)',
     '85,934', '401,577', 'European'],
    ['Parkinson disease', 'Outcome and reverse exposure', 'Nalls et al. 2019, GCST009325 (ref. 14)',
     '33,674', '449,056', 'European'],
    ['Multiple sclerosis', 'Outcome', 'FinnGen release 12, G6_MS (ref. 15)',
     '2,926', '495,931', 'European (Finnish)'],
    ['Multiple sclerosis', 'Reverse exposure', 'IMSGC 2013 ImmunoChip, GCST005531 (ref. 16)',
     '14,498', '24,091', 'European'],
]
add_table(['Trait', 'Role in the analysis', 'Source and accession', 'Cases', 'Controls', 'Ancestry'],
          T1, [1.4, 1.5, 3.0, 0.8, 0.9, 1.1], align=['l', 'l', 'l', 'c', 'c', 'l'])
footnote('IMSGC, International Multiple Sclerosis Genetics Consortium; MVP, Million Veteran Program. '
         'Reference numbers are those of the manuscript reference list. Case and control counts are the '
         'study-level totals reported with each GWAS Catalog study record; they replace the per-variant '
         'counts carried in the harmonized files, which would otherwise deflate the effective sample '
         'size used in Steiger testing. Both dry eye datasets come from the same MVP phenome-wide scan '
         '(reference 12); the multi-ancestry total is the sum of the European, African American, '
         'Hispanic, and East Asian strata. A larger MVP dry eye GWAS based on a chart-validated '
         'algorithm (132,637 cases, 352,201 controls; reference 7) has since been published but '
         'releases summary statistics only through controlled access, so it was not available. The '
         'multiple sclerosis outcome (FinnGen release 12) and the multiple sclerosis instrument set '
         '(IMSGC 2013 ImmunoChip) are drawn from different studies, so no sample overlap arises in the '
         'reverse-direction analysis.')

# ================================================================ Table 2
caption('Table 2. Forward-direction Mendelian randomization estimates: genetic liability to dry eye '
        'disease and risk of neurodegenerative disease')
T2 = []
SETS = [
    ('Primary (P < 5 × 10⁻⁷, 39 loci)', [
        ('DED → Alzheimer disease', 'DED_EUR39 -> AD'),
        ('DED → Parkinson disease', 'DED_EUR39 -> PD'),
        ('DED → multiple sclerosis', 'DED_EUR39 -> MS_FinnGen')]),
    ('Strict (P < 5 × 10⁻⁸, 9 loci)', [
        ('DED → Parkinson disease', '__STRICT__')]),
    ('Relaxed (P < 5 × 10⁻⁶, 130 loci)', [
        ('DED → Alzheimer disease', 'DED_EUR130 -> AD'),
        ('DED → Parkinson disease', 'DED_EUR130 -> PD'),
        ('DED → multiple sclerosis', 'DED_EUR130 -> MS_FinnGen')]),
    ('Multi-ancestry (29 loci)', [
        ('DED → Alzheimer disease', 'DED_meta29 -> AD'),
        ('DED → Parkinson disease', 'DED_meta29 -> PD'),
        ('DED → multiple sclerosis', 'DED_meta29 -> MS_FinnGen')]),
    ('Sensitivity: excluding SLC39A8', [
        ('DED → Parkinson disease', 'DED_EUR38 -> PD [剔SLC39A8]')]),
    ('Sensitivity: excluding 3 pleiotropic loci', [
        ('DED → Parkinson disease', 'DED_EUR36 -> PD [剔3个多效性位点]')]),
]
for setname, pairs in SETS:
    for label, key in pairs:
        src = STRICT if key == '__STRICT__' else PAIRS[key]
        for mkey, mlabel in METHODS:
            r = src.get(mkey)
            if r is None:
                continue
            eg_p = r['egger_pval'] if r['egger_pval'] else '—'
            nsnp = r.get('nsnp_used') or r.get('nsnp')
            T2.append([setname, label, mlabel, nsnp,
                       fci(r['OR'], r['LCI'], r['UCI']), fp(r['pval']),
                       fp(r['Q_pval']), fp(eg_p)])
setattr(doc, '_t2', T2)
add_table(['Instrument set', 'Direction', 'Method', 'SNPs', 'OR (95% CI)', 'P',
           'Cochran Q P', 'MR-Egger intercept P'],
          T2, [1.9, 1.5, 1.3, 0.5, 1.3, 0.8, 0.8, 0.9],
          align=['l', 'l', 'l', 'c', 'c', 'c', 'c', 'c'])
footnote('OR, odds ratio per unit increase in log odds of dry eye disease; CI, confidence interval. '
         'Cochran Q P and MR-Egger intercept P are repeated within each instrument set because they are '
         'properties of the set of instruments rather than of the estimator. The MR-Egger intercept tests '
         'directional pleiotropy; MR-Egger estimates are not used as primary causal estimates because the '
         'InSIDE assumption cannot be verified. DED, dry eye disease; IVW, inverse-variance weighted.')

# ================================================================ Table 3
caption('Table 3. Reverse-direction Mendelian randomization estimates: genetic liability to '
        'neurodegenerative disease and risk of dry eye disease')
T3 = []
for label, key in [('Alzheimer disease → DED', 'AD -> DED  [反向]'),
                   ('Parkinson disease → DED', 'PD -> DED  [反向]'),
                   ('Multiple sclerosis → DED', 'MS_IMSGC -> DED  [反向]')]:
    for mkey, mlabel in METHODS:
        r = PAIRS[key][mkey]
        T3.append([label, mlabel, r['nsnp_used'], fci(r['OR'], r['LCI'], r['UCI']),
                   fp(r['pval']), fp(r['Q_pval']), fp(r['egger_pval'] or '—')])
add_table(['Direction', 'Method', 'SNPs', 'OR (95% CI)', 'P', 'Cochran Q P', 'MR-Egger intercept P'],
          T3, [1.8, 1.4, 0.5, 1.4, 0.8, 0.9, 1.0],
          align=['l', 'l', 'c', 'c', 'c', 'c', 'c'])
footnote('OR, odds ratio per unit increase in log odds of the neurodegenerative disease. '
         'All estimates use non-overlapping exposure and outcome samples. DED, dry eye disease.')

# ================================================================ Table 4
caption('Table 4. Phenome-wide pleiotropy screen of the 25 dry eye instruments used in the primary '
        'dry eye to Parkinson disease analysis')
T4 = [
    ['rs13107325', 'SLC39A8', '552',
     'Lipids (172), body mass index (24), alcohol use (19), blood pressure (16), height (11), '
     'educational attainment (10), sleep apnoea (6), urate and fatty-acid metabolites (6), '
     'depression (3), schizophrenia (3), type 2 diabetes (3), allergy and asthma (3), smoking (3), '
     'physical activity (3), social factors (1), coronary artery disease (1)',
     'Removed in sensitivity analysis'],
    ['rs2305280', 'MADD', '5',
     'HDL lipid composition traits, systolic blood pressure (MTAG)',
     'Removed in sensitivity analysis'],
    ['rs2280381', 'LINC02132 – LINC01082', '2',
     'Systemic lupus erythematosus (MTAG), height',
     'Removed in sensitivity analysis'],
    ['rs10828249', '—', '4',
     'Albumin levels, breast cancer, endometriosis — no trait on a candidate confounding pathway',
     'Retained'],
    ['rs62244186', '—', '1', 'Depressive symptoms — not a candidate confounding pathway', 'Retained'],
    ['rs9401452', '—', '2',
     'Attention deficit hyperactivity disorder / autism spectrum disorder / intelligence (pleiotropy), '
     'intelligence — not a candidate confounding pathway', 'Retained'],
    ['19 further instruments', '—', '0',
     'No previous genome-wide significant association with any trait recorded in the GWAS Catalog',
     'Retained'],
]
add_table(['Instrument', 'Nearest gene', 'GWAS Catalog associations',
           'Traits with a recorded association', 'Action taken'],
          T4, [1.1, 1.1, 0.8, 3.2, 1.2], align=['l', 'l', 'c', 'l', 'l'])
footnote('Each of the 25 instruments was queried against 1,203,840 association records in the GWAS '
         'Catalog (accessed 24 October 2025). Counts are numbers of distinct association records at '
         'P < 5 × 10⁻⁸, so a single trait measured in several cohorts contributes several records. '
         'Parenthesised counts are records assigned to each prespecified confounding category; for '
         'rs13107325 these account for 284 of the 552 records, the remainder falling in traits outside '
         'the prespecified categories. Removing rs13107325 alone gave OR 0.62 (95% CI 0.47–0.83, '
         'P = 0.0014); additionally removing rs2305280 and rs2280381 gave OR 0.64 (95% CI 0.47–0.87, '
         'P = 0.0039). HDL, high-density lipoprotein; MTAG, multi-trait analysis of GWAS.')

# ================================================================ Table 5
caption('Table 5. Colocalization of dry eye disease and Parkinson disease signals at nine loci')
COLOC = [
    ('rs10828249', '10:20,535,798–22,535,798', '7,295', 0.13582, 0.86175, 0.00001, 0.00000, 0.00242),
    ('rs11690184', '2:59,302,564–61,302,564', '7,412', 0.99804, 0.00182, 0.00007, 0.00000, 0.00007),
    ('rs76850872', '3:31,348,964–33,348,964', '9,284', 0.99981, 0.00013, 0.00005, 0.00000, 0.00001),
    ('rs117589607', '8:62,198,361–64,198,361', '7,870', 0.99963, 0.00029, 0.00006, 0.00000, 0.00001),
    ('rs9401452', '6:97,504,804–98,504,804', '3,090', 0.99711, 0.00275, 0.00005, 0.00000, 0.00008),
    ('rs139927001', '11:94,687,023–95,687,023', '3,695', 0.99133, 0.00752, 0.00007, 0.00000, 0.00107),
    ('rs74734235', '3:65,190,264–66,190,264', '4,689', 0.99977, 0.00017, 0.00006, 0.00000, 0.00001),
    ('rs2280381', '16:85,485,027–86,485,027', '5,895', 0.99977, 0.00016, 0.00006, 0.00000, 0.00000),
    ('rs62244186', '3:44,198,291–45,198,291', '3,860', 0.98895, 0.01086, 0.00010, 0.00000, 0.00009),
]
T5 = [[s, reg, n, '%.5f' % a, '%.5f' % b, '%.5f' % c, '0.00000', '%.5f' % e]
      for s, reg, n, a, b, c, d, e in COLOC]
add_table(['Lead SNP', 'Region (GRCh37)', 'SNPs in region', 'PP.H0', 'PP.H1', 'PP.H2', 'PP.H3', 'PP.H4'],
          T5, [1.1, 1.7, 0.9, 0.75, 0.75, 0.75, 0.75, 0.75],
          align=['l', 'l', 'c', 'c', 'c', 'c', 'c', 'c'])
footnote('PP.H0, posterior probability of no association with either trait; PP.H1, association with dry '
         'eye disease only; PP.H2, association with Parkinson disease only; PP.H3, association with both '
         'traits through two distinct variants; PP.H4, association with both traits through a shared '
         'variant. Priors were p1 = p2 = 1 × 10⁻⁴ and p12 = 5 × 10⁻⁶. Posterior mass concentrates on '
         'PP.H1 at every locus, which is the expected pattern when the outcome association is too weak '
         'to be resolved: Parkinson disease z statistics at the candidate variants ranged from 1.6 to 2.9. '
         'Colocalization is therefore uninformative rather than negative.')

# ================================================================ Table 6
caption('Table 6. Minimum detectable effects at 80% power and two-sided α = 0.05')
T6 = []
NAMES = {
    'DED_EUR39 -> AD': 'DED → Alzheimer disease (39 loci)',
    'DED_EUR39 -> PD': 'DED → Parkinson disease (39 loci)',
    'DED_EUR39 -> MS_FinnGen': 'DED → multiple sclerosis (39 loci)',
    'DED_EUR130 -> AD': 'DED → Alzheimer disease (130 loci)',
    'DED_EUR130 -> PD': 'DED → Parkinson disease (130 loci)',
    'DED_EUR130 -> MS_FinnGen': 'DED → multiple sclerosis (130 loci)',
    'DED_meta29 -> AD': 'DED → Alzheimer disease (multi-ancestry)',
    'DED_meta29 -> PD': 'DED → Parkinson disease (multi-ancestry)',
    'DED_meta29 -> MS_FinnGen': 'DED → multiple sclerosis (multi-ancestry)',
    'AD -> DED  [反向]': 'Alzheimer disease → DED',
    'PD -> DED  [反向]': 'Parkinson disease → DED',
    'MS_IMSGC -> DED  [反向]': 'Multiple sclerosis → DED',
}
with open(os.path.join(RES, 'mde.csv'), encoding='utf-8') as f:
    for r in csv.DictReader(f):
        T6.append([NAMES.get(r['pair'], r['pair']), r['nsnp_used'],
                   fci(r['OR'], r['LCI'], r['UCI']), fp(r['pval']),
                   '%.3f' % float(r['MDE_OR_upper']), '%.3f' % float(r['MDE_OR_lower'])])
add_table(['Analysis', 'SNPs', 'Observed OR (95% CI)', 'Observed P',
           'MDE: OR above which detectable', 'MDE: OR below which detectable'],
          T6, [2.3, 0.5, 1.4, 0.8, 1.2, 1.2], align=['l', 'c', 'c', 'c', 'c', 'c'])
footnote('The minimum detectable effect is |β| = 2.8016 × SE, the value at which a two-sided test at '
         'α = 0.05 has 80% power, exponentiated. Observed estimates are from random-effects '
         'inverse-variance weighting, the estimator from which the standard error and hence the '
         'minimum detectable effect were derived. OR, odds ratio; CI, confidence interval; MDE, '
         'minimum detectable effect; DED, dry eye disease.')

doc.save(OUT)
print('已保存:', OUT)
for i, t in enumerate(doc.tables, 1):
    print('  Table %d: %d 行 × %d 列' % (i, len(t.rows), len(t.columns)))
