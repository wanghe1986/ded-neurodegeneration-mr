# Bidirectional Mendelian randomization of dry eye disease and neurodegenerative disease

Analysis code for the manuscript *"Genetic liability to dry eye disease and neurodegenerative
disease: a bidirectional two-sample Mendelian randomization study"* (submitted to *Eye and Vision*).

The study tests, in both directions, whether genetic liability to dry eye disease is causally
related to Alzheimer disease, Parkinson disease, and multiple sclerosis.

**Canonical repository:** https://github.com/wanghe1986/ded-neurodegeneration-mr

## What is in this repository

| File | Purpose |
|---|---|
| `mr_core.R` | Core estimators, implemented from published formulae: allele harmonization, Wald ratio, inverse-variance weighting (fixed and random effects), MR-Egger regression, weighted median, Cochran *Q*, Steiger directionality |
| `_mr_prep.py` | Downloads and formats GWAS summary statistics; builds instrument sets at three *P* thresholds |
| `_mr_pair.R` | Runs one exposure–outcome pair through `mr_core.R` and writes estimates |
| `_loo.py`/`_loo.R` | Leave-one-out exclusion |
| `_wald.py` | Single-SNP Wald ratios |
| `_power.py` | Minimum detectable effect at 80% power |
| `_confound.py` | Phenome-wide pleiotropy screen against the GWAS Catalog |
| `04_coloc.R` | Colocalization (approximate Bayes factor, five hypotheses) |
| `_strict_mr.py` | Recomputes the strict (*P* < 5 × 10⁻⁸) instrument set |
| `_fig_data.py`, `_make_figdata.py` | Data behind Figures 1–4 |
| `_make_tables.py` | Builds Tables 1–6 as a Word document |
| `_fixn.py` | Replaces per-variant sample sizes with study-level totals (see below) |

Files beginning `00_`–`08_`, `_diag`, `_region`, `_forest`, `_probe_headers`, and `_build_exposure`
are earlier exploratory or intermediate scripts retained so that every reported number can be traced.

## Data sources

All inputs are public summary statistics; no individual-level data are used.

| Trait | Role | Source and accession | Cases / controls |
|---|---|---|---|
| Dry eye disease | Exposure (primary) | GWAS Catalog GCST90475895 (Million Veteran Program, PheCode 375.1, European ancestry) | 84,262 / 314,374 |
| Dry eye disease | Exposure (sensitivity) | GWAS Catalog GCST90481119 (Million Veteran Program, multi-ancestry meta-analysis) | 124,106 / 438,620 |
| Alzheimer disease | Outcome and reverse exposure | GWAS Catalog GCST90027158 (Bellenguez et al. 2022) | 85,934 / 401,577 |
| Parkinson disease | Outcome and reverse exposure | GWAS Catalog GCST009325 (Nalls et al. 2019) | 33,674 / 449,056 |
| Multiple sclerosis | Outcome | FinnGen release 12, endpoint G6_MS | 2,926 / 495,931 |
| Multiple sclerosis | Reverse exposure | GWAS Catalog GCST005531 (IMSGC 2013 ImmunoChip) | 14,498 / 24,091 |

Because the GWAS Catalog harmonized files are large, they are not redistributed here. Download them
from https://www.ebi.ac.uk/gwas/ (summary statistics) and https://www.finngen.fi/ (FinnGen release 12)
before running the pipeline.

## Environment

- Python 3.13 with `numpy`, `scipy`, `python-docx`
- R 4.x, base R only (`stats`); no external Mendelian randomization package is required, because all
  estimators are implemented in `mr_core.R`

Two environment-specific constraints shaped the code and are retained deliberately:

1. **Linkage disequilibrium pruning by physical distance.** No PLINK binary or 1000 Genomes reference
   panel was available, so instruments were pruned by retaining a SNP only if no already retained SNP
   lay within 10 Mb on the same chromosome. This is more conservative than *r*² < 0.001 over 1 Mb.
   Replace `_mr_prep.py`'s pruning step if a reference panel is available.
2. **Study-level sample sizes.** The `n_cases` and `n_controls` fields in GWAS Catalog harmonized
   files record per-variant rather than study-level counts. Left uncorrected, this deflates the
   effective sample size used in Steiger testing (for the Parkinson disease file, by a factor of about
   38) and produces spurious indications of reverse causation. `_fixn.py` substitutes the study-level
   totals listed above throughout.

## Reproducing the results

```bash
# 1  fetch and format summary statistics, build instrument sets
python _mr_prep.py iv      --gwas GCST90475895.h.tsv.gz ...
python _mr_prep.py outcome --gwas GCST90027158.h.tsv.gz ...

# 2  one exposure-outcome pair
Rscript _mr_pair.R --exposure iv_DED_eur_p5e7.tsv --outcome out_AD_for_iv_DED_eur_p5e7.tsv

# 3  sensitivity analyses
python _strict_mr.py        # strict P < 5e-8 instrument set
Rscript _loo.R               # leave-one-out
python _power.py             # minimum detectable effects
python _confound.py          # GWAS Catalog pleiotropy screen

# 4  colocalization
Rscript 04_coloc.R

# 5  figures and tables
python _fig_data.py && python _make_figdata.py
python _make_tables.py
```

Every estimator follows `mr_core.R::mr_estimate()`: inverse-variance weights are `1/sy²`; MR-Egger is
a weighted regression of outcome on exposure coefficients with *t*-distributed coefficient tests on
`n − 2` degrees of freedom; the weighted median uses weights `w·bx²` with a bootstrap over SNPs.

## Licence

MIT. See `LICENSE`.
