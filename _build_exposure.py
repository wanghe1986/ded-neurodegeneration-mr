# -*- coding: utf-8 -*-
"""合并 eQTL 与 MAF, 换算 beta/se, 构建工具变量

eQTLGen 只给 Zscore, 换算公式:
    SE   = 1 / sqrt(2 * f * (1 - f) * N)     f = 效应等位基因(AssessedAllele)频率
    beta = Zscore * SE
    F    = (beta/SE)^2 = Zscore^2
"""
import csv, math, os

PROC = "data/processed"
P_THRESH = 5e-8
F_MIN = 10

# ---------- 1. 读 MAF, 建 SNP -> (AlleleA, AlleleB, freqB) ----------
maf = {}
with open(os.path.join(PROC, "maf_target.tsv"), encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        try:
            maf[row["SNP"]] = (row["AlleleA"].upper(), row["AlleleB"].upper(),
                                float(row["AlleleB_freq"]))
        except (ValueError, KeyError):
            continue
print(f"MAF 表: {len(maf)} 个 SNP")

# ---------- 2. 读 eQTL, 对齐等位基因, 换算 beta ----------
rows = []
no_maf = 0
allele_mismatch = 0
with open(os.path.join(PROC, "eqtl_target.tsv"), encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        snp = row["SNP"]
        gene = row["GeneSymbol"]
        try:
            pval = float(row["Pvalue"]); z = float(row["Zscore"])
            n = float(row["NrSamples"])
        except (ValueError, KeyError):
            continue
        if pval >= P_THRESH:
            continue
        ea = row["AssessedAllele"].upper(); oa = row["OtherAllele"].upper()

        m = maf.get(snp)
        if m is None:
            no_maf += 1
            continue
        aA, aB, freqB = m
        # 效应等位基因频率 (AssessedAllele 的频率)
        if ea == aB:
            fq = freqB
        elif ea == aA:
            fq = 1.0 - freqB
        else:
            allele_mismatch += 1
            continue
        if not (0.0 < fq < 1.0):
            allele_mismatch += 1
            continue

        se = 1.0 / math.sqrt(2.0 * fq * (1.0 - fq) * n)
        beta = z * se
        Fstat = z * z          # 等价于 (beta/se)^2
        if Fstat <= F_MIN:
            continue
        rows.append(dict(SNP=snp, gene=gene, effect_allele=ea, other_allele=oa,
                         beta=beta, se=se, pval=pval, F_stat=Fstat, n=int(n),
                         freq=fq, zscore=z))

print(f"缺 MAF 被剔除: {no_maf}")
print(f"等位基因无法对齐被剔除: {allele_mismatch}")
print(f"通过 P<{P_THRESH} 且 F>{F_MIN}: {len(rows)} 条")

# ---------- 3. 写出全部候选 IV ----------
out_all = os.path.join(PROC, "exposure_all_sig.csv")
with open(out_all, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"已写出 {out_all}")

# ---------- 4. 每基因取 top SNP (主分析用, 单 IV 无 LD 问题) ----------
best = {}
for r in rows:
    g = r["gene"]
    if g not in best or r["pval"] < best[g]["pval"]:
        best[g] = r
top = sorted(best.values(), key=lambda x: x["pval"])
out_top = os.path.join(PROC, "exposure_top_snp.csv")
with open(out_top, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(top)
print(f"\n每基因 top SNP: {len(top)} 个基因")

# ---------- 5. 输出供 AD GWAS 过滤用的 SNP 列表 ----------
snps = sorted({r["SNP"] for r in rows})
with open(os.path.join(PROC, "snps_for_ad.txt"), "w", encoding="utf-8") as f:
    for s in snps:
        f.write(s + "\n")
print(f"待匹配 AD 的 SNP 数: {len(snps)}")

# ---------- 6. 报告 ----------
print("\n=== 各基因可用 IV 情况 (top SNP) ===")
by_gene = {}
for r in rows:
    by_gene.setdefault(r["gene"], []).append(r)
for g, lst in sorted(by_gene.items(), key=lambda kv: -len(kv[1])):
    t = best[g]
    print(f"  {g:<12} 显著 SNP 数={len(lst):<4} top={t['SNP']:<14} "
          f"P={t['pval']:.2e}  Z={t['zscore']:+.2f}  F={t['F_stat']:.0f}")
