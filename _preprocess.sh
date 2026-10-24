#!/usr/bin/env bash
# A1 预处理: 流式解压 + 只保留目标基因/SNP 的行
# 目的: 避免把 2-3 GB 的原始文件全量展开, 只留下几十 MB 的有用数据
set -e
cd "D:/workbuddy working/A1_眼表营养轴MR_分析管线"

PY="C:/Users/wangh/.workbuddy/binaries/python/versions/3.13.12/python.exe"
PROC="data/processed"
mkdir -p "$PROC"

echo "=== [1/5] 生成目标基因列表 ==="
"$PY" -c "
import csv
genes=set()
with open('$PROC/gene_table.csv', encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        genes.add(r['gene'])
for g in sorted(genes): print(g)
" > "$PROC/gene_list.txt"
wc -l < "$PROC/gene_list.txt" | xargs echo "  目标基因数:"

echo "=== [2/5] 过滤 eQTLGen (只保留目标基因) ==="
# 列序: 1 Pvalue 2 SNP 3 SNPChr 4 SNPPos 5 AssessedAllele 6 OtherAllele
#       7 Zscore 8 Gene 9 GeneSymbol 10 GeneChr 11 GenePos 12 NrCohorts 13 NrSamples
gzip -dc data/raw/eqtlgen_cis_fdr005.txt.gz | \
awk -F'\t' -v OFS='\t' '
  NR==FNR { g[$1]=1; next }
  FNR==1  { print "SNP","GeneSymbol","AssessedAllele","OtherAllele","Zscore","Pvalue","NrSamples"; next }
  ($9 in g) { print $2,$9,$5,$6,$7,$1,$13 }
' "$PROC/gene_list.txt" - > "$PROC/eqtl_target.tsv"
wc -l < "$PROC/eqtl_target.tsv" | xargs echo "  保留行数:"

echo "=== [3/5] 提取 eQTL 涉及的 SNP 列表 ==="
awk -F'\t' 'NR>1 {print $1}' "$PROC/eqtl_target.tsv" | sort -u > "$PROC/snp_list.txt"
wc -l < "$PROC/snp_list.txt" | xargs echo "  SNP 数:"

echo "=== [4/5] 过滤 MAF (只保留这些 SNP) ==="
# 列序: 1 SNP 2 hg19_chr 3 hg19_pos 4 AlleleA 5 AlleleB 6 allA_total
#       7 allAB_total 8 allB_total 9 AlleleB_all
gzip -dc data/raw/eqtlgen_maf.txt.gz | \
awk -F'\t' -v OFS='\t' '
  NR==FNR { s[$1]=1; next }
  FNR==1  { print "SNP","AlleleA","AlleleB","AlleleB_freq"; next }
  ($1 in s) { print $1,$4,$5,$9 }
' "$PROC/snp_list.txt" - > "$PROC/maf_target.tsv"
wc -l < "$PROC/maf_target.tsv" | xargs echo "  保留行数:"

echo "=== [5/5] 完成 ==="
ls -la "$PROC"/eqtl_target.tsv "$PROC"/maf_target.tsv | awk '{printf "  %-40s %.1f MB\n", $9, $5/1048576}'
