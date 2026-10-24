# 方案 B · 共定位分析 (coloc, 单因果变异假设)
#
# 目的: 对 MR 显著的位点, 判断干眼与神经退行是否共享同一个因果变异
#       PP.H4 > 0.8 视为强共定位证据
#
# 输入: 由 _region.py 生成的两个 GWAS 的区域文件
#       列: SNP pos effect_allele other_allele beta se pval eaf
#
# 用法: Rscript 04_coloc.R <trait1_region.tsv> <trait2_region.tsv> <标签>
#   例: Rscript 04_coloc.R data/processed/region_DED_BIN1.tsv \
#                          data/processed/region_AD_BIN1.tsv "chr2_BIN1"
#
# 说明: 本机未安装 coloc 包 (GitHub 不通), 这里直接实现 Giambartolomei 2013
#       的近似贝叶斯因子, 公式与 coloc::coloc.abf 一致。

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  files <- list.files("data/processed", pattern = "^region_", full.names = TRUE)
  if (length(files) == 0) stop("无区域文件, 请先用 _region.py 提取")
  cat("可用区域文件:\n"); print(files)
  stop("用法: Rscript 04_coloc.R <trait1.tsv> <trait2.tsv> <标签>")
}
f1 <- args[1]; f2 <- args[2]; label <- args[3]

# ---- 先验参数 (coloc 默认) ----
P1 <- 1e-4      # SNP 与性状1 关联的先验概率
P2 <- 1e-4      # SNP 与性状2 关联的先验概率
P12 <- 5e-6     # 同时与两者关联的先验概率
W <- 0.4^2      # 效应量先验方差 (log OR 尺度)

lsum <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

d1 <- fread(f1)[, .(SNP, effect_allele, other_allele, beta, se)]
d2 <- fread(f2)[, .(SNP, effect_allele, other_allele, beta, se)]
m <- merge(d1, d2, by = "SNP", suffixes = c(".1", ".2"))
if (nrow(m) == 0) stop("两个区域文件无共同 SNP")

ea1 <- toupper(m$effect_allele.1); oa1 <- toupper(m$other_allele.1)
ea2 <- toupper(m$effect_allele.2); oa2 <- toupper(m$other_allele.2)
same    <- ea1 == ea2 & oa1 == oa2
swapped <- !same & ea1 == oa2 & oa1 == ea2
comp <- function(a) { mp <- c(A = "T", T = "A", G = "C", C = "G")
                      ifelse(a %in% names(mp), mp[a], a) }
keep <- same | swapped
dropped <- sum(!keep)
fl <- swapped[keep]
m <- m[keep, ]
m$beta.2[fl] <- -m$beta.2[fl]
cat(sprintf("共同 SNP: %d (剔除方向不可判定 %d 个)\n", nrow(m), dropped))
if (nrow(m) < 10) warning("SNP 数偏少 (<10), 共定位结果不稳定")

V1 <- m$se.1^2; V2 <- m$se.2^2
z1 <- m$beta.1 / m$se.1; z2 <- m$beta.2 / m$se.2
r1 <- W / (W + V1); r2 <- W / (W + V2)
# 近似贝叶斯因子 (Wakefield 2009)
l1 <- 0.5 * log(1 - r1) + 0.5 * z1^2 * r1
l2 <- 0.5 * log(1 - r2) + 0.5 * z2^2 * r2

n <- nrow(m)
lH0 <- 0
lH1 <- log(P1) + lsum(l1) - log(n)
lH2 <- log(P2) + lsum(l2) - log(n)
lH3 <- log(P1) + log(P2) + lsum(l1 + l2) - log(n)
lH4 <- log(P12) + lsum(l1 + l2) - log(n)
den <- lsum(c(lH0, lH1, lH2, lH3, lH4))
pp <- exp(c(lH0, lH1, lH2, lH3, lH4) - den)
names(pp) <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")

cat("\n===== 共定位: ", label, " =====\n", sep = "")
cat(sprintf("SNP 数 %d | p1=%g p2=%g p12=%g W=%g\n", n, P1, P2, P12, W))
print(round(pp, 4))
cat(sprintf("\n结论: PP.H4 = %.3f -> %s\n", pp["PP.H4"],
            ifelse(pp["PP.H4"] > 0.8, "强共定位证据 (共享因果变异)",
            ifelse(pp["PP.H4"] > 0.5, "中等共定位证据", "无共定位证据"))))

top <- m[which.max(l1 + l2), ]
cat(sprintf("最可能共享因果变异: %s (性状1 z=%.2f, 性状2 z=%.2f)\n",
            top$SNP, top$beta.1 / top$se.1, top$beta.2 / top$se.2))

dir.create("results", showWarnings = FALSE)
out <- data.table(locus = label, nsnp = n, as.list(round(pp, 5)))
fn <- "results/coloc_results.csv"
if (file.exists(fn)) {
  fwrite(rbindlist(list(fread(fn), out), fill = TRUE), fn)
} else fwrite(out, fn)
cat("已追加到 results/coloc_results.csv\n")
