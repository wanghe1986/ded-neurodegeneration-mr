# 通用配对 MR: 给定暴露工具变量文件与结局提取文件, 输出全部估计量
# 用法: Rscript _mr_pair.R <暴露iv.tsv> <结局out.tsv> <标签> [ncx nkx ncy nky]
#
# !! 样本量陷阱 !!  GWAS Catalog harmonised 文件的 n_cases/n_controls 常是
#   「逐 SNP 分批样本量」(PD 文件首行仅 1169/968), 不是研究总样本量。
#   直接用会把 Steiger 的 r2 放大数十倍, 造成"方向错误"的假警报。
#   务必用论文报告的总样本量以第 4-7 个参数覆盖。
library(data.table)
source("mr_core.R")
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3) stop("用法: Rscript _mr_pair.R <iv.tsv> <out.tsv> <标签>")

e <- fread(a[1]); o <- fread(a[2])
eh <- e[, .(SNP, effect_allele, other_allele, beta, se, pval)]
oh <- o[, .(SNP, effect_allele, other_allele, beta, se, pval)]
h <- harmonize_pair(as.data.frame(eh), as.data.frame(oh))
if (is.null(h) || nrow(h) == 0) stop("harmonize 后为空")

r <- mr_estimate(h)
setDT(r)
r[, OR := exp(b)]; r[, LCI := exp(b - 1.96 * se)]; r[, UCI := exp(b + 1.96 * se)]
cat("\n===== ", a[3], " =====\n", sep = "")
cat("工具变量: 暴露 ", nrow(e), " 个, 结局命中 ", nrow(o),
    " 个, harmonize 后可用 ", nrow(h), " 个\n", sep = "")
print(r[, .(method, nsnp, b, se, OR, LCI, UCI, pval, Q_pval, egger_intercept, egger_pval)])

# 暴露与结局都是二分类时, 报告 Steiger 方向
# 论文报告的总样本量 (可用命令行第 4-7 参数覆盖)
ncase_x <- if (length(a) >= 5) as.numeric(a[4]) else NA
nctrl_x <- if (length(a) >= 6) as.numeric(a[5]) else NA
ncase_y <- if (length(a) >= 7) as.numeric(a[6]) else NA
nctrl_y <- if (length(a) >= 8) as.numeric(a[7]) else NA
if (!is.na(ncase_x) && !is.na(nctrl_x)) {
  e$n_cases <- ncase_x; e$n_controls <- nctrl_x
}
if (!is.na(ncase_y) && !is.na(nctrl_y)) {
  o$n_cases <- ncase_y; o$n_controls <- nctrl_y
}

if (all(c("n_cases", "n_controls") %in% names(e)) &&
    all(c("n_cases", "n_controls") %in% names(o))) {
  nx <- 4 * e$n_cases[1] * e$n_controls[1] / (e$n_cases[1] + e$n_controls[1])
  ny <- 4 * o$n_cases[1] * o$n_controls[1] / (o$n_cases[1] + o$n_controls[1])
  cat(sprintf("样本量(用于Steiger): 暴露 %d/%d (Neff=%.0f), 结局 %d/%d (Neff=%.0f)\n",
              e$n_cases[1], e$n_controls[1], nx,
              o$n_cases[1], o$n_controls[1], ny))
  r2x <- h$beta.exp^2 / (h$beta.exp^2 + h$se.exp^2 * nx)
  r2y <- h$beta.out^2 / (h$beta.out^2 + h$se.out^2 * ny)
  cat(sprintf("Steiger: 方向正确比例 %.2f (r2暴露=%.3g, r2结局=%.3g)\n",
              mean(r2x > r2y, na.rm = TRUE), mean(r2x, na.rm = TRUE), mean(r2y, na.rm = TRUE)))
}
dir.create("results", showWarnings = FALSE)
fn <- "results/mr_pairs.csv"
row <- data.table(pair = a[3], nsnp_iv = nrow(e), nsnp_used = nrow(h), r)
if (file.exists(fn)) fwrite(rbindlist(list(fread(fn), row), fill = TRUE), fn) else fwrite(row, fn)
cat("已追加 results/mr_pairs.csv\n")
