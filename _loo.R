# 留一法: 逐个剔除工具变量, 检查 MR 估计是否被单个 SNP 驱动
# 用法: Rscript _loo.R <iv.tsv> <out.tsv> <标签>
library(data.table)
source("mr_core.R")
a <- commandArgs(trailingOnly = TRUE)
e <- fread(a[1]); o <- fread(a[2])
h <- harmonize_pair(
  as.data.frame(e[, .(SNP, effect_allele, other_allele, beta, se, pval)]),
  as.data.frame(o[, .(SNP, effect_allele, other_allele, beta, se, pval)]))
cat("\n===== 留一法: ", a[3], " (", nrow(h), " 个 IV) =====\n", sep = "")

base <- mr_estimate(h); base <- base[base$method == "ivw_fe", ]
cat(sprintf("全样本 IVW: b=%+.4f  OR=%.3f  p=%.4g\n\n", base$b, exp(base$b), base$pval))
cat(sprintf("%-16s %-10s %-10s %-10s %-10s\n", "剔除的SNP", "b", "OR", "p", "变化"))
rows <- list()
for (i in seq_len(nrow(h))) {
  r <- mr_estimate(h[-i, ]); r <- r[r$method == "ivw_fe", ]
  cat(sprintf("%-16s %+-10.4f %-10.3f %-10.4g %+-9.1f%%\n",
              h$SNP[i], r$b, exp(r$b), r$pval, 100 * (r$b - base$b) / abs(base$b)))
  rows[[i]] <- data.table(dropped = h$SNP[i], b = r$b, se = r$se,
                          OR = exp(r$b), pval = r$pval)
}
d <- rbindlist(rows)
cat(sprintf("\n剔除后 OR 范围: %.3f - %.3f\n", min(d$OR), max(d$OR)))
cat(sprintf("剔除后 p 范围 : %.4g - %.4g\n", min(d$pval), max(d$pval)))
cat(sprintf("方向是否全部一致: %s\n",
            ifelse(all(sign(d$b) == sign(base$b)), "是 (稳健)", "否 (存在方向翻转, 结果不稳健)")))
if (any(d$pval > 0.05) && base$pval < 0.05)
  cat("警示: 存在剔除某个 SNP 后 p 值跨过 0.05 的情况, 结论对个别 SNP 敏感\n")
