# 方案 B · 结果可视化 (base graphics, 不依赖 ggplot2)
#
# 输出 (600 dpi PNG, 符合 Eye and Vision 等期刊的图件规格):
#   results/fig_forest.png       双向 MR 森林图 (IVW)
#   results/fig_scatter_<分析>.png   SNP 效应散点图
#   results/fig_loo_<分析>.png       留一法
#
# 用法: Rscript 08_plots.R

library(data.table)
source("mr_core.R")
DIR <- "data/processed"
dir.create("results", showWarnings = FALSE)

# ---------- 载入主结果 ----------
res <- fread("results/bidirectional_mr.csv")
prim <- res[method == "ivw_fe"][order(pval)]

# ---------- 图1 森林图 ----------
png("results/fig_forest.png", width = 2400, height = 1800, res = 600)
op <- par(mar = c(5, 22, 4, 2))
labs <- prim$analysis
n <- nrow(prim)
y <- rev(seq_len(n))
xlim <- range(c(prim$OR_lci, prim$OR_uci, 1), na.rm = TRUE)
xlim <- c(max(0.5, min(xlim)), min(2.0, max(xlim)))
plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5), yaxt = "n", ylab = "",
     xlab = "OR per 1-unit increase in genetically predicted liability (log OR scale)",
     main = "Bidirectional MR: dry eye disease vs neurodegenerative diseases",
     log = "x")
abline(v = 1, lty = 2, col = "grey50")
grid(nx = NULL, ny = NA, col = "grey85")
cols <- ifelse(prim$pval < 0.05, "#c0392b", "#2c3e50")
segments(prim$OR_lci, y, prim$OR_uci, y, lwd = 2, col = cols)
points(prim$OR, y, pch = 18, cex = 1.6, col = cols)
axis(2, at = y, labels = labs, las = 1, cex.axis = 0.62)
for (i in seq_len(n)) {
  text(xlim[2], y[i], sprintf("%.3f (%.2f-%.2f), p=%.2g",
       prim$OR[i], prim$OR_lci[i], prim$OR_uci[i], prim$pval[i]),
       adj = c(1.05, 0.5), cex = 0.5, col = cols[i])
}
legend("bottomleft", c("p < 0.05", "p >= 0.05"), col = c("#c0392b", "#2c3e50"),
       pch = 18, bty = "n", cex = 0.7)
par(op); dev.off()
cat("图1 -> results/fig_forest.png\n")

# ---------- 图2 散点图 / 图3 留一法 ----------
analyses <- list(
  list(name = "DED -> AD", exp = "DED", out = "AD"),
  list(name = "DED -> PD", exp = "DED", out = "PD"),
  list(name = "DED -> MS", exp = "DED", out = "MS"),
  list(name = "AD -> DED", exp = "AD",  out = "DED"),
  list(name = "PD -> DED", exp = "PD",  out = "DED"),
  list(name = "MS -> DED", exp = "MS",  out = "DED")
)

for (a in analyses) {
  ef <- file.path(DIR, sprintf("iv_%s.tsv", a$exp))
  of <- file.path(DIR, sprintf("out_%s_for_%s.tsv", a$out, a$exp))
  if (!file.exists(ef) || !file.exists(of)) next
  e <- fread(ef); o <- fread(of)
  h <- harmonize_pair(
    as.data.frame(e[, .(SNP, effect_allele, other_allele, beta, se, pval)]),
    as.data.frame(o[, .(SNP, effect_allele, other_allele, beta, se, pval)]))
  if (is.null(h) || nrow(h) == 0) next
  tag <- gsub("[^A-Za-z0-9]", "", a$name)

  # 散点图
  png(sprintf("results/fig_scatter_%s.png", tag), width = 1800, height = 1600, res = 600)
  plot(h$beta.exp, h$beta.out, pch = 19, col = "#2c3e50",
       xlab = sprintf("SNP effect on %s (log OR)", a$exp),
       ylab = sprintf("SNP effect on %s (log OR)", a$out),
       main = a$name)
  abline(h = 0, v = 0, col = "grey80")
  r <- mr_estimate(h)
  for (mth in c("ivw_fe", "mr_egger", "weighted_median")) {
    rr <- r[method == mth]
    if (nrow(rr) == 0) next
    abline(a = 0, b = rr$b, lwd = 2,
           col = c(ivw_fe = "#c0392b", mr_egger = "#2980b9",
                   weighted_median = "#27ae60")[mth])
  }
  legend("topleft", c("IVW", "MR-Egger", "Weighted median"),
         col = c("#c0392b", "#2980b9", "#27ae60"), lwd = 2, bty = "n", cex = 0.8)
  dev.off()

  # 留一法
  loo <- rbindlist(lapply(seq_len(nrow(h)), function(i) {
    rr <- mr_estimate(h[-i, ]); rr <- rr[method == "ivw_fe"]
    rr$dropped <- h$SNP[i]; rr
  }))
  png(sprintf("results/fig_loo_%s.png", tag), width = 1800, height = 1600, res = 600)
  plot(seq_len(nrow(loo)), loo$b, pch = 19, col = "#2c3e50", xaxt = "n",
       xlab = "Instrument excluded", ylab = "IVW estimate (beta)",
       main = paste("Leave-one-out:", a$name))
  arrows(seq_len(nrow(loo)), loo$b - 1.96 * loo$se,
         seq_len(nrow(loo)), loo$b + 1.96 * loo$se, code = 3, angle = 90,
         length = 0.03, col = "#7f8c8d")
  abline(h = 0, col = "grey70", lty = 2)
  dev.off()
  cat("图 -> fig_scatter_", tag, ".png / fig_loo_", tag, ".png\n", sep = "")
}
