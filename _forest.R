# 双向 MR 森林图 (读 results/mr_pairs.csv 的 IVW 结果)
# 规格: 170 mm 宽 (期刊全页宽), 600 dpi
library(data.table)
d <- fread("results/mr_pairs.csv")
d <- d[method == "ivw_fe"]
if (nrow(d) == 0) stop("mr_pairs.csv 中没有 ivw_fe 结果")
setorder(d, pval)
# 简化标签
d[, lab := pair]
# 用不含中文的关键字匹配, 避免 Windows 下 R 的编码问题
d[grepl("^DED", lab) & grepl("AD$", lab),                    lab := "Dry eye -> AD (97 IV)"]
d[grepl("^DED", lab) & grepl("MS", lab),                     lab := "Dry eye -> MS (70 IV, FinnGen)"]
d[grepl("^DED", lab) & grepl("PD", lab) & grepl("130", lab), lab := "Dry eye -> PD (84 IV, extended)"]
d[grepl("^DED", lab) & grepl("PD", lab) & grepl("9 IV", lab), lab := "Dry eye -> PD (6 IV, strict)"]
d[grepl("^AD",  lab), lab := "AD -> dry eye (40 IV)"]
d[grepl("^PD",  lab), lab := "PD -> dry eye (24 IV)"]
d[grepl("^MS",  lab), lab := "MS -> dry eye (14 IV, IMSGC)"]

dir.create("results", showWarnings = FALSE)
png("results/fig_forest.png", width = 4016, height = 2600, res = 600)
op <- par(mar = c(5, 21, 4, 12), mgp = c(3.5, 1, 0))

n <- nrow(d); y <- rev(seq_len(n))
xlim <- range(c(d$LCI, d$UCI, 1), finite = TRUE)
xlim <- c(min(0.30, xlim[1] * 0.9), max(1.40, xlim[2] * 1.1))

plot(NA, xlim = xlim, ylim = c(0.4, n + 0.6), yaxt = "n", ylab = "", log = "x",
     xlab = "Odds ratio (95% CI) per 1-unit increase in genetically predicted liability",
     main = "Bidirectional Mendelian randomization: dry eye disease and neurodegenerative diseases",
     cex.main = 1.5, cex.lab = 1.15, cex.axis = 1.05)
abline(v = 1, lty = 2, lwd = 1.5, col = "grey45")
abline(v = c(0.5, 0.8, 1.25, 2), col = "grey88", lty = 3)

colv <- ifelse(d$pval < 0.05, "#c0392b", "#2c3e50")
cexv <- ifelse(d$pval < 0.05, 2.4, 2.0)
segments(d$LCI, y, d$UCI, y, lwd = 3, col = colv)
points(d$OR, y, pch = 18, cex = cexv, col = colv)

axis(2, at = y, labels = d$lab, las = 1, cex.axis = 1.05, line = -1)
mtext(sprintf("%.3f (%.2f-%.2f)", d$OR, d$LCI, d$UCI), side = 4, at = y,
      las = 1, cex = 1.0, line = 1, col = colv, adj = 0)
mtext(sprintf("p=%.3g", d$pval), side = 4, at = y, las = 1, cex = 1.0,
      line = 6, col = colv, adj = 0)
mtext(sprintf("n=%d", d$nsnp_used), side = 2, at = y, las = 1, cex = 0.9,
      line = 17, col = "grey35", adj = 1)

legend("bottomleft", c("p < 0.05", "p >= 0.05"), col = c("#c0392b", "#2c3e50"),
       pch = 18, pt.cex = 2, bty = "n", cex = 1.15)
par(op)
dev.off()
cat("图已写出 results/fig_forest.png\n")
print(d[, .(pair, nsnp_used, OR = round(OR, 3), LCI = round(LCI, 3),
            UCI = round(UCI, 3), pval = round(pval, 4))])
