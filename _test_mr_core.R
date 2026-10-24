source("mr_core.R")
set.seed(20261022)

n <- 20
true_b <- 0.5
bx <- rnorm(n, 0.3, 0.08);  sx <- abs(rnorm(n, 0.02, 0.005))
by <- true_b * bx + rnorm(n, 0, 0.01); sy <- abs(rnorm(n, 0.02, 0.005))
dat <- data.frame(SNP = paste0("rs", 1:n), beta.exp = bx, se.exp = sx,
                  beta.out = by, se.out = sy)

cat("=== 多 SNP 估计 (真实因果效应 = 0.5) ===\n")
print(mr_estimate(dat), digits = 4)

cat("\n=== 单 SNP: Wald ratio ===\n")
print(mr_estimate(dat[1, ]), digits = 4)

cat("\n=== harmonization 测试 ===\n")
e <- data.frame(SNP = c("rs1","rs2","rs3"),
                effect_allele = c("A","A","A"), other_allele = c("G","G","G"),
                beta = c(0.1,0.2,0.3), se = 0.01)
o <- data.frame(SNP = c("rs1","rs2","rs3"),
                effect_allele = c("A","G","T"), other_allele = c("G","A","C"),
                beta = c(0.05,0.10,0.15), se = 0.01)
h <- harmonize_pair(e, o)
cat("保留", nrow(h), "个 SNP;  状态:", paste(h$harmonize_status, collapse = ", "), "\n")
cat("rs2 属 swapped, 结局 beta 应为 -0.10 :", h$beta.out[h$SNP == "rs2"], "\n")
cat("rs3 属回文/不匹配, 应被剔除\n")

cat("\n=== 无偏性检查: 重复 200 次看 IVW 均值 ===\n")
est <- replicate(200, {
  bx <- rnorm(n, 0.3, 0.08); sx <- abs(rnorm(n, 0.02, 0.005))
  by <- 0.5 * bx + rnorm(n, 0, 0.01); sy <- abs(rnorm(n, 0.02, 0.005))
  d <- data.frame(SNP = 1:n, beta.exp = bx, se.exp = sx, beta.out = by, se.out = sy)
  mr_estimate(d)$b[mr_estimate(d)$method == "ivw_fe"]
})
cat("IVW 均值 =", round(mean(est), 4), "  (真实 0.5)\n")
cat("IVW 标准差 =", round(sd(est), 4), "\n")
