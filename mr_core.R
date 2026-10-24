# A1 方案 · MR 核心算法自实现
#
# 目的: 不依赖 TwoSampleMR (本机 GitHub API 不通, 安装失败两次)
# MR 的数学是公开的, 这里直接实现标准估计量:
#   Wald ratio / IVW (固定效应+随机效应) / MR-Egger / 加权中位数
#   外加 Cochran Q 异质性、Egger 截距多效性、Steiger 定向检验
#
# 参考公式与 TwoSampleMR / MendelianRandomization 包一致

# ---------------- 等位基因 harmonization ----------------
harmonize_pair <- function(exp_dat, out_dat) {
  # exp_dat/out_dat: data.frame, 至少含 SNP, effect_allele, other_allele, beta, se, pval
  d <- merge(exp_dat, out_dat, by = "SNP", suffixes = c(".exp", ".out"))
  if (nrow(d) == 0) return(d)

  ea1 <- toupper(d$effect_allele.exp); oa1 <- toupper(d$other_allele.exp)
  ea2 <- toupper(d$effect_allele.out); oa2 <- toupper(d$other_allele.out)

  comp <- function(a) {
    a <- toupper(a)
    map <- c(A = "T", T = "A", G = "C", C = "G")
    ifelse(a %in% names(map), map[a], a)
  }

  status <- character(nrow(d))
  flip   <- logical(nrow(d))

  # 情况1: 效应等位基因一致
  same <- ea1 == ea2 & oa1 == oa2
  status[same] <- "direct"

  # 情况2: 完全互换 (效应等位 = 结局的非效应等位)
  swapped <- !same & ea1 == oa2 & oa1 == ea2
  status[swapped] <- "swapped"
  flip[swapped] <- TRUE

  # 情况3: 回文 (A/T 或 G/C), 方向不可判定
  pal <- comp(ea1) == oa1
  ambiguous <- !same & !swapped & pal
  status[ambiguous] <- "ambiguous"

  # 情况4: 其他不匹配
  status[status == ""] <- "mismatch"

  d$harmonize_status <- status

  # 对 swapped 翻转结局效应 (等价于翻转暴露效应)
  d$beta.out[flip] <- -d$beta.out[flip]

  # 只保留可判定方向的
  d <- d[d$harmonize_status %in% c("direct", "swapped"), ]
  d
}

# ---------------- 估计量 ----------------
# dat: harmonize_pair() 的输出, 含 beta.exp, se.exp, beta.out, se.out
mr_estimate <- function(dat, method = c("ivw_fe", "ivw_re", "egger", "wmedian", "wald")) {

  bx <- dat$beta.exp; by <- dat$beta.out
  sx <- dat$se.exp;   sy <- dat$se.out
  n  <- nrow(dat)

  # 单 SNP 时用 Wald ratio
  if (n == 1) {
    b <- by / bx
    se <- sy / abs(bx)
    return(data.frame(method = "wald_ratio", b = b, se = se,
                      pval = 2 * pnorm(-abs(b / se)), nsnp = 1,
                      Q = NA, Q_pval = NA, egger_intercept = NA, egger_pval = NA))
  }

  w <- 1 / (sy^2)          # IVW 权重 (标准做法: 仅用结局 SE)

  out <- list()

  # ---- IVW (固定效应) ----
  b_ivw <- sum(w * bx * by) / sum(w * bx^2)
  se_ivw <- sqrt(1 / sum(w * bx^2))
  Q <- sum(w * (by - b_ivw * bx)^2)
  Q_df <- n - 1
  Q_p <- pchisq(Q, Q_df, lower.tail = FALSE)
  out$ivw_fe <- data.frame(method = "ivw_fe", b = b_ivw, se = se_ivw,
                           pval = 2 * pnorm(-abs(b_ivw / se_ivw)),
                           nsnp = n, Q = Q, Q_pval = Q_p,
                           egger_intercept = NA, egger_pval = NA)

  # ---- IVW (随机效应): SE 乘以 sqrt(max(1,Q/df)) ----
  se_ivw_re <- se_ivw * sqrt(max(1, Q / Q_df))
  out$ivw_re <- data.frame(method = "ivw_re", b = b_ivw, se = se_ivw_re,
                           pval = 2 * pnorm(-abs(b_ivw / se_ivw_re)),
                           nsnp = n, Q = Q, Q_pval = Q_p,
                           egger_intercept = NA, egger_pval = NA)

  # ---- MR-Egger: by ~ bx 的加权线性回归, 权重 1/sy^2 ----
  # 斜率 = 因果估计; 截距 != 0 提示定向多效性
  fit <- lm(by ~ bx, weights = w)
  cf <- summary(fit)$coefficients
  b_egg <- cf["bx", "Estimate"]; se_egg <- cf["bx", "Std. Error"]
  b_int <- cf["(Intercept)", "Estimate"]; se_int <- cf["(Intercept)", "Std. Error"]
  out$egger <- data.frame(method = "mr_egger", b = b_egg, se = se_egg,
                          pval = cf["bx", "Pr(>|t|)"], nsnp = n,
                          Q = Q, Q_pval = Q_p,
                          egger_intercept = b_int,
                          egger_pval = cf["(Intercept)", "Pr(>|t|)"])

  # ---- 加权中位数 ----
  # 按 |beta| 加权的中位数 (Bowden 2016 的简化实现)
  ratio <- by / bx
  wm_w <- w * bx^2                        # 权重
  ord <- order(ratio)
  ratio_s <- ratio[ord]; wm_w_s <- wm_w[ord]
  cw <- cumsum(wm_w_s) / sum(wm_w_s)
  idx <- which(cw >= 0.5)[1]
  b_wm <- ratio_s[idx]
  # bootstrap 估计 SE
  B <- 1000
  bs <- replicate(B, {
    s <- sample(seq_len(n), n, replace = TRUE)
    rr <- (by / bx)[s]; ww <- wm_w[s]
    o <- order(rr); cc <- cumsum(ww[o]) / sum(ww[o])
    rr[o][which(cc >= 0.5)[1]]
  })
  se_wm <- sd(bs)
  out$wmedian <- data.frame(method = "weighted_median", b = b_wm, se = se_wm,
                            pval = 2 * pnorm(-abs(b_wm / se_wm)),
                            nsnp = n, Q = Q, Q_pval = Q_p,
                            egger_intercept = NA, egger_pval = NA)

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

# ---------------- Steiger 定向检验 ----------------
# 比较暴露与结局各自被 SNP 解释的方差, 判断因果方向是否 X -> Y
steiger_test <- function(dat, n_exp, n_out) {
  r2_x <- 2 * (dat$beta.exp^2) * dat$maf * (1 - dat$maf) /
          (2 * (dat$beta.exp^2) * dat$maf * (1 - dat$maf) +
           2 * n_exp * (dat$se.exp^2) * dat$maf * (1 - dat$maf))
  r2_y <- 2 * (dat$beta.out^2) * dat$maf * (1 - dat$maf) /
          (2 * (dat$beta.out^2) * dat$maf * (1 - dat$maf) +
           2 * n_out * (dat$se.out^2) * dat$maf * (1 - dat$maf))
  # 方向正确 (X->Y) 要求 r2_x > r2_y
  correct <- mean(r2_x > r2_y, na.rm = TRUE)
  data.frame(frac_correct_direction = correct,
             mean_r2_exposure = mean(r2_x, na.rm = TRUE),
             mean_r2_outcome  = mean(r2_y, na.rm = TRUE))
}
