# 方案 B · 敏感性分析
# 1) 留一法 (leave-one-out)
# 2) 区域剔除: APOE (AD 的单基因驱动风险) / HLA (MS 的多效性重灾区)
# 3) 仅保留强工具 (F > 100)
# 4) Steiger 定向检验 (确认 SNP 先影响暴露而非结局)
#
# 依赖 05_ded_mr.R 已生成的 data/processed/iv_*.tsv 与 out_*_for_*.tsv

library(data.table)
source("mr_core.R")

DIR <- "data/processed"
OUTS <- c("AD", "PD", "MS")

# 有效样本量 (二分类 GWAS 的 Steiger 需要)
NEFF <- function(nc, nk) 4 * nc * nk / (nc + nk)

# ---- 区域定义 (GRCh38) ----
REGIONS <- list(
  APOE = data.frame(chr = "19", start = 43.9e6, end = 45.9e6),
  HLA   = data.frame(chr = "6",  start = 25.0e6, end = 36.0e6)
)

load_pair <- function(exp_name, out_name) {
  ef <- file.path(DIR, sprintf("iv_%s.tsv", exp_name))
  of <- file.path(DIR, sprintf("out_%s_for_%s.tsv", out_name, exp_name))
  if (!file.exists(ef) || !file.exists(of)) return(NULL)
  e <- fread(ef); o <- fread(of)
  # 剔除区域会用到 chr/pos, 若暴露文件没有就从结局补
  eh <- e[, .(SNP, effect_allele, other_allele, beta, se, pval)]
  oh <- o[, .(SNP, effect_allele, other_allele, beta, se, pval)]
  h <- harmonize_pair(as.data.frame(eh), as.data.frame(oh))
  if (is.null(h) || nrow(h) == 0) return(NULL)
  keep <- e[SNP %in% h$SNP, .(SNP, chr, pos, F, eaf, beta, se, n_cases, n_controls)]
  h <- merge(h, keep[, .(SNP, chr, pos, F, eaf, n_cases, n_controls)], by = "SNP")
  list(h = h, exp_iv = e)
}

ivw <- function(d) {
  bx <- d$beta.exp; by <- d$beta.out; sy <- d$se.out
  w <- 1 / sy^2
  b <- sum(w * bx * by) / sum(w * bx^2)
  se <- sqrt(1 / sum(w * bx^2))
  data.frame(b = b, se = se, pval = 2 * pnorm(-abs(b / se)), nsnp = nrow(d))
}

steiger <- function(d, n_exp, n_out) {
  # r2 = beta^2 / (beta^2 + se^2 * Neff); 方向正确要求 r2_x > r2_y
  r2x <- d$beta.exp^2 / (d$beta.exp^2 + d$se.exp^2 * n_exp)
  r2y <- d$beta.out^2 / (d$beta.out^2 + d$se.out^2 * n_out)
  data.frame(frac_correct = mean(r2x > r2y, na.rm = TRUE),
             mean_r2_exp = mean(r2x, na.rm = TRUE),
             mean_r2_out = mean(r2y, na.rm = TRUE))
}

analyses <- list(
  list(name = "DED -> AD", exp = "DED", out = "AD"),
  list(name = "DED -> PD", exp = "DED", out = "PD"),
  list(name = "DED -> MS", exp = "DED", out = "MS"),
  list(name = "AD -> DED", exp = "AD",  out = "DED"),
  list(name = "PD -> DED", exp = "PD",  out = "DED"),
  list(name = "MS -> DED", exp = "MS",  out = "DED")
)

rows <- list(); loo_rows <- list()
for (a in analyses) {
  p <- load_pair(a$exp, a$out)
  if (is.null(p)) next
  d <- p$h
  n_exp <- if (all(c("n_cases", "n_controls") %in% names(d)) && !any(is.na(d$n_cases)))
             NEFF(d$n_cases[1], d$n_controls[1]) else NA
  n_out <- if ("n_cases" %in% names(d)) NA else NA
  message("\n===== ", a$name, " (", nrow(d), " 个工具变量) =====")

  base <- ivw(d); message(sprintf("  主 IVW: b=%+.4f se=%.4f p=%.3g", base$b, base$se, base$pval))

  # --- Steiger ---
  if (!is.na(n_exp)) {
    # 结局有效样本量: 从结局文件取
    of <- file.path(DIR, sprintf("out_%s_for_%s.tsv", a$out, a$exp))
    if (file.exists(of)) {
      oo <- fread(of)
      if (all(c("n_cases", "n_controls") %in% names(oo)))
        n_out <- NEFF(oo$n_cases[1], oo$n_controls[1])
    }
  }
  if (!is.na(n_exp) && !is.na(n_out)) {
    s <- steiger(d, n_exp, n_out)
    message(sprintf("  Steiger: 方向正确比例 %.2f  (r2暴露=%.3g, r2结局=%.3g)",
                     s$frac_correct, s$mean_r2_exp, s$mean_r2_out))
    rows[[length(rows) + 1]] <- data.frame(analysis = a$name, check = "steiger",
                                           nsnp = nrow(d), b = NA, se = NA, pval = NA,
                                           note = sprintf("correct=%.2f", s$frac_correct))
  }

  # --- 留一法 ---
  for (i in seq_len(nrow(d))) {
    r <- ivw(d[-i, ]); r$analysis <- a$name; r$dropped <- d$SNP[i]
    loo_rows[[length(loo_rows) + 1]] <- r
  }
  lo <- rbindlist(loo_rows[grepl(a$name, sapply(loo_rows, function(x) x$analysis))])
  if (nrow(lo) > 0)
    message(sprintf("  留一法: b 范围 [%+.4f, %+.4f], p 范围 [%.3g, %.3g]",
                    min(lo$b), max(lo$b), min(lo$pval), max(lo$pval)))

  # --- 区域剔除 ---
  for (rg in names(REGIONS)) {
    r0 <- REGIONS[[rg]]
    sub <- d[!(as.character(chr) == r0$chr & pos >= r0$start & pos <= r0$end), ]
    if (nrow(sub) == 0 || nrow(sub) == nrow(d)) next
    r <- ivw(sub)
    message(sprintf("  剔除 %-5s: n=%d  b=%+.4f  p=%.3g", rg, nrow(sub), r$b, r$pval))
    rows[[length(rows) + 1]] <- data.frame(analysis = a$name,
                                           check = paste0("drop_", rg),
                                           nsnp = nrow(sub), b = r$b, se = r$se,
                                           pval = r$pval, note = "")
  }

  # --- 强工具子集 ---
  if ("F" %in% names(d)) {
    sub <- d[!is.na(F) & F > 100, ]
    if (nrow(sub) >= 2) {
      r <- ivw(sub)
      message(sprintf("  仅 F>100  : n=%d  b=%+.4f  p=%.3g", nrow(sub), r$b, r$pval))
      rows[[length(rows) + 1]] <- data.frame(analysis = a$name, check = "F_gt_100",
                                            nsnp = nrow(sub), b = r$b, se = r$se,
                                            pval = r$pval, note = "")
    }
  }
}

if (length(rows)) {
  res <- rbindlist(rows)
  res[, OR := exp(b)]
  dir.create("results", showWarnings = FALSE)
  fwrite(res, "results/sensitivity.csv")
  message("\n敏感性分析已写出 results/sensitivity.csv")
}
if (length(loo_rows)) {
  loo <- rbindlist(loo_rows)
  fwrite(loo, "results/leave_one_out.csv")
  message("留一法已写出 results/leave_one_out.csv")
}
