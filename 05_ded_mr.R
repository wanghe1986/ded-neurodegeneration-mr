# 方案 B · 干眼 x 神经退行 双向两样本 MR
#
# 暴露: 干眼 (Dry eyes, PheCode 375.1) —— MVP, GCST90475895 (EUR, 84,262 / 314,374)
# 结局: AD  Bellenguez 2022 (GCST90027158, 85,934 / 401,577)
#       PD  Nalls 2019      (GCST009325,   33,674 / 449,056)
#       MS  IMSGC 2019      (GCST005531,   47,429 / 68,374)
#
# 样本重叠: MVP(美国退伍军人) 与三个欧洲神经退行 GWAS 均无重叠 -> 两样本假设成立
#
# 输入: data/processed/iv_DED.tsv, out_<OUT>_for_DED.tsv     (正向)
#       data/processed/iv_<OUT>.tsv, out_DED_for_<OUT>.tsv   (反向)
# 输出: results/bidirectional_mr.csv

library(data.table)
source("mr_core.R")

DIR <- "data/processed"
OUTCOMES <- c(AD = "AD", PD = "PD", MS = "MS")

# MS 的病例/对照数来自文献 (GWAS Catalog harmonised 文件不含该列)
N_OUT <- list(AD = c(85934, 401577), PD = c(33674, 449056), MS = c(47429, 68374))
N_DED <- c(84262, 314374)

run_mr <- function(exp_file, out_file, label) {
  if (!file.exists(exp_file) || !file.exists(out_file)) {
    warning("跳过 ", label, " (缺文件)"); return(NULL)
  }
  e <- fread(exp_file); o <- fread(out_file)
  eh <- e[, .(SNP, effect_allele, other_allele, beta, se, pval)]
  oh <- o[, .(SNP, effect_allele, other_allele, beta, se, pval)]
  h  <- harmonize_pair(as.data.frame(eh), as.data.frame(oh))
  if (is.null(h) || nrow(h) == 0) { warning("跳过 ", label, " (harmonize 后为空)"); return(NULL) }
  r <- mr_estimate(h)
  r$analysis <- label
  r$nsnp_used <- nrow(h)
  r$nsnp_dropped_ambiguous <- nrow(e) - nrow(h)
  r
}

rows <- list()

# ---------- 正向: 干眼 -> 神经退行 ----------
for (nm in names(OUTCOMES)) {
  r <- run_mr(file.path(DIR, "iv_DED.tsv"),
              file.path(DIR, sprintf("out_%s_for_DED.tsv", nm)),
              sprintf("DED -> %s", nm))
  if (!is.null(r)) rows[[length(rows) + 1]] <- r
}

# ---------- 反向: 神经退行 -> 干眼 ----------
for (nm in names(OUTCOMES)) {
  r <- run_mr(file.path(DIR, sprintf("iv_%s.tsv", nm)),
              file.path(DIR, sprintf("out_DED_for_%s.tsv", nm)),
              sprintf("%s -> DED", nm))
  if (!is.null(r)) rows[[length(rows) + 1]] <- r
}

if (length(rows) == 0) stop("无任何分析完成")

res <- rbindlist(rows)
res[, direction := fifelse(grepl("^DED ->", analysis), "forward", "reverse")]
res[, outcome := fifelse(direction == "forward",
                         sub("DED -> ", "", analysis), sub(" -> DED", "", analysis))]
# 效应量换算: 暴露 log OR 每增加 1 个单位 -> 结局 OR
res[, OR := exp(b)]
res[, OR_lci := exp(b - 1.96 * se)]
res[, OR_uci := exp(b + 1.96 * se)]
# 每种方法单独做 BH 校正
res[, FDR := p.adjust(pval, method = "BH"), by = .(method)]

dir.create("results", showWarnings = FALSE)
fwrite(res, "results/bidirectional_mr.csv")

message("\n=========================================")
message(" 干眼 x 神经退行 双向 MR 结果")
message("=========================================")

prim <- res[method == "ivw_fe"]
if (nrow(prim) > 0) {
  message("\n----- 主估计量 IVW (固定效应) -----")
  print(prim[order(pval), .(analysis, nsnp_used, b, se, OR, OR_lci, OR_uci, pval, FDR, Q_pval)])
}

message("\n----- 敏感性分析: MR-Egger / 加权中位数 -----")
sens <- res[method %in% c("mr_egger", "weighted_median")]
if (nrow(sens) > 0) {
  print(sens[order(analysis, method),
        .(analysis, method, nsnp_used, b, se, OR, pval, FDR, egger_intercept, egger_pval)])
}

message("\n----- 异质性 (Cochran Q p<0.05 建议改用随机效应) -----")
q <- res[method == "ivw_fe" & !is.na(Q_pval)]
if (nrow(q) > 0) print(q[, .(analysis, nsnp_used, Q, Q_pval)])
ivwre <- res[method == "ivw_re" & !is.na(pval)]
if (nrow(ivwre) > 0) {
  message("\n----- IVW 随机效应 (异质性存在时的稳健估计) -----")
  print(ivwre[, .(analysis, b, se, OR, pval, FDR)])
}

message("\n判定要点:")
message("  1) 主估计 IVW p<0.05 且 FDR<0.05 -> 有因果证据")
message("  2) MR-Egger 与加权中位数方向一致 -> 结果稳健")
message("  3) Egger 截距 p>0.05 -> 无明显定向多效性")
message("  4) Steiger 方向: 需确认 SNP 对暴露的解释度大于对结局")
message("\n已写出 results/bidirectional_mr.csv")
