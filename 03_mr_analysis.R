# A1 方案 · MR 主分析
#
# 使用 mr_core.R 的自实现估计量 (本机 TwoSampleMR 装不上, 已验证自实现无偏)
# 输入: data/processed/exposure_top_snp.csv  (每基因 top cis-eQTL, 单 IV)
#       data/processed/AD_target.tsv         (Bellenguez 2022 AD GWAS, 已按目标 SNP 过滤)
# 输出: results/mr_results.csv
#
# 说明: 未做 LD clumping (本机无 PLINK + 1000G 参考面板), 故主分析每基因只取
#       最显著 cis-eQTL 作单工具变量 (Wald ratio), 规避连锁偏倚。

# ================= CONFIG =================
EXPOSURE_FILE <- "data/processed/exposure_top_snp.csv"
OUTCOME_FILES <- c(
  AD = "data/processed/AD_target.tsv"
  # PD = "data/processed/PD_target.tsv",   # 待获取
  # MS = "data/processed/MS_target.tsv"    # 待获取
)
# 结局 GWAS 的列名 -> 标准列名映射
OUT_COLS <- c(SNP = "SNP", effect_allele = "effect_allele",
              other_allele = "other_allele", beta = "beta", se = "se", pval = "pval")
# ==========================================

library(data.table)
source("mr_core.R")
source("01_gene_sets.R")

if (!file.exists(EXPOSURE_FILE)) {
  stop("缺少 ", EXPOSURE_FILE, "。请先运行 02_extract_instruments.R")
}
exp <- fread(EXPOSURE_FILE)
message("暴露: ", nrow(exp), " 行, ", uniqueN(exp$gene), " 个基因")

# ---- 基因覆盖诊断 ----
all_genes <- unique(c(main_set, control_inflammation, positive_control))
have <- unique(exp$gene)
missing <- setdiff(all_genes, have)
message("\n===== 基因覆盖 =====")
message("设计基因数: ", length(all_genes), " | 有 cis-eQTL 工具变量: ", length(have),
        " | 缺失: ", length(missing))
if (length(missing) > 0) {
  message("\n缺失工具变量的基因 (按分层):")
  for (nm in c("main_set", "control_inflammation", "positive_control")) {
    g <- intersect(get(nm), missing)
    if (length(g) > 0) message("  ", nm, " (", length(get(nm)), "): ", paste(g, collapse = ", "))
  }
}

read_outcome <- function(path, name) {
  if (!file.exists(path)) { warning("缺失结局文件: ", path); return(NULL) }
  dt <- fread(path, showProgress = FALSE)
  src <- names(OUT_COLS)
  miss <- setdiff(src, names(dt))
  if (length(miss) > 0) {
    warning(name, " 缺少列: ", paste(miss, collapse = ", "), "  (请调整 OUT_COLS)")
    return(NULL)
  }
  dt <- dt[, ..src]
  setnames(dt, src, unname(OUT_COLS))
  dt[, outcome := name]
  dt
}

all_rows <- list()
for (nm in names(OUTCOME_FILES)) {
  out <- read_outcome(OUTCOME_FILES[[nm]], nm)
  if (is.null(out)) next
  message("\n--- 结局: ", nm, " (", nrow(out), " 个 SNP) ---")

  for (g in unique(exp$gene)) {
    e <- exp[gene == g]
    if (nrow(e) == 0) next
    o <- out[SNP %in% e$SNP]
    if (nrow(o) == 0) next

    eh <- e[, .(SNP, effect_allele, other_allele, beta, se, pval)]
    oh <- o[, .(SNP, effect_allele, other_allele, beta, se, pval)]
    h  <- harmonize_pair(as.data.frame(eh), as.data.frame(oh))
    if (is.null(h) || nrow(h) == 0) next

    r <- mr_estimate(h)
    if (is.null(r) || nrow(r) == 0) next
    r$gene <- g; r$outcome <- nm; r$nsnp_used <- nrow(h)
    r$top_snp <- paste(h$SNP, collapse = ";")
    r$harmonize_status <- paste(unique(h$harmonize_status), collapse = ";")
    all_rows[[length(all_rows) + 1]] <- r
  }
}

if (length(all_rows) == 0) stop("没有任何基因-结局组合完成分析")

res <- rbindlist(all_rows)
res[, set_role := fifelse(gene %in% main_set, "main",
                  fifelse(gene %in% control_inflammation, "control", "positive_control"))]
res[, FDR := p.adjust(pval, method = "BH"), by = .(outcome, set_role, method)]

dir.create("results", showWarnings = FALSE)
fwrite(res, "results/mr_results.csv")
message("\n已写出 results/mr_results.csv  (", nrow(res), " 行)")

# ---- 概览 ----
message("\n===== 各基因集可用基因数 =====")
print(res[, .(n_genes = uniqueN(gene)), by = .(outcome, set_role)][order(outcome, set_role)])

# ---- 阳性对照检查 (管线试金石) ----
pc <- res[set_role == "positive_control"]
if (nrow(pc) > 0) {
  message("\n===== 阳性对照 (Tier5, 全部方法) =====")
  print(pc[order(pval)][1:min(15, .N),
        .(gene, outcome, method, b, se, pval, FDR, nsnp_used, harmonize_status)])
  message("\n判定: LRRK2/SNCA (PD) 与 PSEN1/PSEN2/APP (AD) 等已知基因应有强信号;")
  message("      若无信号, 说明列名映射/效应方向/harmonize 有误, 先修管线。")
} else {
  warning("阳性对照为空")
}

# ---- 主分析 ----
ms <- res[set_role == "main"]
if (nrow(ms) > 0) {
  message("\n===== 主分析 前 20 行 (按 p 排序) =====")
  print(ms[order(pval)][1:min(20, .N),
        .(gene, outcome, method, b, se, pval, FDR, nsnp_used, harmonize_status)])
  sig <- res[set_role == "main" & FDR < 0.05]
  if (nrow(sig) > 0) {
    fwrite(sig, "results/main_significant.csv")
    message("\nFDR < 0.05 的主分析结果: ", nrow(sig), " 条, 已写出 results/main_significant.csv")
  } else {
    message("\n主分析暂无 FDR < 0.05 结果。阴性结果同样可发表, 但须先确认阳性对照有信号。")
  }
}

# ---- 异质性与多效性 ----
het <- res[!is.na(Q_pval) & Q_pval < 0.05]
if (nrow(het) > 0) message("\nCochran Q p < 0.05 (异质性): ", nrow(het), " 条")
ple <- res[!is.na(egger_pval) & egger_pval < 0.05]
if (nrow(ple) > 0) message("Egger 截距 p < 0.05 (定向多效性): ", nrow(ple), " 条")
