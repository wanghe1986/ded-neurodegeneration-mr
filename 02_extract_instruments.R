# A1 方案 · 从 eQTL / pQTL 提取工具变量
#
# 支持两类输入:
#   格式 A (推荐): 有 beta/se  —— GTEx v8 (slope/slope_se), UKB-PPP, deCODE pQTL
#   格式 B:        只有 Zscore —— eQTLGen 血液 cis-eQTL, 需用 MAF 近似换算
#
# 输出: data/processed/exposure_clumped.csv (已格式化为 TwoSampleMR 暴露格式)

# ================= CONFIG =================
# 选择一个数据源, 取消对应注释

# --- 格式 A1: GTEx v8 (推荐, 有 slope/slope_se, 且多组织可选) ---
QTL_FILE  <- "data/raw/GTEx_v8_Brain_Cortex.signif_variant_gene_pairs.txt.gz"
FMT       <- "gtex"
# GTEx 列名: variant_id, gene_id, tss_distance, maf, slope, slope_se, pval_nominal ...
# GTEx 基因是 Ensembl ID, 需 GENE_MAP 转成 Symbol (见下)

# --- 格式 A2: pQTL (UKB-PPP / deCODE), 有 beta/se ---
# QTL_FILE <- "data/raw/ukbppp_pqtl.txt.gz"
# FMT      <- "pqtl"

# --- 格式 B: eQTLGen (只有 Zscore) ---
# QTL_FILE <- "data/raw/eqtlgen_cis_eqtl.txt.gz"
# FMT      <- "eqtlgen"

P_THRESH   <- 5e-8      # 全基因组显著
CIS_WINDOW <- 1e6       # cis 窗口 (GTEx/pQTL 已限定 cis, eQTLGen 需过滤)
F_MIN      <- 10        # 弱工具变量阈值
CLUMP_R2   <- 0.001
CLUMP_KB   <- 10000
# ==========================================

library(data.table)
source("01_gene_sets.R")

target <- unique(c(main_set, control_inflammation, positive_control))
message("目标基因数: ", length(target))

dat <- fread(QTL_FILE, showProgress = FALSE)
message("读入: ", nrow(dat), " 行")

# ---------- 格式 A1: GTEx ----------
if (FMT == "gtex") {
  # GTEx 的 gene_id 是 Ensembl ID (带版本号), 需转 Symbol
  # 方案1: 用 biomaRt 在线转换 (需网络)
  # 方案2: 用本地 gencode 注释文件 (推荐, 稳定)
  if (!file.exists("data/raw/gencode_gene_id_to_symbol.tsv")) {
    stop("缺少 data/raw/gencode_gene_id_to_symbol.tsv\n",
         "从 GENCODE 下载 gene_id->gene_name 映射, 或用 biomaRt 生成:\n",
         "  mart <- biomaRt::useMart('ensembl','hsapiens_gene_ensembl')\n",
         "  bm <- biomaRt::getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart=mart)\n",
         "  data.table::fwrite(as.data.table(bm), 'data/raw/gencode_gene_id_to_symbol.tsv', sep='\\t')")
  }
  gmap <- fread("data/raw/gencode_gene_id_to_symbol.tsv")
  setnames(gmap, names(gmap), c("gene_id", "gene"))
  # GTEx gene_id 带版本后缀 (.x), 去掉
  dat[, gene_id := sub("\\..*$", "", gene_id)]
  dat <- merge(dat, unique(gmap), by = "gene_id", all.x = TRUE)

  dat <- dat[gene %in% target & pval_nominal < P_THRESH]
  setnames(dat, c("slope", "slope_se", "pval_nominal"),
               c("beta", "se", "pval"), skip_absent = TRUE)
  dat[, n := 670]   # GTEx 各组织样本量不同, 按实际填 (73-706)
  message("GTEx 命中: ", nrow(dat), " 行, ", uniqueN(dat$gene), " 个基因")
}

# ---------- 格式 A2: pQTL ----------
if (FMT == "pqtl") {
  # 按实际文件列名调整
  setnames(dat, c("gene", "effect_allele", "other_allele", "beta", "se", "pval", "n"),
           c("gene", "effect_allele", "other_allele", "beta", "se", "pval", "n"),
           skip_absent = TRUE)
  dat <- dat[gene %in% target & pval < P_THRESH]
  message("pQTL 命中: ", nrow(dat), " 行, ", uniqueN(dat$gene), " 个基因")
}

# ---------- 格式 B: eQTLGen ----------
if (FMT == "eqtlgen") {
  # eQTLGen 列名 (2019-12-11 版):
  # Gene, GeneSymbol, SNP, AssessedAllele, OtherAllele, Zscore, FDR, Pvalue, NrSamples
  dat <- dat[GeneSymbol %in% target & Pvalue < P_THRESH]
  # Zscore -> beta 近似换算 (需要 MAF; 若文件无 MAF, 需从 1000G 补齐)
  if (!"MAF" %in% names(dat)) {
    stop("eQTLGen 文件缺少 MAF 列, 无法从 Zscore 换算 beta。\n",
         "解决办法: (1) 从 1000G 按 rsID 补齐 MAF; (2) 改用 GTEx 数据源 (自带 slope/slope_se)。")
  }
  # SE = 1 / sqrt(2 * MAF * (1-MAF) * N);  beta = Z * SE   (标准化效应量近似)
  dat[, se   := 1 / sqrt(2 * MAF * (1 - MAF) * NrSamples)]
  dat[, beta := Zscore * se]
  setnames(dat, c("GeneSymbol", "AssessedAllele", "OtherAllele", "Pvalue", "NrSamples"),
               c("gene", "effect_allele", "other_allele", "pval", "n"))
  message("eQTLGen 命中: ", nrow(dat), " 行, ", uniqueN(dat$gene), " 个基因")
  message("提示: 此处 beta 为标准化效应量近似值, 方法学部分须写明换算公式。")
}

if (nrow(dat) == 0) stop("没有任何基因命中, 请检查基因名或数据源格式")

# ---------- F 统计量筛选 ----------
dat[, F_stat := (beta^2) / (se^2)]
dat_strong <- dat[F_stat > F_MIN]
message("F > ", F_MIN, " 后: ", nrow(dat_strong), " 行, ",
        uniqueN(dat_strong$gene), " 个基因")
message("未命中基因: ", paste(setdiff(target, unique(dat_strong$gene)), collapse = ", "))

# ---------- LD clumping ----------
# 推荐本地 PLINK (不依赖 MR-Base API):
#   plink --bfile 1000G_EUR --clump-p1 1 --clump-r2 0.001 --clump-kb 10000 \
#         --clump data/processed/for_plink_clump.txt --out data/processed/clumped
# 然后把 clumped.clumped 里的 SNP 列表取回
CLUMPED_FILE <- "data/processed/clumped.clumped"
if (file.exists(CLUMPED_FILE)) {
  cl <- fread(CLUMPED_FILE)
  keep <- unique(cl$SNP)
  dat_clumped <- dat_strong[SNP %in% keep]
  message("PLINK clumping 后: ", nrow(dat_clumped), " 行")
} else {
  warning("未找到 ", CLUMPED_FILE, "\n",
          "请先本地运行 PLINK clumping, 或改用 ieugwasr::ld_clump()。\n",
          "当前输出未 clumping 的结果 (仅用于流程验证, 不可用于投稿)。")
  dat_clumped <- dat_strong
}

# ---------- 格式化为 TwoSampleMR 暴露格式 ----------
library(TwoSampleMR)
exposure_fmt <- format_data(
  as.data.frame(dat_clumped),
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  effect_allele_col = "effect_allele",
  other_allele_col  = "other_allele",
  pval_col          = "pval",
  gene_col          = "gene",
  samplesize_col    = "n",
  id_col            = "gene"
)

fwrite(as.data.table(exposure_fmt), "data/processed/exposure_clumped.csv")
message("\n已写出 data/processed/exposure_clumped.csv")
message("工具变量数 (按基因):")
print(as.data.table(exposure_fmt)[, .(n_iv = .N), by = id.exposure][order(-n_iv)])
