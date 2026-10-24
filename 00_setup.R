# A1 方案 · 环境搭建
# 眼表神经营养轴 -> AD/PD/MS 的 SMR 与两样本 MR
# 运行方式: Rscript 00_setup.R

# ---- 需要安装的包 ----
cran <- c("data.table", "dplyr", "tidyr", "ggplot2", "readr", "stringr")
bioc <- c("GenomicRanges", "ensembldb")   # 可选，用于基因坐标注释

# TwoSampleMR 与 ieugwasr 来自 MR-Base，需从 GitHub 安装
github <- c("MRCIEU/TwoSampleMR", "MRCIEU/ieugwasr")

# coloc 与 susieR（共定位）
cran <- c(cran, "coloc", "susieR")

for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("安装 CRAN 包: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  } else {
    message("已存在: ", p)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (p in github) {
  nm <- basename(p)
  if (!requireNamespace(nm, quietly = TRUE)) {
    message("安装 GitHub 包: ", p)
    remotes::install_github(p)
  } else {
    message("已存在: ", nm)
  }
}

# ---- 关键提醒 ----
# ieugwasr 依赖 MR-Base 的在线 API (gwas-api.mrcieu.ac.uk)。
# 该接口偶发 SSL 证书过期或限流。本管线主推【本地文件】方式，
# 在线 API 仅用于 LD clumping 等可选步骤。
#
# 若在线 API 不可用, 用本地 PLINK + 1000G EUR 参考面板做 clumping。

# ---- 目录结构 ----
dirs <- c("data/raw", "data/processed", "results", "figs")
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
message("目录已就绪: ", paste(dirs, collapse = ", "))

# ---- 验证 ----
library(TwoSampleMR)
library(data.table)
message("TwoSampleMR 版本: ", as.character(packageVersion("TwoSampleMR")))
message("环境检查完成")
