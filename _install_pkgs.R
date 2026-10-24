# A1 方案 · 依赖安装（本机 R 4.6.0 专用）
# 用法: "D:/program filies/r/R-4.6.0/bin/Rscript.exe" _install_pkgs.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

need_cran <- c("remotes", "coloc", "susieR", "readr", "stringr", "tidyr")
for (p in need_cran) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat("[SKIP] ", p, "\n", sep = "")
  } else {
    cat("[INSTALL] ", p, " ...\n", sep = "")
    try(install.packages(p, dependencies = TRUE))
  }
}

need_github <- c("MRCIEU/ieugwasr", "MRCIEU/TwoSampleMR")
for (p in need_github) {
  nm <- basename(p)
  if (requireNamespace(nm, quietly = TRUE)) {
    cat("[SKIP] ", nm, "\n", sep = "")
  } else {
    cat("[INSTALL-GH] ", p, " ...\n", sep = "")
    try(remotes::install_github(p, upgrade = "never", dependencies = TRUE))
  }
}

cat("\n===== 最终检查 =====\n")
for (p in c("data.table", "dplyr", "coloc", "susieR", "TwoSampleMR", "ieugwasr")) {
  cat(sprintf("%-14s %s\n", p,
      ifelse(requireNamespace(p, quietly = TRUE),
             as.character(packageVersion(p)), "MISSING")))
}
