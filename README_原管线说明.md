# A1 分析管线 · 执行说明（本机实测版）

配套文档：`A1_眼表营养轴MR_基因集与研究方案_20261022.md`

---

## 一、本机环境（已实测确认）

| 项目 | 状态 |
|---|---|
| R 版本 | **4.6.0 (2026-04-24 ucrt)** |
| R 路径 | `D:\program filies\r\R-4.6.0\bin\Rscript.exe` |
| data.table | 1.18.4 已装 |
| dplyr | 1.2.1 已装 |
| ggplot2 | 4.0.3 已装 |
| **coloc** | 5.2.3 已装 |
| **susieR** | 0.14.2 已装 |
| ieugwasr | 1.2.0.9000 已装 |
| **TwoSampleMR** | **装不上**（GitHub API 不通，重试两次均失败） |

**应对方案**：MR 估计量已自实现于 `mr_core.R`，不依赖 TwoSampleMR。

---

## 二、mr_core.R 自实现的正确性验证

在 R 4.6.0 下用模拟数据验证（真实因果效应设为 0.5）：

| 检验项 | 结果 |
|---|---|
| IVW 固定效应（单次） | 0.5022 |
| **IVW 200 次重复均值** | **0.4997**（真实 0.5，标准差 0.0101）→ 无偏 |
| MR-Egger | 0.5270，截距 -0.0075，p = 0.258 → 正确识别无多效性 |
| 加权中位数 | 0.4935 |
| Wald ratio（单 SNP） | 正常 |
| 等位基因 harmonize | direct 保留 / swapped 正确翻转 / 回文**正确剔除** |

实现内容：Wald ratio、IVW（固定+随机效应）、MR-Egger、加权中位数（bootstrap SE）、Cochran Q、Egger 截距多效性、Steiger 定向。

---

## 三、数据源可达性实测结果

| 数据源 | 状态 | 说明 |
|---|---|---|
| eQTLGen 主文件 | **可下载 322.8 MB** | 域名已从 molgenis26 迁至 `download.gcc.rug.nl` |
| eQTLGen MAF | **可下载 240 MB** | 必需，用于 Z-score → beta 换算 |
| GTEx v8 | 403 Forbidden | 当前不可达，待换镜像或换入口 |
| UKB-PPP 门户 | 502 | 当前不可达 |
| EBI GWAS Catalog FTP | 可达 | AD 目录存在，具体文件需再定位 |
| pdgenetics.org | 可达 | PD GWAS 来源 |
| imsgc.net | 502 | MS GWAS 需换源 |

---

## 四、eQTLGen 的关键限制（重要）

官网确认 eQTLGen **只提供 Z-score，不提供 Beta/SE**，列名为：
`SNP, Gene, P-value, SNP Chr/Pos(hg19), Gene Symbol, Z-score, Assessed, Other, Nr Cohorts, Nr Samples, FDR`

因此必须走换算：

```
SE   = 1 / sqrt(2 * MAF * (1 - MAF) * N)
beta = Zscore * SE
```

MAF 需从单独文件下载：
`https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz`

⚠️ 这是标准化效应量的近似换算，**方法学部分必须写明公式**。若追求更严谨，应改用 GTEx（自带 slope/slope_se），但 GTEx 当前不可达。

正确文件名（已实测）：
`2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz`

---

## 五、执行顺序

```bash
cd "D:/workbuddy working/A1_眼表营养轴MR_分析管线"
RS="D:/program filies/r/R-4.6.0/bin/Rscript.exe"

$RS 01_gene_sets.R             # 基因集（已跑通）
#   手动下载 data/raw/ 数据（见第三节）
$RS 02_extract_instruments.R   # 提取工具变量
$RS _test_mr_core.R            # 算法自检（可选，建议跑）
$RS 03_mr_analysis.R           # MR 主分析
$RS 04_coloc.R                 # 共定位
```

---

## 六、两处必须人工核对

1. **EFNA5 还是 EPHA5**：EBioMedicine 正文（PMID 42790339）写 EFNA5，medRxiv 预印本（PMID 40166553）写 EPHA5。配体-受体关系，结论方向不同，投稿前按正式版 Table 核对。基因集两者都保留。

2. **ECP 基因归属**：ECP 通常 RNASE3，EDN 通常 RNASE2。Tier 4 两者都保留，按试剂盒说明书确认。

---

## 七、管线验证（不要跳过）

`03_mr_analysis.R` 会输出 Tier 5 阳性对照结果。**APOE 等已知基因应在 AD 上显示信号；若无，说明列名映射、效应方向或 harmonize 有误，先修管线。**

---

## 八、已知弱点与应对

| 弱点 | 应对 |
|---|---|
| TwoSampleMR 装不上 | 已用 mr_core.R 自实现替代，验证无偏 |
| eQTLGen 只有 Z-score | 用 MAF 换算 beta，方法学写明公式 |
| 血液 eQTL 代理眼表表达 | GTEx 组织验证（当前不可达，需解决）+ Limitations 说明 |
| 部分基因可能无强 cis-eQTL | 02 脚本会输出未命中基因，需改 pQTL 或剔除 |
| MR-Base API 不稳定 | clumping 走本地 PLINK |
| 结局 GWAS 部分源不可达 | MS 需换源（imsgc 502），可考虑 EBI GWAS Catalog |
