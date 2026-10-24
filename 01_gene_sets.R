# A1 方案 · 基因集定义
# 与《A1_眼表营养轴MR_基因集与研究方案_20261022.md》保持一致
# 分层设计: Tier1 锚点 / Tier2-3 眼表营养轴 / Tier4 炎症对照 / Tier5 阳性对照

library(data.table)

gene_sets <- list(

  # ---- Tier 1: MVP 干眼 GWAS 直接命中位点 ----
  # 依据 PMID 42790339 (EBioMedicine 2026;132:106490)
  # 注意: 同研究预印本 (PMID 40166553) 此处写作 EPHA5, 与正文 EFNA5 不一致, 投稿前须核对
  Tier1_DED_GWAS = c("EFNA5", "EPHA5", "GRIA1", "SYNGAP1", "BLK"),

  # ---- Tier 2: 角膜神经再生与维持 ----
  # 依据 PMID 40940733 (Cells 2025 综述), 41186355 (IOVS 2025), 34132748, 28910445
  Tier2a_neurotrophin = c(
    "NGF", "NTRK1", "NGFR", "BDNF", "NTF3", "NTF4", "GDNF", "CNTF", "SERPINF1", "WNT5A"
  ),
  Tier2b_axon_guidance = c(
    "NTN1", "SLIT3", "SEMA3A", "SEMA3F", "DCC", "UNC5C", "ROBO1", "EPHA4", "EPHA5", "EFNA5"
  ),
  Tier2c_neuropeptide = c(
    "TAC1", "TACR1", "CALCA", "CALCRL", "VIP", "NPY", "SST", "GAL"
  ),
  Tier2d_sensory_channel = c(
    "TRPV1", "TRPM8", "TRPA1", "SCN9A", "SCN10A", "P2RX3"
  ),
  Tier2e_other = c("AQP5"),

  # ---- Tier 3: 眼表神经肽 (Tier2c 的强化子集) ----
  # 依据 PMID 35149082 (Exp Eye Res 2022)
  Tier3_ocular_neuropeptide = c(
    "TAC1", "CALCA", "VIP", "NPY", "SST", "GAL", "NTS", "TACR1", "CALCRL", "RAMP1"
  ),

  # ---- Tier 4: 泪液炎症轴 (对照集) ----
  # 对应课题组实测面板: IL-4 / IFN-gamma / IgE / ECP / MMP-9
  # 注意: ECP 通常由 RNASE3 编码, EDN 由 RNASE2 编码, 请按试剂盒说明书确认
  Tier4_inflammation = c(
    "IL4", "IL4R", "IFNG", "IFNGR1", "MMP9", "RNASE2", "RNASE3",
    "IL1B", "IL6", "IL6R", "TNF", "IL10", "ICAM1", "IL13"
  ),

  # ---- Tier 5: 神经退行核心基因 (阳性对照, 用于验证管线) ----
  Tier5_neurodegeneration = c(
    "MAPT", "SNCA", "APP", "PSEN1", "PSEN2", "APOE", "LRRK2", "GBA",
    "PINK1", "PRKN", "C9orf72", "SOD1", "TARDBP", "FUS", "HTT", "FXN"
  )
)

# ---- 合并为主分析集与对照集 ----
# 主分析: Tier1 + Tier2 + Tier3 (去重)
main_set <- unique(c(
  gene_sets$Tier1_DED_GWAS,
  unlist(gene_sets[grep("^Tier2", names(gene_sets))]),
  gene_sets$Tier3_ocular_neuropeptide
))

# 对照: Tier4 炎症轴
control_inflammation <- gene_sets$Tier4_inflammation

# 阳性对照: Tier5
positive_control <- gene_sets$Tier5_neurodegeneration

# ---- 输出 ----
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)

gene_table <- data.table::rbindlist(lapply(names(gene_sets), function(t) {
  data.table::data.table(tier = t, gene = gene_sets[[t]])
}))
gene_table[, set_role := fifelse(
  tier %like% "^Tier[123]", "main",
  fifelse(tier == "Tier4_inflammation", "control", "positive_control")
)]
data.table::fwrite(gene_table, "data/processed/gene_table.csv")

message("基因集统计:")
print(gene_table[, .(n_genes = .N), by = .(tier, set_role)])
message("主分析集: ", length(main_set), " 个基因")
message("已写出 data/processed/gene_table.csv")
