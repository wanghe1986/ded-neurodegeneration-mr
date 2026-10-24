#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 Figure 1-4 的绘图 JSON（喂给 academic-figures skill）。"""
import json, os

OUTDIR = r"D:\workbuddy working\figures_MR"
os.makedirs(OUTDIR, exist_ok=True)

# ══════════ Figure 1：研究设计与工具变量构建（diagram）══════════
# 布局铁律：框宽必须按 Arial 9pt 实测墨迹宽 + 2×PAD 估算，先手算再出图，
# 否则文字溢出（2026-10-24 第一版 w=3.4 全部溢出的教训）。
COL = [0.0, 7.6, 15.2]
ROW = [6.6, 4.9, 3.2, 1.4]
W, H = 6.8, 1.4

src = [
    ("src_ded", "Dry eye disease", "MVP · 84,262 / 314,374", "#0072B2"),
    ("src_ad", "Alzheimer disease", "85,934 / 401,577", "#009E73"),
    ("src_pd", "Parkinson disease", "33,674 / 449,056", "#D55E00"),
    ("src_ms", "Multiple sclerosis", "FinnGen · 2,926 / 495,931", "#CC79A7"),
]
ivs = [
    ("iv_ded", "DED instruments", "39 loci · F 21.7-45.4", "#56B4E9"),
    ("iv_ad", "AD instruments", "43 loci", "#56B4E9"),
    ("iv_pd", "PD instruments", "25 loci", "#56B4E9"),
    ("iv_ms", "MS instruments", "35 loci", "#56B4E9"),
]
blocks = []
for i, (bid, lab, sub, col) in enumerate(src):
    blocks.append({"id": bid, "label": lab, "sublabel": sub,
                   "x": COL[0], "y": ROW[i], "w": W, "h": H, "color": col})
for i, (bid, lab, sub, col) in enumerate(ivs):
    blocks.append({"id": bid, "label": lab, "sublabel": sub,
                   "x": COL[1], "y": ROW[i], "w": W, "h": H, "color": col})
blocks.append({"id": "fwd", "label": "Forward MR", "sublabel": "DED -> AD / PD / MS",
               "x": COL[2], "y": 4.9, "w": W, "h": 1.5, "color": "#E69F00"})
blocks.append({"id": "rev", "label": "Reverse MR", "sublabel": "AD / PD / MS -> DED",
               "x": COL[2], "y": 2.8, "w": W, "h": 1.5, "color": "#E69F00"})
blocks.append({"id": "sens", "label": "Sensitivity",
               "sublabel": "LOO · pleiotropy · coloc",
               "x": COL[2], "y": 0.2, "w": W, "h": 1.5, "color": "#8a8a8a"})

arrows = [
    {"from": "src_ded", "to": "iv_ded"},
    {"from": "src_ad", "to": "iv_ad"},
    {"from": "src_pd", "to": "iv_pd"},
    {"from": "src_ms", "to": "iv_ms"},
    {"from": "iv_ded", "to": "fwd"},
    {"from": "iv_ad", "to": "rev"},
    {"from": "iv_pd", "to": "rev", "rad": -0.15},
    {"from": "iv_ms", "to": "rev", "rad": -0.25},
    {"from": "fwd", "to": "sens"},
    {"from": "rev", "to": "sens"},
]
# 列标题用 annotations（groups 的 label 会相互重叠）
CX = [COL[0] + W / 2, COL[1] + W / 2, COL[2] + W / 2]
annotations = [
    {"text": "GWAS data sources", "x": CX[0], "y": 8.6,
     "fontsize": 10, "weight": "bold", "ha": "center", "color": "#444444"},
    {"text": "Instruments", "x": CX[1], "y": 8.6,
     "fontsize": 10, "weight": "bold", "ha": "center", "color": "#444444"},
    {"text": "Two-sample MR", "x": CX[2], "y": 8.6,
     "fontsize": 10, "weight": "bold", "ha": "center", "color": "#444444"},
    {"text": " ", "x": COL[0] - 0.4, "y": -0.4, "fontsize": 8, "color": "#ffffff"},
    {"text": " ", "x": COL[2] + W + 0.4, "y": 8.9, "fontsize": 8, "color": "#ffffff"},
]
fig1 = {"background": "light", "blocks": blocks, "arrows": arrows,
        "annotations": annotations}


# ══════════ Figure 2：双向 MR 森林图（IVW 固定效应）══════════
# 数值取自 results/mr_pairs.csv 的 ivw_fe 行，与稿件 Table 2 / Table 3 一致
fwd = [
    ("DED -> Parkinson disease", 0.6198, 0.4681, 0.8207, 25),
    ("DED -> Alzheimer disease", 1.0348, 0.9417, 1.1372, 32),
    ("DED -> Multiple sclerosis", 1.2815, 0.8987, 1.8274, 25),
]
rev = [
    ("Parkinson disease -> DED", 1.0039, 0.9849, 1.0232, 24),
    ("Alzheimer disease -> DED", 0.9754, 0.9482, 1.0034, 40),
    ("Multiple sclerosis -> DED", 1.0018, 0.9827, 1.0212, 14),
]
order = [fwd[0], fwd[1], fwd[2], rev[0], rev[1], rev[2]]
fig2 = {
    "labels": [f"{n}  (n={k} SNPs)" for n, _, _, _, k in order],
    "estimates": [e for _, e, _, _, _ in order],
    "ci_low": [l for _, _, l, _, _ in order],
    "ci_high": [u for _, _, _, u, _ in order],
    "ref_line": 1.0,
    "measure": "Odds ratio (95% CI)",
}


# ══════════ Figure 3 / 4：单 SNP Wald ratio 与留一法 ═══════════
d = json.load(open("results/figdata_wald_loo.json", encoding="utf-8"))

w = d["wald"]                       # 已按 OR 升序
# log 尺度森林图（MR 标准做法）：25 个 CI 完整无截断，ref_line=0 即 OR=1
import math as _m
fig3 = {
    "labels": [x["snp"] for x in w],
    "estimates": [round(_m.log(x["OR"]), 3) for x in w],
    "ci_low": [round(_m.log(x["lci"]), 3) for x in w],
    "ci_high": [round(_m.log(x["uci"]), 3) for x in w],
    "ref_line": 0.0,
    "measure": "Wald ratio odds ratio (log scale)",
}
# 完整 CI（自然尺度）供补充表使用
json.dump([{"snp": x["snp"], "OR": round(x["OR"], 3),
            "lci": round(x["lci"], 3), "uci": round(x["uci"], 3),
            "p": float("%.3g" % x["p"])} for x in w],
          open(os.path.join(OUTDIR, "fig3_wald_full.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
print("Fig3: log 尺度 CI 范围 %.2f 至 %.2f，无截断"
      % (min(fig3["ci_low"]), max(fig3["ci_high"])))
l = d["loo"]                        # 已按 OR 升序
ov = d["overall"]
fig4 = {
    "labels": [f"Omitting {x['snp']}" for x in l],
    "estimates": [round(x["OR"], 3) for x in l],
    "ci_low": [round(x["lci"], 3) for x in l],
    "ci_high": [round(x["uci"], 3) for x in l],
    "ref_line": 1.0,
    "measure": "Odds ratio (95% CI) after exclusion",
}

for name, obj in [("fig1_diagram", fig1), ("fig2_forest", fig2),
                  ("fig3_wald", fig3), ("fig4_loo", fig4)]:
    p = os.path.join(OUTDIR, name + ".json")
    json.dump(obj, open(p, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("written", p)

print("\n核对：全样本 IVW OR=%.4f (%.4f-%.4f), p=%.5g"
      % (ov["OR"], ov["lci"], ov["uci"], ov["p"]))
print("Wald 点估计范围 %.3f - %.3f ；留一法范围 %.3f - %.3f"
      % (fig3["estimates"][0], fig3["estimates"][-1],
         min(fig4["estimates"]), max(fig4["estimates"])))
print("留一法 CI 范围 %.3f - %.3f" % (min(fig4["ci_low"]), max(fig4["ci_high"])))
