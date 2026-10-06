# 化合物清洗与筛选规则

原始检索结果（HERB + TCMSP）包含大量重复、非中药特有成分、脏名与混合物，
需经三级过滤 + 去重后得到可信的中药单体候选集。

数据流：
```
02_raw_compounds.json
  → clean_name()       名称修正 + 别名归并
  → L1 硬黑名单/关键词  剔除内源物、合成药、溶剂、氨基酸…
  → 名称归一去重        同名合并（累加来源证据）
  → PubChem SMILES 补全
  → L2 物性过滤(RDKit)  元素/分子量/原子数/成环性
  → InChIKey 去重       结构级去重
  → 排序 → (可选)精简
03_clean_compounds.json
```

---

## 名称修正 (`clean_name`)

HERB 名称脏数据示例 → 修正：

| 原始 | 修正后 |
|---|---|
| `1- coffee acid` | caffeic acid |
| `apigeni` | apigenin |
| `adeninenucleoside` | adenosine |
| `meletin` | quercetin |
| `campherol` | kaempferol |
| `backuchiol` | bakuchiol |
| `(-)-epigallocatechin-3-gallate` = `epigallocatechin 3-gallate(egcg)` | 归一合并 |

- **拼写修正**：`NAME_FIXES`（脚本内）
- **别名归并**：`assets/compound_aliases.json`
- **去重键 `dedup_key`**：小写 → 去括号 → 仅留字母数字。
  例如 `epigallocatechin-3-gallate` 与 `epigallocatechin 3-gallate (egcg)`
  都归一为 `epigallocatechin3gallate`，从而合并。

---

## L1 硬黑名单 + 关键词

黑名单文件：`assets/filter_blacklist.json`，分类：

| 类别 | 示例 | 理由 |
|---|---|---|
| `endogenous` | 雌二醇、睾酮、皮质醇、葡萄糖 | 内源激素/代谢物，非中药特有 |
| `amino_acids` | 甘氨酸、脯氨酸、色氨酸… | 单氨基酸 |
| `ions_solvents` | 水、乙醇、DMSO、Na⁺、Cl⁻ | 溶剂/离子 |
| `synthetic_drugs` | **异烟肼**、阿司匹林、二甲双胍… | 化学合成药物 |
| `vitamins_cofactors` | 抗坏血酸、NAD、ATP、腺苷… | 维生素/辅因子 |
| `fatty_acids_simple` | 棕榈酸、硬脂酸、山梨醇 | 简单脂/糖醇 |
| `voc_general` | 甲醛、乙酸、乳酸、柠檬酸 | 通用小分子 |
| `heavy_metals_toxins` | 砷、铅、汞、黄曲霉毒素 | 重金属/毒素 |
| `common_metabolites` | 水杨酸、苯甲酸、吲哚、腐胺 | 常见代谢物 |

**关键词子串匹配**（`EXCLUDE_KEYWORDS`，命中即剔除）补充黑名单不足：
- 神经递质：noradrenaline, dopamine, serotonin, histamine…
- 成瘾/非法：tetrahydrocannabinol, cannabidiol, nicotine, morphine…
- 简单醛：benzaldehyde, cinnamaldehyde
- 混合物/多糖：polyphenol, tannin, extract, polysaccharide…
- 单氨基酸：hydroxyproline, valine…
- 醇类溶剂残余：ethyl alcohol, ethanol, glycerol…

> 匹配前统一 `normalize_name`（小写 + 去 `- _ 空格`）。

---

## L2 物性过滤 (RDKit)

对能拿到 SMILES 的结构：

| 规则 | 阈值 | 目的 |
|---|---|---|
| 允许元素 | C H O N S P F Cl Br I | 剔除金属/异常元素 |
| 分子量 | 60 ≤ MW ≤ 900 | 剔除碎片与过大分子 |
| 原子数 | ≥ 5 | 剔除极小分子 |
| 成环性 | 有环 **或** 有芳香环 | 剔除无环脂肪链小分子 |

> **SMILES 拿不到 ≠ 剔除**：保留该化合物但标记 `dockable=False`、`no_smiles=True`，
> 仍列入候选展示、供人工核查，只是跳过对接。避免误杀名字冷门但真实的成分。

---

## 排序

`rank_key`（越小越靠前）：
```
(0 if ob_dl_pass else 1,      # TCMSP OB≥30 且 DL≥0.18 优先
 -db_evidence_count,          # 多库交叉印证优先
 -herb_count,                 # HERB 中出现频次
 -qed,                        # 类药性
 name)
```

## 精简（可选，`--top N`）

> 默认 **不精简**（`--top 0`，全部对接）。用户若要求精简：

1. 先按排序取前 `2N`；
2. 用 **Murcko 骨架** 去重（同骨架系列只留代表）；
3. 不足 `N` 时按序补齐。

兼顾「排序靠前」与「骨架多样性」。

---

## 输出结构

`03_clean_compounds.json`：
```json
{
  "protein": {...},
  "total_raw": 91, "total_final": 29,
  "compounds": [ {"name","smiles","dockable","source_db","cross_sources","db_evidence_count","mw","xlogp","qed","inchikey","rank", ...} ],
  "removed":   [ {"name","reason":"L1_blacklist / L1_keyword:xx / L2_mw_gt_900 / dup_inchikey"} ],
  "sources":   {"HERB": 74, "TCMSP": 17}
}
```
