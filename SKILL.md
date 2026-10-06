---
name: herb-target-screening
description: 输入一个蛋白质名称（基因名/蛋白名），从 HERB、TCMSP 等中药数据库反向检索可靶向该蛋白的中药单体化合物，经清洗去重（剔除内源物、合成药、离子溶剂等非中药特有成分）后，用 AutoDock Vina 批量分子对接筛选结合能最低的化合物，并对最优化合物查询 ADMET 性质、检索相关临床研究，最终输出 Markdown/HTML 报告与数据表。Use when the user wants to find traditional Chinese medicine (TCM) monomer compounds targeting a given protein, perform reverse screening of herbs/herbal ingredients against a target, run molecular docking between herbal compounds and a protein, or obtain ADMET and clinical trial information for hit compounds.
license: MIT
---

# 中药单体反向筛选与分子对接

给定一个**蛋白质名称**，反向找到能靶向它的**中药单体化合物**，并通过分子对接、
ADMET、临床研究逐层筛选，产出可读报告。

> 适用场景：中药/天然产物领域的反向靶点筛选、网络药理学前期命中化合物发现、
> 「某蛋白 → 哪些中药成分可能作用」的探索。**结果为计算预测，不构成实验结论。**

## 什么时候用这个 skill

- 用户给出蛋白名（如 `MMP1`、`COX2`、`TNF`），想找靶向它的中药成分
- 需要「中药单体和某蛋白的分子对接」并排出结合能顺序
- 需要某个命中化合物的 ADMET 性质、相关临床研究

## 全流程（8 步）

```
用户提供蛋白名
   │
 ① 蛋白名标准化        UniProt → gene / recommended_name / accession   → 01_protein.json
   │
 ② 多库检索中药单体    HERB API + TCMSP 抓取                          → 02_raw_compounds.json
   │
 ③ 清洗去重           L1黑名单/关键词 → 名称归一 → PubChem补SMILES
   │                   → L2物性过滤 → InChIKey去重 → (可选)精简       → 03_clean_compounds.json
   │
 ④ 受体准备            RCSB找结构(AlphaFold兜底) → 去水去配体(留金属)
   │                   → 定对接盒子 → 转PDBQT                          → 04_receptor.json
   │
 ⑤ 批量分子对接        SMILES→3D→PDBQT → Vina → 结合能排序            → 05_docking.json
   │
 ⑥ ADMET(最优化合物)   PubChem物性 + RDKit规则                        → 06_admet.json
   │
 ⑦ 临床研究(最优化合物) ClinicalTrials.gov v2                         → 07_clinical.json
   │
 ⑧ 报告生成           Markdown + HTML(Plotly柱状图)                   → report.md / report.html
```

## 依赖

```bash
# Python 包
pip install requests rdkit meeko vina gemmi numpy
```

- `rdkit` — SMILES/物性/骨架
- `meeko` — 配体 PDBQT（0.8.0：`prepare()` 返回 list）
- `vina` — AutoDock Vina **Python API**（不是命令行）
- `gemmi` — PDB 结构处理

必需数据源：HERB、TCMSP、PubChem、UniProt、RCSB/AlphaFold、ClinicalTrials.gov（均需联网）。

## 用法

### 方式一：一键全流程（推荐）

```bash
python scripts/run_pipeline.py --protein MMP1 --outdir ./tcm_run
```

参数：
- `--protein` （必填）蛋白名/基因名
- `--outdir`  输出目录（默认 `./tcm_run`）
- `--top N`   精简到 N 个化合物（默认 0 = 全部对接）
- `--pdb XXXX` 手动指定 PDB 结构
- `--workers N` 对接并行进程数（默认自动 ≤8）

### 方式二：分步执行

```bash
python scripts/01_standardize_protein.py --outdir ./tcm_run --protein MMP1
python scripts/02_fetch_compounds.py    --outdir ./tcm_run
python scripts/03_clean_dedupe.py       --outdir ./tcm_run          # 可加 --top 15
python scripts/04_prepare_receptor.py   --outdir ./tcm_run          # 可加 --pdb 1HFC
python scripts/05_batch_docking.py      --outdir ./tcm_run --workers 8
python scripts/06_admet.py              --outdir ./tcm_run
python scripts/07_clinical_trials.py    --outdir ./tcm_run
python scripts/08_report.py             --outdir ./tcm_run
```

## 与用户的交互约定（重要）

运行前，如果用户未明确，可用 `AskUserQuestion` 确认以下选项：

| 选项 | 默认 | 说明 |
|---|---|---|
| 化合物数量 | **全部对接，不精简** | 用户若要精简，用 `--top 15` |
| 数据源 | **HERB + TCMSP** | TCMID/LTM-TCM 当前不可抓取 |
| 受体结构 | **PDB 优先，AlphaFold 兜底** | 可传 `--pdb` 指定 |
| 对接盒子 | **自动（原配体→兜底盲对接）** | 可传 `--center/--size` |

> 若检索到的化合物**过多**（如 >50），主动建议用户用 `--top 15` 精简
> （按排序取前 2N 再按 Murcko 骨架去重，兼顾活性预测与结构多样性）。

## 产出物

| 文件 | 内容 |
|---|---|
| `report.html` | **交互式报告**（暗色主题 + Plotly 结合能横向柱状图） |
| `report.md` | Markdown 版报告（适合复制/二次编辑） |
| `05_docking.json` | 全部化合物的结合能、构象、排名 |
| `03_clean_compounds.json` | 清洗后候选化合物（含剔除原因） |
| `06_admet.json` / `07_clinical.json` | 最优化合物的 ADMET / 临床研究 |
| `docking/pose_*.pdbqt` | 每个化合物的最优对接构象 |

报告包含 5 节：靶点信息 → 检索与清洗 → 对接结果表 → 最优化合物 ADMET → 相关临床研究。

## 关键实现要点（避坑）

1. **TCMSP 内联 JSON**：正则需要 `re.S`（JSON 跨行）。
   `dataSource : { data : [...] }`（检索页）、`var mol_data = [...]`（明细页）。
2. **PubChem 字段名**：请求 `CanonicalSMILES` 时返回体字段实际叫 **`SMILES`**，
   必须 `props.get("SMILES") or props.get("ConnectivitySMILES")`。
3. **Meeko 0.8.0 API**：`prepare()` 返回 **list**，需遍历后 `PDBQTWriterLegacy.write_string()`。
4. **受体 PDBQT 列格式**：`line[:66]` + 4空格 + 电荷(6右对齐) + 1空格 + 原子类型，
   否则 Vina 报 `bad_conversion` / `Charge not valid`。
5. **受体保留金属离子**（ZN/MG/CA/FE/MN…），它们是酶活中心组分。
6. **SMILES 失败不致命**：保留候选并标 `dockable=False`，跳过对接但仍在报告展示。
7. **并行文件名用 MD5(name)[:8]**，勿用内置 `hash()`（跨进程随机化会撞名）。
8. **TCMID / LTM-TCM 不可用**：仅保留 HERB 返回的交叉引用 ID，不主动抓取。

## 参考资料

- `references/api_notes.md` — 各数据库 API / 抓取方式与可用性
- `references/filter_rules.md` — 化合清洗三级过滤与精简规则
- `references/docking_protocol.md` — AutoDock Vina 对接协议与参数
- `assets/filter_blacklist.json` — 非中药成分黑名单
- `assets/compound_aliases.json` — 化合物别名归并表

## 免责声明

所有结合能（对接打分）、ADMET 与临床信息均为**计算/数据库预测**，属于研究线索，
需要湿实验验证，不构成任何诊断或用药建议。
