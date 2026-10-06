# herb-target-screening

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Claude%20Code%20%7C%20CodeBuddy-purple.svg)](#%EF%B8%8F-安装)
[![Engine](https://img.shields.io/badge/Docking-AutoDock%20Vina%201.2-ff69b4.svg)](https://vina.scripps.edu/)

> 🌿 **中药单体反向筛选与分子对接 Agent Skill** —— 给定一个蛋白质名称，自动从 HERB / TCMSP
> 检索可靶向它的中药单体化合物，清洗去重后用 AutoDock Vina 批量分子对接，输出结合能排名、
> 最优化合物 ADMET 性质与相关临床研究报告。

An Agent Skill for [Claude Code](https://claude.com/claude-code) / [CodeBuddy](https://www.codebuddy.cn/):
reverse-screen traditional Chinese medicine monomer compounds against a protein target,
with batch molecular docking (AutoDock Vina), ADMET and clinical-trial reporting.

---

## ✨ 功能特性

- 🔍 **多库检索**：HERB API + TCMSP 双数据源，保留 TCMID / SymMap 交叉引用
- 🧹 **智能清洗**：三级过滤（黑名单 / 关键词 / RDKit 物性）+ 别名归并 + InChIKey 结构去重，
  自动剔除异烟肼、酒精等非中药成分
- 🎯 **受体自动准备**：RCSB 择优（优先含真实配体结构）→ AlphaFold 兜底；
  自动排除结晶添加剂与糖基，以真实口袋配体定盒；保留金属酶活性中心
- ⚗️ **批量对接**：AutoDock Vina Python API，8 进程并行，进程池崩溃自动降级串行
- 💊 **ADMET 查询**：PubChem 理化性质 + Lipinski / Veber / PAINS / BBB 经验规则
- 🏥 **临床研究**：ClinicalTrials.gov v2 检索（自动去手性前缀、过滤假阳性）
- 📊 **可视化报告**：暗色主题 HTML（Plotly 交互图）+ Markdown 双格式

## 🧬 全流程

```
蛋白名 → ①标准化(UniProt) → ②多库检索(HERB+TCMSP) → ③清洗去重 → ④受体准备(RCSB/AlphaFold)
       → ⑤批量对接(Vina) → ⑥ADMET → ⑦临床研究 → ⑧报告(HTML+MD)
```

## 🛠️ 环境要求

- Python **3.9+**
- 联网（HERB / TCMSP / PubChem / UniProt / RCSB / ClinicalTrials.gov）

## 📦 安装

### 1. 克隆到 Agent Skills 目录

```bash
# Claude Code
git clone https://github.com/<your-username>/herb-target-screening.git \
    ~/.claude/skills/herb-target-screening

# CodeBuddy
git clone https://github.com/<your-username>/herb-target-screening.git \
    ~/.codebuddy/skills/herb-target-screening
```

### 2. 安装 Python 依赖

```bash
pip install -r herb-target-screening/requirements.txt
```

<details>
<summary>依赖明细</summary>

| 包 | 用途 |
|---|---|
| `requests` | 数据库 HTTP 请求 |
| `rdkit` | SMILES / 物性 / 骨架去重 |
| `meeko` | 配体转 PDBQT |
| `vina` | AutoDock Vina Python API |
| `gemmi` | PDB 结构处理 |
| `numpy` | 数值计算 |

</details>

## 🚀 使用

对 Agent 说一句话即可（如 *「帮我筛选靶向 COX2 的中药单体」*），
或直接命令行一键运行：

```bash
python scripts/run_pipeline.py --protein COX2 --outdir ./tcm_run
```

| 参数 | 说明 | 默认 |
|---|---|---|
| `--protein` | （必填）蛋白名 / 基因名，如 `MMP1`、`COX2` | — |
| `--outdir` | 输出目录 | `./tcm_run` |
| `--top N` | 精简到 N 个（热门靶点推荐） | `0` 全部对接 |
| `--pdb` | 手动指定 PDB 结构 | 自动检索 |
| `--workers` | 并行进程数 | 自动 ≤8 |

## 📊 产出物

| 文件 | 内容 |
|---|---|
| `report.html` | **交互式报告**（Plotly 结合能排名图 + ADMET + 临床研究表） |
| `report.md` | Markdown 版报告 |
| `05_docking.json` | 全部化合物结合能与构象 |
| `03_clean_compounds.json` | 候选化合物（含剔除原因） |
| `docking/pose_*.pdbqt` | 各化合物最优对接构象 |

## 🧪 验证示例

| 靶点 | 原始记录 | 最终候选 | 最优化合物 | 结合能 |
|---|---|---|---|---|
| MMP1 (P03956) | 91 条 | 27 个 | hypericin（贯叶连翘素） | **-9.51 kcal/mol** |
| COX2 (P35354) | 6841 条 | 15 个 | (S)-Stylopine（四氢小檗碱类） | **-9.49 kcal/mol** |

> COX2 场景：热门靶点自动启用「预排序 + 头部补 SMILES」快速通道，3325 个候选 20 秒内完成清洗精简。

## 📁 项目结构

```
herb-target-screening/
├── SKILL.md            # Agent 加载的 skill 入口
├── README.md           # 本文件
├── LICENSE             # MIT
├── CITATION.cff        # 引用元数据（GitHub 自动渲染引用按钮）
├── scripts/            # 8 步流程脚本 + run_pipeline.py
├── references/         # 数据库 API / 清洗规则 / 对接协议文档
└── assets/             # 黑名单与别名表
```

## 📄 引用

如果您在研究中使用了本项目，欢迎引用：

```bibtex
@software{herb_target_screening_2026,
  title   = {herb-target-screening: TCM monomer reverse screening
             and molecular docking Agent Skill},
  author  = {Your Name},                 % ← 改成你的名字
  year    = {2026},
  url     = {https://github.com/<your-username>/herb-target-screening},
  license = {MIT}
}
```

也可以直接点击仓库首页右侧的 **「Cite this repository」** 按钮获取引用格式
（由 `CITATION.cff` 自动生成）。

## 🤝 贡献

欢迎提交 Issue 与 Pull Request！请先阅读
[CONTRIBUTING.md](CONTRIBUTING.md) 与
[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)。

## ⚠️ 免责声明

所有结合能（对接打分）、ADMET 与临床信息均为**计算 / 数据库预测**，属于研究线索，
需湿实验验证，不构成任何诊断或用药建议。HERB、TCMSP、PubChem 等数据库内容
版权归各自数据库所有，二次分发或商用请遵循其数据使用条款。

## 📄 许可证

[MIT](LICENSE) © 2026 <Your Name>
