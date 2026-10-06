# 分子对接协议 (AutoDock Vina)

## 工具链

| 环节 | 工具 | 版本要点 |
|---|---|---|
| 配体 3D 生成 | RDKit ETKDGv3 + MMFF | `AllChem.EmbedMolecule` |
| 配体 PDBQT | Meeko | 0.8.0：`prepare()` 返回 **list** → `PDBQTWriterLegacy.write_string()` |
| 对接引擎 | AutoDock Vina | **Python API** `from vina import Vina`（非 CLI） |
| 受体处理 | gemmi | 去水/去配体/保留金属 |
| 并行 | `ProcessPoolExecutor` | 多进程，避免 GIL |

## 1. 受体准备 (`04_prepare_receptor.py`)

1. **检索结构**：RCSB Search API 按基因名 → 蛋白名依次尝试。
2. **择优**：RCSB entry API 取分辨率与配体数，打分 `resolution - 0.5 * n_ligands`，
   取最优（分辨率低、有共晶配体者优先）。
3. **兜底**：RCSB 无结果 → AlphaFold DB（`AF-{ac}-F1-model_v4.pdb`）。
4. **清理**：`st[0].remove_waters()` 去水；删除 HETATM 原配体，**但保留金属离子**
   （ZN/MG/CA/FE/MN/CU/NI/CO/NA/K）——它们是许多酶活性中心的必需组分。
5. **转 PDBQT**（Vina 列格式，见下）。

### 盒子（grid box）计算

| 模式 | 触发条件 | 计算 |
|---|---|---|
| `ligand` | 有共晶原配体 | 中心=原配体质心；尺寸=max−min+2×8 Å，最小 22 Å |
| `blind` | 无原配体 | 整蛋白边界框 + 2×4 Å padding |
| `manual` | 用户传 `--center/--size` | 直接使用 |

> 用户可选 **自动+智能兜底**：默认按原配体定盒，无配体则盲对接。

## 2. 配体准备与对接 (`05_batch_docking.py`)

```
SMILES
  → RDKit MolFromSmiles → AddHs
  → ETKDGv3 EmbedMolecule (randomSeed=42) → MMFF 优化
  → Meeko MoleculePreparation.prepare() → PDBQTWriterLegacy.write_string()
  → Vina.set_ligand_from_file()
```

### Vina 参数

| 参数 | 值 | 说明 |
|---|---|---|
| `sf_name` | `vina` | Vina 打分函数 |
| `seed` | 42 | 可复现 |
| `exhaustiveness` | 8 | 计算量/精度平衡 |
| `n_poses` | 9 | 输出构象数 |

### 关键：PDBQT 列格式（受体转换）

Vina 对受体 PDBQT 的列对齐**非常严格**，格式错误的典型报错是
`"bad_conversion"` 或 `"Charge not valid"`。

正确格式 = **66 列**（沿用原 PDB 的 `line[:66]` 并 `ljust(66)`）
+ **4 空格**
+ **电荷（6 右对齐）**
+ **1 空格**
+ **原子类型**

```python
base   = line[:66].ljust(66)
charge = "0.000".rjust(6)
lines.append(base + "    " + charge + " " + atype + "\n")
```

元素→类型映射见 `04_prepare_receptor.py` 的 `ELEM_MAP`
（O→OA, S→SA, Zn→Zn, …）。

### 常见报错与修复

| 报错 | 原因 | 修复 |
|---|---|---|
| `ligand outside grid box` | 配体坐标不在盒子内 | 确认盒子中心/尺寸；必要时平移到中心 |
| `bad_conversion` / `Charge not valid` | 受体 PDBQT 列未对齐 | 按上表 66+4+6+1 格式重排 |
| Meeko `prepare()` 返回不是 str | 0.8.0 API 变更 | 遍历 list，逐个 `write_string()` |
| 配体 `embed_failed` | 3D 嵌入失败 | **重试 `useRandomCoords=True`** |

### 并行

`ProcessPoolExecutor(max_workers=8)`；文件名用 **MD5(name)[:8]** 做稳定 ID
（勿用内置 `hash()`——跨进程随机化会撞名）。

**性能实测**：单个约 4 s；15 个约 1 min；29 个 8 并行约 2 min。

## 3. 结果

`05_docking.json` 结构：
```json
{
  "receptor_source": "PDB:1HFC",
  "box_center": [...], "box_size": [...], "box_mode": "ligand",
  "exhaustiveness": 8,
  "n_docked": 24, "n_failed": 0,
  "best_compound": "dieckol", "best_affinity": -8.417,
  "results": [ {"name","affinity","energies","pose_file","rank","source_db","mw","xlogp","cross_sources"} ],
  "failed": []
}
```
按 `affinity` 升序（越低结合越强）排序，`rank=1` 为最优。
