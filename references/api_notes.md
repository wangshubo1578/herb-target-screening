# 数据库 API / 抓取笔记

本 skill 检索中药单体化合物的数据源。下表为各库可用性实测结论。

## 可用性总览

| 数据库 | 状态 | 接入方式 | 备注 |
|---|---|---|---|
| **HERB** (herb.ac.cn) | ✅ 可用 | RPC HTTP API | 主数据源，返回较规范 |
| **TCMSP** (old.tcmsp-e.com) | ✅ 可用 | HTML 内联 JSON 抓取 | 次数据源，含 OB/DL |
| TCMID | ❌ 不可用 | — | 旧域名重定向，新 IP 连接失败 |
| LTM-TCM (cloud.tasly.com) | ❌ 不可用 | — | 无响应 / 需登录 |
| PubChem | ✅ 可用 | PUG-REST | 补全 SMILES 与物性 |
| UniProt | ✅ 可用 | REST | 蛋白名标准化 |
| RCSB PDB | ✅ 可用 | Search API v2 | 受体结构首选 |
| AlphaFold DB | ✅ 可用 | 文件直链 | 受体结构兜底 |
| ClinicalTrials.gov | ✅ 可用 | API v2 | 临床研究检索 |

> **TCMID / LTM-TCM 目前不可抓取。** skill 在检索时仍会在 `cross_sources`
> 中保留 HERB 返回的 TCMID 交叉引用 ID（如 `TCMID:5483`），但不会主动抓取这两个库。

---

## 1. HERB API

- **端点**: `POST http://herb.ac.cn/chedi/api/`
- **Content-Type**: `application/json`
- **请求体**: `{"func_name": "<函数名>", ...参数}`

### 1.1 靶点检索 → 成分

```
func_name = "search_api"
payload   = {"func_name": "search_api", "search_type": "target", "keyword": "<gene_name>", ...}
```

返回靶点匹配列表。取 `target_id` 后再调 `detail_api` 拿成分-靶点关系。

### 1.2 取成分明细

```
func_name = "detail_api"
```

返回 `ingredient_target`（成分→靶点）、`herb_target`（药材→靶点）等键。

### 1.3 关键字段

- `ingredient_name` / `ingredient_name_cn` — 成分名（英文/中文）
- `herb_name` — 来源药材
- `symmap_id` / `tcmid_id` — 交叉引用 ID
- 分子量、XLogP 等物性

> ⚠️ HERB 成分名脏数据较多，如 `"1- coffee acid"`、`"apigeni"`、`"adeninenucleoside"`。
> 清洗步骤 `03_clean_dedupe.py` 的 `clean_name()` + `compound_aliases.json` 负责修正。

---

## 2. TCMSP 抓取

- **检索 URL**: `https://old.tcmsp-e.com/tcmspsearch.php`
- **靶点 URL**: `https://old.tcmsp-e.com/target.php`
- **Token**: `fa6966e547446646375e7c2a176163ab`（实测恒定，首页不含 token）

### 2.1 靶点检索页

页面内联 JS 变量：
```js
dataSource : { data : [ { ... }, ... ] }
```
正则（**必须 `re.S`**，JSON 跨行）：
```python
re.search(r"dataSource\s*:\s*\{\s*data\s*:\s*(\[.*?\])\s*,", html, re.S)
```

### 2.2 成分明细页

页面内联：
```js
var mol_data = [ { ... }, ... ] ;
```
正则：
```python
re.search(r"var\s+mol_data\s*=\s*(\[.*?\])\s*;", html, re.S)
```

### 2.3 关键字段

- `MOL_ID` (如 MOL000006)
- `molecule_name`
- `ob` — 口服生物利用度 (%)
- `dl` — 类药性
- 用于排序过滤：`OB ≥ 30` 且 `DL ≥ 0.18` 为经典筛选阈值

---

## 3. PubChem PUG-REST

- **Base**: `https://pubchem.ncbi.nlm.nih.gov/rest/pug`

```
GET /compound/name/{name}/property/CanonicalSMILES,IsomericSMILES,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON
```

> ⚠️ **坑**：请求 `IsomericSMILES`/`CanonicalSMILES` 时，返回体里的字段名实际是
> **`SMILES`** 或 **`ConnectivitySMILES`**。必须按此读取：
> ```python
> smi = props.get("SMILES") or props.get("ConnectivitySMILES") or props.get("IsomericSMILES")
> ```

名称支持同义词匹配（如 `meletin` → quercetin），但**返回的可能是同义名对应的结构**，
因此仍需本地别名归并做去重。

---

## 4. UniProt REST

```
GET https://rest.uniprot.org/uniprotkb/search?query={name}&format=json&size=5
```

用于蛋白名标准化：
- `gene` — 基因名（HERB 检索用）
- `recommended_name` — 推荐全名（TCMSP 检索用）
- `aliases` — 别名列表

例：`MMP1` → P03956，gene=MMP1，recommended_name=`Interstitial collagenase`。

---

## 5. RCSB PDB Search API v2

```
POST https://search.rcsb.org/rcsbsearch/v2/query
```

按 UniProt accession 检索实验结构：
```json
{"query":{"type":"terminal","service":"text","parameters":{"attribute":"rcsb_polymer_entity_container_identifiers.uniprot_accession","operator":"exact_match","value":"P03956"}},"return_type":"entry","request_options":{"results_content_type":["experimental"],"paginate":{"start":0,"rows":25}}}
```
按分辨率、配体数等排序，取最佳 PDB。失败时兜底 AlphaFold。

---

## 6. AlphaFold DB

```
GET https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb
```

RCSB 无合适结构时使用。注意 AlphaFold 模型无共晶配体 → 盒子需按全蛋白或
预测口袋计算（`04_prepare_receptor.py` 的自动兜底逻辑）。

---

## 7. ClinicalTrials.gov API v2

```
GET https://clinicaltrials.gov/api/v2/studies?query.intr={intervention}&pageSize=20
```

按干预物（化合物名）检索临床研究，取：
- `protocolSection.identificationModule.nctId`
- `briefTitle` / `officialTitle`
- `overallStatus`
- `phase`
- `conditions`
