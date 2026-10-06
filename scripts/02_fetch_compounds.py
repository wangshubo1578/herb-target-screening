#!/usr/bin/env python3
"""
02_fetch_compounds.py — 多库检索靶向某蛋白的中药单体
输入: 01_protein.json
输出: 02_raw_compounds.json — 各库命中的化合物原始列表

数据源（实测可用）:
  - HERB  : POST http://herb.ac.cn/chedi/api/  (search_api -> detail_api)
  - TCMSP : GET  https://old.tcmsp-e.com/tcmspsearch.php (正则提取内联 JSON)
"""
import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (HERB_API, TCMSP_SEARCH, TCMSP_TARGET, TCMSP_TOKEN,
                    http_get, http_post_json, write_json, read_json,
                    get_outdir, log, normalize_name)


# ============ HERB ============
def herb_search_target(gene, aliases):
    """search_api 找靶点 ID；用基因名和别名依次尝试"""
    candidates = [gene] + [a for a in aliases if a and a != gene]
    for kw in candidates:
        try:
            r = http_post_json(HERB_API, {"func_name": "search_api",
                                          "keyword": kw, "label": "Target"}, timeout=30)
            data = r.json().get("res_data", [])
            if len(data) > 1:
                header, rows = data[0], data[1:]
                # 优先精确匹配 Gene name 列
                for row in rows:
                    gene_name = row[1] if len(row) > 1 else ""
                    if normalize_name(gene_name) == normalize_name(gene):
                        tid = row[0]["title"] if isinstance(row[0], dict) else row[0]
                        log(f"[HERB] 靶点命中: {gene_name} -> {tid} (keyword={kw})")
                        return tid, gene_name
        except Exception as e:
            log(f"[HERB] search_api 失败 (keyword={kw}): {e}", level="warning")
    return None, None


def herb_fetch_ingredients(target_id):
    """detail_api 取关联成分（含跨库来源）"""
    try:
        r = http_post_json(HERB_API, {"func_name": "detail_api",
                                      "key_id": target_id, "label": "Target"}, timeout=40)
        d = r.json()
    except Exception as e:
        log(f"[HERB] detail_api 失败: {e}", level="warning")
        return [], []

    def parse_table(key):
        tbl = d.get(key, [])
        if len(tbl) <= 1:
            return [], []
        header, rows = tbl[0], tbl[1:]
        return header, rows

    ing_header, ing_rows = parse_table("ingredient_target")
    herb_header, herb_rows = parse_table("herb_target")

    ingredients = []
    # ingredient_target 列: [ {link,title}=ID, name, [ {link,title}=跨库来源 ] ]
    for row in ing_rows:
        if not row or len(row) < 2:
            continue
        idcell = row[0]
        ing_id = idcell.get("title") if isinstance(idcell, dict) else str(idcell)
        name = row[1]
        sources = []
        if len(row) > 2 and isinstance(row[2], list):
            for s in row[2]:
                if isinstance(s, dict):
                    sources.append(s.get("title", ""))
                else:
                    sources.append(str(s))
        ingredients.append({
            "name": name, "source_db": "HERB", "herb_id": ing_id,
            "cross_sources": sources,
        })

    herbs = []
    for row in herb_rows:
        if not row or len(row) < 2:
            continue
        idcell = row[0]
        hid = idcell.get("title") if isinstance(idcell, dict) else str(idcell)
        herbs.append({"herb_id": hid, "herb_name": row[1],
                      "pvalue": row[2] if len(row) > 2 else None,
                      "fdr": row[3] if len(row) > 3 else None})

    log(f"[HERB] 靶点 {target_id}: {len(ingredients)} 个成分, {len(herbs)} 味中药")
    return ingredients, herbs


# ============ TCMSP ============
def tcmsp_search_target(names):
    """用候选名（推荐全名/别名）精确匹配 TCMSP 靶点，返回数字ID"""
    pat = re.compile(r"dataSource\s*:\s*\{\s*data\s*:\s*(\[.*?\])\s*,", re.S)
    for q in names:
        try:
            r = http_get(TCMSP_SEARCH, params={"q": q, "qs": "target_name",
                                               "token": TCMSP_TOKEN}, timeout=30)
            m = pat.search(r.text)
            if not m:
                continue
            try:
                data = json.loads(m.group(1))
            except Exception:
                continue
            if data:
                row = data[0]
                tid = row.get("target_ID") or row.get("TAR_ID", "").replace("TAR", "")
                log(f"[TCMSP] 靶点命中: '{q}' -> ID={tid} ({row.get('target_name')})")
                return str(tid), row.get("target_name", q)
        except Exception as e:
            log(f"[TCMSP] 搜索失败 (q={q}): {e}", level="warning")
    return None, None


def tcmsp_fetch_ingredients(target_id):
    """靶点页正则提取 var mol_data"""
    try:
        r = http_get(TCMSP_TARGET, params={"qt": target_id}, timeout=30)
        m = re.search(r"var\s+mol_data\s*=\s*(\[.*?\])\s*;", r.text, re.S)
        if not m:
            log(f"[TCMSP] 靶点页无 mol_data (qt={target_id})", level="warning")
            return []
        mols = json.loads(m.group(1))
    except Exception as e:
        log(f"[TCMSP] 靶点页失败: {e}", level="warning")
        return []

    ingredients = []
    for mol in mols:
        ingredients.append({
            "name": mol.get("molecule_name"),
            "source_db": "TCMSP",
            "herb_id": mol.get("molecule_ID"),
            "ob": mol.get("OB"),          # 口服生物利用度(如提供)
            "dl": mol.get("DL"),          # 类药性
            "cross_sources": [],
        })
    log(f"[TCMSP] 靶点 {target_id}: {len(ingredients)} 个成分")
    return ingredients


# ============ 主流程 ============
def main():
    ap = argparse.ArgumentParser(description="多库检索中药单体")
    ap.add_argument("--outdir", default="./tcm_run")
    args = ap.parse_args()
    out = get_outdir(args.outdir)

    info = read_json(out / "01_protein.json")
    gene = info.get("name_for_herb") or info.get("gene")
    tcmsp_name = info.get("name_for_tcmsp")
    aliases = info.get("name_aliases", [])

    result = {"protein": {"gene": gene, "uniprot_ac": info.get("uniprot_ac"),
                          "name_for_tcmsp": tcmsp_name},
              "HERB": {"target_id": None, "ingredients": [], "herbs": []},
              "TCMSP": {"target_id": None, "ingredients": []},
              "errors": []}

    # --- HERB ---
    try:
        tid, gname = herb_search_target(gene, aliases)
        if tid:
            ings, herbs = herb_fetch_ingredients(tid)
            result["HERB"] = {"target_id": tid, "target_name": gname,
                              "ingredients": ings, "herbs": herbs}
        else:
            result["errors"].append("HERB: 未找到该靶点")
    except Exception as e:
        result["errors"].append(f"HERB 失败: {e}")
        log(f"HERB 整体失败: {e}", level="warning")

    # --- TCMSP ---
    try:
        names = [tcmsp_name] + [a for a in aliases if a != tcmsp_name]
        tid, tname = tcmsp_search_target(names)
        if tid:
            ings = tcmsp_fetch_ingredients(tid)
            result["TCMSP"] = {"target_id": tid, "target_name": tname,
                               "ingredients": ings}
        else:
            result["errors"].append("TCMSP: 未找到该靶点（可尝试其它别名）")
    except Exception as e:
        result["errors"].append(f"TCMSP 失败: {e}")
        log(f"TCMSP 整体失败: {e}", level="warning")

    total = len(result["HERB"]["ingredients"]) + len(result["TCMSP"]["ingredients"])
    write_json(out / "02_raw_compounds.json", result)
    log(f"[02] 完成：HERB {len(result['HERB']['ingredients'])} + "
        f"TCMSP {len(result['TCMSP']['ingredients'])} = 共 {total} 条原始记录")

    if total == 0:
        log("未检索到任何成分，请检查靶点名称或稍后重试", level="warning")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
