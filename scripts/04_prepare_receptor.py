#!/usr/bin/env python3
"""
04_prepare_receptor.py — 受体准备：找结构、去水去配体、定对接盒子、转 PDBQT
输入: 01_protein.json
输出: 04_receptor.json + receptor.pdb + receptor.pdbqt

策略:
  1. RCSB 搜索 API 按基因名/蛋白名找结构（优先含配体、分辨率低）
  2. 下载 PDB，用 gemmi 去水、去原配体（保留金属离子 ZN/MG/CA/FE/MN）
  3. 定盒子：有原配体→原配体质心 + max-min + 2*padding；无→盲对接(整蛋白)
  4. 转 PDBQT（Vina 列格式）
"""
import argparse
import json
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (RCSB_SEARCH, ALPHAFOLD_FILES, http_get, write_json, read_json,
                    get_outdir, log)

try:
    import gemmi
    GEMMI = True
except Exception:
    GEMMI = False
    log("gemmi 不可用，受体处理将受限", level="warning")

METALS = {"ZN", "MG", "CA", "FE", "MN", "CU", "NI", "CO", "NA", "K"}
WATER = {"HOH", "WAT", "DOD"}

# 结晶添加剂/缓冲液成分：不作为原配体（不用于定盒），也不应保留在受体中
# （但 SO4/PO4 等有时是活性中心组分，故仅排除在"配体定盒"之外，仍按 HETATM 存留与否走金属逻辑）
CRYO_ADDITIVES = {
    # 缓冲液/盐/小溶剂
    "SO4", "PO4", "CL", "BR", "IOD", "NO3", "ACT", "ACY", "FMT", "MES", "EPE",
    "TRS", "GOL", "EDO", "PEG", "PG4", "PGE", "DMS", "MPD", "BME", "DTT",
    "NHE", "HEP", "IMD", "CIT", "TAR", "MLA", "SIN", "OGA", "1PE", "P6G",
    "MOH", "EOH", "IPA", "ACE", "DOD", "AZI", "SCN", "BCT", "CO3",
    # 糖基/糖链（N-糖基化修饰，非活性口袋配体）
    "NAG", "NDG", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA", "FUC", "FUL",
    "XYS", "XYP", "SIA", "NGA", "A2G", "BM3", "RAM", "RIB", "MMA",
    # 常见脂类/去垢剂/载体
    "BOG", "LMT", "LDA", "DDM", "LMN", "PLM", "MYR", "OLA", "STE",
}

# 定盒时盒子单边最大尺寸（Å）——避免误将多个分散配体合并导致超大盒子
MAX_BOX_EDGE = 30.0

ELEM_MAP = {
    "C": "C", "N": "N", "O": "OA", "S": "SA", "H": "H",
    "ZN": "Zn", "MG": "Mg", "CA": "Ca", "FE": "Fe", "MN": "Mn",
    "CU": "Cu", "NI": "Ni", "CO": "Co", "NA": "Na", "K": "K",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
}


# ---------------- 结构检索 ----------------
def rcsb_search(query, rows=10):
    """返回候选结构列表 [(pdb_id, score)]"""
    payload = {
        "query": {"type": "terminal", "service": "full_text",
                  "parameters": {"value": query}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows}},
    }
    try:
        r = http_get(RCSB_SEARCH, params={"json": json.dumps(payload)}, timeout=30)
        data = r.json()
        out = []
        for rs in data.get("result_set", []):
            out.append({"pdb_id": rs["identifier"], "score": rs.get("score")})
        return out
    except Exception as e:
        log(f"[04] RCSB 搜索失败: {e}", level="warning")
        return []


def pick_best_pdb(candidates):
    """优先选含真实配体、分辨率低的结构（用 RCSB entry API 查详情）

    打分：分辨率越低越好；含"真实药物样配体"额外加分（用于定盒）；
    含金属离子（金属酶活性中心）小幅加分。
    """
    scored = []
    for c in candidates[:10]:
        pid = c["pdb_id"]
        try:
            r = http_get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}", timeout=20)
            d = r.json()
            res = (d.get("rcsb_entry_info", {}).get("resolution_combined") or [9.9])[0]
            n_lig = d.get("rcsb_entry_info", {}).get("nonpolymer_entity_count", 0)
            # 非聚合物实体组成（残基名），用于识别真实配体 vs 添加剂
            comp_ids = (d.get("rcsb_entry_info", {})
                        .get("nonpolymer_bound_components") or [])
            real_lig = [x for x in comp_ids
                        if x.upper() not in CRYO_ADDITIVES and x.upper() not in METALS]
            has_real = 1 if real_lig else 0
            score = res - 1.5 * has_real - 0.2 * n_lig
            scored.append((score, pid, res, n_lig, has_real, real_lig))
        except Exception:
            scored.append((9.9, pid, None, 0, 0, []))
    if not scored:
        return None
    scored.sort()
    best = scored[0]
    log(f"[04] 选定结构: {best[1]} (分辨率={best[2]}, 配体数={best[3]}, "
        f"真实配体={best[5] or '无'})")
    return best[1]


# ---------------- 结构处理 ----------------
def process_receptor(pdb_path, out_pdb, keep_metals=True):
    """去水去配体（保留金属），返回 (原配体残基名列表, 原配体坐标)

    说明：仅将"真实药物样配体"（非结晶添加剂、非金属）视为原配体，
    用于确定对接盒子，避免 SO4/CL/EPE 等分散添加剂把盒子撑大。
    """
    st = gemmi.read_structure(str(pdb_path))
    st.setup_entities()
    if len(st) > 0:
        st[0].remove_waters()

    # 收集原配体：按残基名分组统计原子数，排除添加剂/金属/水
    from collections import defaultdict
    groups = defaultdict(list)   # resname -> [coords...]
    for model in st:
        for chain in model:
            for res in chain:
                het = res.het_flag == "H"
                if res.name in WATER or res.name in METALS:
                    continue
                if not het:
                    continue
                if res.name in CRYO_ADDITIVES:
                    continue
                for at in res:
                    groups[res.name].append([at.pos.x, at.pos.y, at.pos.z])

    ligand_resnames = sorted(groups.keys(), key=lambda k: -len(groups[k]))
    ligand_coords = []
    if groups:
        # 优先选"药物样"配体：碳/杂原子丰富的口袋配体，而非长链脂/糖链
        # 用原子数近似，但排除已知脂类/载体（已在 CRYO_ADDITIVES 中）
        main_res = ligand_resnames[0]
        ligand_coords = groups[main_res]
        log(f"[04] 原配体候选 {dict((k, len(v)) for k, v in groups.items())}；"
            f"定盒用 {main_res}（{len(ligand_coords)} 原子）")

    # 删除原配体（保留金属，删除添加剂）
    st.remove_ligands_and_waters()
    if keep_metals:
        st2 = gemmi.read_structure(str(pdb_path))
        st2[0].remove_waters()
        for model in st2:
            for chain in model:
                i = 0
                while i < len(chain):
                    if chain[i].het_flag == "H" and chain[i].name not in METALS:
                        del chain[i]
                    else:
                        i += 1
        st = st2

    st.remove_hydrogens()
    st.write_pdb(str(out_pdb))
    log(f"[04] 受体写出 {out_pdb}；原配体 {ligand_resnames}")
    return ligand_resnames, np.array(ligand_coords) if ligand_coords else None


def compute_box(ligand_coords, padding=8.0, min_size=22.0, max_edge=MAX_BOX_EDGE):
    if ligand_coords is not None and len(ligand_coords) > 0:
        center = ligand_coords.mean(axis=0)
        size = ligand_coords.max(axis=0) - ligand_coords.min(axis=0) + 2 * padding
        size = np.maximum(size, min_size)
        size = np.minimum(size, max_edge)   # 上限，防止超大盒子
        return center.tolist(), size.tolist(), "ligand"
    return None, None, None


def blind_box(pdb_path, padding=6.0, max_edge=MAX_BOX_EDGE):
    """无原配体：兜底定盒策略

    1) 若有催化金属离子（Zn/Ca/Mg/Fe/Mn…）→ 以其质心为盒心（金属酶活性中心）
    2) 否则用整蛋白几何中心
    盒子单边做上限裁剪，避免超大盒子导致内存/时间爆炸。
    """
    st = gemmi.read_structure(str(pdb_path))
    if len(st) > 0:
        st[0].remove_waters()
    all_coords, metal_coords = [], []
    for model in st:
        for chain in model:
            for res in chain:
                for at in res:
                    p = [at.pos.x, at.pos.y, at.pos.z]
                    all_coords.append(p)
                    if res.het_flag == "H" and res.name in METALS:
                        metal_coords.append(p)

    anchor = np.array(metal_coords) if metal_coords else np.array(all_coords)
    center = ((anchor.max(axis=0) + anchor.min(axis=0)) / 2).tolist()
    # 尺寸：以锚点范围 + padding，并限制上下限
    span = anchor.max(axis=0) - anchor.min(axis=0) + 2 * padding
    size = np.minimum(np.maximum(span, 22.0), max_edge).tolist()
    if metal_coords:
        log(f"[04] 无原配体，以金属离子({len(metal_coords)}个)为盒心（金属酶活性口袋）")
    else:
        log("[04] 无原配体、无金属，采用蛋白中心 + 上限裁剪的盲对接盒子", level="warning")
    return center, size


def pdb_to_pdbqt(pdb_in, pdbqt_out):
    """按 Vina 列格式转 PDBQT"""
    lines = []
    n = 0
    for line in open(pdb_in):
        if line.startswith(("ATOM", "HETATM")):
            elem = line[76:78].strip().upper()
            if not elem:
                name = line[12:16].strip()
                elem = "".join(c for c in name if c.isalpha())[:2].upper()
            atype = ELEM_MAP.get(elem, elem.capitalize())
            base = line[:66].ljust(66)
            charge = "0.000".rjust(6)
            lines.append(base + "    " + charge + " " + atype + "\n")
            n += 1
    Path(pdbqt_out).write_text("".join(lines))
    log(f"[04] PDBQT 写出 {pdbqt_out}（{n} 原子）")


def main():
    ap = argparse.ArgumentParser(description="受体准备")
    ap.add_argument("--outdir", default="./tcm_run")
    ap.add_argument("--pdb", help="手动指定 PDB ID（跳过自动检索）")
    ap.add_argument("--center", help="手动指定盒子中心 'x,y,z'")
    ap.add_argument("--size", help="手动指定盒子尺寸 'x,y,z'")
    args = ap.parse_args()
    out = get_outdir(args.outdir)

    info = read_json(out / "01_protein.json")
    gene = info.get("name_for_herb") or info.get("gene")
    pname = info.get("recommended_name")
    ac = info.get("uniprot_ac")

    result = {"gene": gene, "pdb_id": None, "source": None,
              "box_center": None, "box_size": None, "box_mode": None,
              "receptor_pdb": None, "receptor_pdbqt": None, "errors": []}

    # --- 1. 找 PDB 结构 ---
    pdb_id = args.pdb
    if not pdb_id:
        cands = rcsb_search(gene, rows=10)
        if not cands and pname:
            cands = rcsb_search(pname, rows=10)
        if cands:
            pdb_id = pick_best_pdb(cands)

    pdb_file = out / "receptor_raw.pdb"
    if pdb_id:
        try:
            r = http_get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=40)
            pdb_file.write_text(r.text)
            result["pdb_id"] = pdb_id
            result["source"] = f"PDB:{pdb_id}"
            log(f"[04] 下载 PDB {pdb_id} ({len(r.text)} 字节)")
        except Exception as e:
            result["errors"].append(f"PDB 下载失败: {e}")
            pdb_id = None

    # --- 2. 兜底 AlphaFold ---
    if not pdb_id and ac:
        try:
            url = f"{ALPHAFOLD_FILES}/AF-{ac}-F1-model_v4.pdb"
            r = http_get(url, timeout=40)
            if r.status_code == 200 and r.text.startswith(("HEADER", "ATOM", "PARENT", "MODEL", "data")):
                pdb_file.write_text(r.text)
                result["source"] = f"AlphaFold:{ac}"
                log(f"[04] 下载 AlphaFold 结构 {ac}")
        except Exception as e:
            result["errors"].append(f"AlphaFold 失败: {e}")

    if not result["source"]:
        result["errors"].append("未能获取任何蛋白结构，跳过对接")
        write_json(out / "04_receptor.json", result)
        log("[04] 未获得结构，终止", level="error")
        return 1

    # --- 3. 处理受体 + 定盒子 ---
    rec_pdb = out / "receptor.pdb"
    ligands, lig_coords = process_receptor(pdb_file, rec_pdb)

    if args.center and args.size:
        center = [float(x) for x in args.center.split(",")]
        size = [float(x) for x in args.size.split(",")]
        mode = "manual"
    else:
        center, size, mode = compute_box(lig_coords)
        if center is None:
            center, size = blind_box(rec_pdb)
            mode = "blind"
            log("[04] 无原配体，采用盲对接（整蛋白边界框）", level="warning")

    result.update({
        "box_center": center, "box_size": size, "box_mode": mode,
        "receptor_pdb": str(rec_pdb),
        "ligand_resnames": ligands,
    })
    log(f"[04] 盒子中心={[round(x,2) for x in center]} 尺寸={[round(x,1) for x in size]} 模式={mode}")

    # --- 4. 转 PDBQT ---
    rec_pdbqt = out / "receptor.pdbqt"
    pdb_to_pdbqt(rec_pdb, rec_pdbqt)
    result["receptor_pdbqt"] = str(rec_pdbqt)

    write_json(out / "04_receptor.json", result)
    log(f"[04] 完成：受体={result['source']} 盒子模式={mode}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
