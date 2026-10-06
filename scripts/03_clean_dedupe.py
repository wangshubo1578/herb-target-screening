#!/usr/bin/env python3
"""
03_clean_dedupe.py — 化合物清洗、去重、SMILES 补全、（可选）精简
输入: 02_raw_compounds.json
输出: 03_clean_compounds.json — 清洗后的候选化合物（含 SMILES、来源、排序）

清洗三级:
  L1 硬黑名单（内源物/氨基酸/离子溶剂/合成药/维生素辅因子）
  L2 物性过滤（RDKit：非有机元素、MW范围、原子数、无环非芳香）
  L3 来源证据（至少出现在一个中药库）
去重: 名称归一 + InChIKey
SMILES: PubChem name->property 补全
精简: 默认全保留；--top N 时排序取 2N 后 MurckoScaffold 去重到 N
"""
import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (PUBCHEM_BASE, http_get, http_post_json, write_json, read_json,
                    get_outdir, log, normalize_name, safe_step)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, QED, inchi
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDKIT = True
except Exception:
    RDKIT = False
    log("RDKit 不可用，L2 过滤将跳过", level="warning")

BLACKLIST_PATH = Path(__file__).parent.parent / "assets" / "filter_blacklist.json"
ALIAS_PATH = Path(__file__).parent.parent / "assets" / "compound_aliases.json"

# 加载别名表（规范化 key -> 标准名）
try:
    _alias_raw = json.loads(ALIAS_PATH.read_text(encoding="utf-8"))
    ALIASES = {k: v for k, v in _alias_raw.items() if not k.startswith("_")}
except Exception:
    ALIASES = {}

ALLOWED_ELEMENTS = {"C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"}

# 名称常见脏数据修正（normalize_name 后的 key -> 标准名）
NAME_FIXES = {
    "apigeni": "apigenin",
    "1 coffee acid": "caffeic acid",
    "1coffee acid": "caffeic acid",
    "adeninenucleoside": "adenosine",
    "backuchiol": "bakuchiol",
    "bakuchiol": "bakuchiol",
    "campsiol": "campesterol",
    "nordihydroguaiareticacid": "nordihydroguaiaretic acid",
}


def dedup_key(name):
    """生成用于去重的名称键：小写、去括号内容、去连字符/空格/标点、去手性前缀。

    例: '(-)-epigallocatechin-3-gallate' 与 'epigallocatechin 3-gallate(egcg)'
        -> 'epigallocatechin3gallate'
    """
    if not name:
        return ""
    n = name.lower()
    n = re.sub(r"\([^)]*\)", " ", n)          # 去括号（含 egcg / -)- 等）
    n = re.sub(r"[^a-z0-9]+", "", n)          # 仅保留字母数字
    return n


def clean_name(name):
    """去除多余空格、修正常见拼写错误、别名归并"""
    if not name:
        return name
    n = re.sub(r"\s+", " ", name.strip())
    n = n.strip(" ,;")
    key = normalize_name(n)
    if key in NAME_FIXES:
        return NAME_FIXES[key]
    # 别名归并
    std = ALIASES.get(key)
    if std:
        return std
    # 去重键归并（处理 EGCG 等各种写法）
    dk = dedup_key(n)
    if dk in ALIASES:
        return ALIASES[dk]
    return n


# 非中药/非典型小分子的关键词（命中即剔除）
EXCLUDE_KEYWORDS = [
    # 神经递质/内源物
    "noradrenaline", "norepinephrine", "adrenaline", "epinephrine", "dopamine",
    "serotonin", "histamine", "acetylcholine", "gaba", "melatonin",
    # 非法/成瘾/非中药
    "tetrahydrocannabinol", "cannabidiol", "cannabinol", "nicotine", "cotinine",
    "morphine", "codeine", "heroin",
    # 简单醛/小分子溶剂
    "benzaldehyde", "benzoic aldehyde", "cinnamaldehyde",
    # 非单体（混合物/多糖/提取物）
    "polyphenol", "tannin", "extract", "coixan", "polysaccharide", "glycoside mixture",
    "tea polyphenols", "vanilloid", "oleovitamin",
    # 单个氨基酸/衍生物
    "hydroxyproline", "proline", "valine", "leucine", "alanine",
    # 醇类/溶剂残余（避免 "xx ethyl alcohol" 这类脏名）
    "ethyl alcohol", "ethanol", "methanol", "acetone", "glycerol",
]


def load_blacklist():
    bl = json.loads(BLACKLIST_PATH.read_text(encoding="utf-8"))
    flat = set()
    for k, v in bl.items():
        if k.startswith("_"):
            continue
        for x in v:
            flat.add(normalize_name(x))
    return flat


def l1_pass(name, bl):
    """L1：黑名单 + 关键词排除"""
    key = normalize_name(name)
    if key in bl:
        return False, "L1_blacklist"
    for kw in EXCLUDE_KEYWORDS:
        if normalize_name(kw) in key:
            return False, f"L1_keyword:{kw}"
    return True, "ok"


def l2_check(smiles):
    """返回 (是否通过, 原因)"""
    if not RDKIT:
        return True, "rdkit_skipped"
    if not smiles:
        return False, "no_smiles"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "invalid_smiles"
    atoms = list(mol.GetAtoms())
    for a in atoms:
        if a.GetSymbol() not in ALLOWED_ELEMENTS:
            return False, f"disallowed_element:{a.GetSymbol()}"
    mw = Descriptors.MolWt(mol)
    if mw < 60:
        return False, "mw_lt_60"
    if mw > 900:
        return False, "mw_gt_900"
    if len(atoms) < 5:
        return False, "too_few_atoms"
    n_rings = mol.GetRingInfo().NumRings()
    has_aromatic = any(a.GetIsAromatic() for a in atoms)
    if n_rings == 0 and not has_aromatic:
        return False, "acyclic_nonaromatic"
    return True, "ok"


def pubchem_smiles(name):
    """通过 PubChem 名称查 SMILES + 基础性质"""
    try:
        url = f"{PUBCHEM_BASE}/compound/name/{name}/property/CanonicalSMILES,IsomericSMILES,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON"
        r = http_get(url, timeout=20)
        props = r.json()["PropertyTable"]["Properties"][0]
        # PubChem 实际返回字段名为 SMILES/ConnectivitySMILES（非请求里的 IsomericSMILES）
        smi = (props.get("SMILES") or props.get("ConnectivitySMILES")
               or props.get("IsomericSMILES") or props.get("CanonicalSMILES"))
        return {
            "smiles": smi,
            "mw": props.get("MolecularWeight"),
            "xlogp": props.get("XLogP"),
            "tpsa": props.get("TPSA"),
            "hbd": props.get("HBondDonorCount"),
            "hba": props.get("HBondAcceptorCount"),
            "cid": props.get("CID"),
        }
    except Exception:
        return None


def murcko_scaffold(smiles):
    if not RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return None


def rank_key(c):
    return (
        0 if c.get("ob_dl_pass") else 1,
        -(c.get("db_evidence_count", 1)),
        -(c.get("herb_count", 0)),
        -round(c.get("qed", 0) or 0, 3),
        c.get("name", ""),
    )


def main():
    ap = argparse.ArgumentParser(description="化合物清洗与去重")
    ap.add_argument("--outdir", default="./tcm_run")
    ap.add_argument("--top", type=int, default=0,
                    help="精简到 N 个（0=不精简，全部保留）")
    ap.add_argument("--no-smiles", action="store_true",
                    help="跳过 PubChem SMILES 补全（离线测试用）")
    args = ap.parse_args()
    out = get_outdir(args.outdir)

    raw = read_json(out / "02_raw_compounds.json")
    bl = load_blacklist()

    # 1) 汇总所有记录
    records = []
    for src in ("HERB", "TCMSP"):
        for ing in raw.get(src, {}).get("ingredients", []):
            rec = dict(ing)
            rec.setdefault("source_db", src)
            records.append(rec)

    log(f"[03] 原始记录 {len(records)} 条")

    # 2) L1 黑名单 + 名称归一去重
    removed = []
    seen_names = {}
    kept = []
    for rec in records:
        name = clean_name(rec.get("name"))
        if not name:
            continue
        rec["name"] = name
        ok, reason = l1_pass(name, bl)
        if not ok:
            removed.append({"name": name, "reason": reason})
            continue
        key = dedup_key(name)
        if key in seen_names:
            # 合并来源
            seen_names[key]["cross_sources"] = list(
                set(seen_names[key].get("cross_sources", []) + [rec.get("source_db", "")]))
            seen_names[key]["db_evidence_count"] = \
                seen_names[key].get("db_evidence_count", 1) + 1
            continue
        rec["db_evidence_count"] = 1
        seen_names[key] = rec
        kept.append(rec)

    log(f"[03] L1+名称去重后 {len(kept)} 个（剔除 {len(removed)} 个）")

    # 2.5) 大候选集 + 精简模式：先按证据预排序，只对头部补 SMILES
    #      避免对上万个候选逐个查 PubChem（COX2 这类热门靶点可达数千个）
    if args.top and len(kept) > args.top * 5:
        def pre_rank(c):
            return (
                0 if (c.get("ob") and float(c.get("ob") or 0) >= 30
                      and c.get("dl") and float(c.get("dl") or 0) >= 0.18) else 1,
                -(c.get("db_evidence_count", 1)),
                c.get("name", ""),
            )
        kept.sort(key=pre_rank)
        head_n = max(args.top * 12, 200)
        log(f"[03] 候选过多（{len(kept)}）→ 按证据预排序，仅对前 {head_n} 个补 SMILES")
        kept = kept[:head_n]

    # 3a) 并发补全 SMILES（PubChem）
    if not args.no_smiles:
        need = [r for r in kept if not r.get("smiles")]
        if need:
            log(f"[03] 并发查询 PubChem 补全 {len(need)} 个化合物的 SMILES ...")
            from concurrent.futures import ThreadPoolExecutor, as_completed
            with ThreadPoolExecutor(max_workers=8) as ex:
                fut = {ex.submit(pubchem_smiles, r["name"]): r for r in need}
                done = 0
                for f in as_completed(fut):
                    rec = fut[f]
                    pc = safe_step(lambda: f.result(), default=None)
                    if pc and pc.get("smiles"):
                        rec.update({k: v for k, v in pc.items() if k != "cid"})
                        rec["pubchem_cid"] = pc.get("cid")
                        rec["smiles"] = pc.get("smiles")
                    done += 1
                log(f"[03] SMILES 补全完成")

    # 3b) L2 物性过滤 + InChIKey 去重
    #    关键：SMILES 拿不到不致命 —— 保留但标记 no_smiles（跳过后续对接）
    final = []
    seen_inchi = set()
    n_no_smiles = 0
    for rec in kept:
        if rec.get("smiles"):
            ok, reason = l2_check(rec["smiles"])
            if not ok:
                removed.append({"name": rec["name"], "reason": f"L2_{reason}"})
                continue
            ikey = None
            if RDKIT:
                m = Chem.MolFromSmiles(rec["smiles"])
                if m:
                    ikey = inchi.MolToInchiKey(m)
            if ikey and ikey in seen_inchi:
                removed.append({"name": rec["name"], "reason": "dup_inchikey"})
                continue
            if ikey:
                seen_inchi.add(ikey)
                rec["inchikey"] = ikey
            if RDKIT:
                m = Chem.MolFromSmiles(rec["smiles"])
                if m:
                    try:
                        rec["qed"] = round(QED.qed(m), 3)
                    except Exception:
                        pass
            rec["dockable"] = True
        else:
            # 无 SMILES：保留但不可对接（仍可作为候选展示、供人工核查）
            rec["dockable"] = False
            rec["no_smiles"] = True
            n_no_smiles += 1
        final.append(rec)

    log(f"[03] 过滤后 {len(final)} 个（其中无 SMILES 待人工核查 {n_no_smiles} 个）")

    # 4) OB/DL 判定（TCMSP 若有）
    for rec in final:
        ob = rec.get("ob")
        dl = rec.get("dl")
        if ob is not None and dl is not None:
            try:
                rec["ob_dl_pass"] = float(ob) >= 30 and float(dl) >= 0.18
            except Exception:
                rec["ob_dl_pass"] = False
        else:
            rec["ob_dl_pass"] = False

    # 5) 排序
    final.sort(key=rank_key)
    for i, rec in enumerate(final, 1):
        rec["rank"] = i

    # 6) 可选精简
    if args.top and len(final) > args.top:
        # 优先保留可对接（有 SMILES）的化合物，最大化对接产出
        dockable_final = [r for r in final if r.get("dockable")]
        nodock = [r for r in final if not r.get("dockable")]
        # 候选池 = 可对接优先 + 少量待核查（补足展示）
        pool = dockable_final + nodock
        head = pool[:max(args.top * 2, 30)]
        picked, scaffolds = [], set()
        for rec in head:
            if not rec.get("dockable") and len(dockable_final) >= args.top:
                continue   # 已凑够可对接者，跳过无 SMILES 的
            sc = murcko_scaffold(rec.get("smiles")) or rec.get("name")
            if sc in scaffolds:
                continue
            scaffolds.add(sc)
            picked.append(rec)
            if len(picked) >= args.top:
                break
        if len(picked) < args.top:
            for rec in head:
                if rec not in picked:
                    picked.append(rec)
                if len(picked) >= args.top:
                    break
        final = picked
        n_dock = len([r for r in final if r.get("dockable")])
        log(f"[03] 精简至 {len(final)} 个（骨架去重，其中可对接 {n_dock} 个）")

    result = {
        "protein": raw.get("protein"),
        "total_raw": len(records),
        "total_final": len(final),
        "compounds": final,
        "removed": removed,
        "sources": {
            "HERB": len(raw.get("HERB", {}).get("ingredients", [])),
            "TCMSP": len(raw.get("TCMSP", {}).get("ingredients", [])),
        },
    }
    write_json(out / "03_clean_compounds.json", result)
    log(f"[03] 完成：最终候选 {len(final)} 个（剔除 {len(removed)} 个）")
    return 0 if final else 1


if __name__ == "__main__":
    sys.exit(main())
