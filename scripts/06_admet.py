#!/usr/bin/env python3
"""
06_admet.py — 对结合能最低的化合物查询 ADMET 信息
输入: 05_docking.json
输出: 06_admet.json

数据源（ADMETlab 不可用，改用）:
  - PubChem : MW / XLogP / TPSA / HBD / HBA
  - RDKit   : Lipinski / Veber / QED / PAINS 警报（本地计算）
"""
import argparse
import json
import sys
from pathlib import Path
from urllib.parse import quote

sys.path.insert(0, str(Path(__file__).parent))
from common import PUBCHEM_BASE, http_get, write_json, read_json, get_outdir, log, safe_step

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, QED, Lipinski, rdMolDescriptors
    from rdkit.Chem import FilterCatalog
    RDKIT = True
except Exception:
    RDKIT = False


def pubchem_props(name, cid=None):
    """查询 PubChem 理化性质。优先用 CID（更稳），失败再退回名称查询。"""
    props = ("MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,"
             "HBondAcceptorCount,RotatableBondCount")
    urls = []
    if cid:
        urls.append(f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{props}/JSON")
    urls.append(f"{PUBCHEM_BASE}/compound/name/{quote(name)}/property/{props}/JSON")
    for url in urls:
        for attempt in range(3):
            try:
                r = http_get(url, timeout=25)
                return r.json()["PropertyTable"]["Properties"][0]
            except Exception as e:
                log(f"[06] PubChem 查询失败 ({name}, 试{attempt+1}): {e}", level="warning")
    return None


def pubchem_cid(name):
    try:
        r = http_get(f"{PUBCHEM_BASE}/compound/name/{quote(name)}/cids/JSON", timeout=20)
        return r.json()["IdentifierList"]["CID"][0]
    except Exception:
        return None


def rdkit_props(smiles):
    if not RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    out = {
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotBonds": Lipinski.NumRotatableBonds(mol),
        "QED": round(QED.qed(mol), 3),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
        "Rings": rdMolDescriptors.CalcNumRings(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }
    # Lipinski 五规则
    violations = 0
    if out["MW"] > 500:
        violations += 1
    if out["LogP"] > 5:
        violations += 1
    if out["HBD"] > 5:
        violations += 1
    if out["HBA"] > 10:
        violations += 1
    out["LipinskiViolations"] = violations
    out["LipinskiPass"] = violations <= 1
    # Veber 规则（口服可及）
    out["VeberPass"] = out["RotBonds"] <= 10 and out["TPSA"] <= 140
    # PAINS 警报
    try:
        params = FilterCatalog.FilterCatalogParams()
        for cat in (FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
                    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
                    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C):
            params.AddCatalog(cat)
        catalog = FilterCatalog.FilterCatalog(params)
        out["PAINS"] = catalog.HasMatch(mol)
    except Exception:
        out["PAINS"] = None
    # 血脑屏障经验判断（TPSA/LogP 粗判）
    if out["TPSA"] < 90 and 1 < out["LogP"] < 3.5:
        out["BBB_hint"] = "可能透过 (TPSA<90, 1<LogP<3.5)"
    elif out["TPSA"] > 120:
        out["BBB_hint"] = "较难透过 (TPSA>120)"
    else:
        out["BBB_hint"] = "不确定"
    return out


def main():
    ap = argparse.ArgumentParser(description="ADMET 查询")
    ap.add_argument("--outdir", default="./tcm_run")
    ap.add_argument("--compound", help="指定化合物名（默认用对接最优者）")
    args = ap.parse_args()
    out = get_outdir(args.outdir)

    dock = read_json(out / "05_docking.json")
    comps = read_json(out / "03_clean_compounds.json")

    name = args.compound or dock.get("best_compound")
    if not name:
        log("[06] 无最优化合物，跳过 ADMET", level="warning")
        write_json(out / "06_admet.json", {"error": "no_best_compound"})
        return 1

    # 找 smiles
    smiles = None
    for c in comps.get("compounds", []):
        if c.get("name") == name:
            smiles = c.get("smiles")
            break

    log(f"[06] 查询最优化合物 ADMET: {name}")
    cid = safe_step(pubchem_cid, name, default=None)
    pc = safe_step(pubchem_props, name, cid=cid, default=None)
    rd = safe_step(rdkit_props, smiles, default=None)

    result = {"compound": name, "smiles": smiles, "pubchem_cid": cid,
              "pubchem": pc, "rdkit": rd,
              "best_affinity": dock.get("best_affinity")}
    write_json(out / "06_admet.json", result)
    log(f"[06] 完成：{name} — Lipinski {'通过' if (rd or {}).get('LipinskiPass') else '不通过'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
