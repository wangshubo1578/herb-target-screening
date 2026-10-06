#!/usr/bin/env python3
"""
05_batch_docking.py — 批量分子对接
输入: 03_clean_compounds.json（可对接的化合物）+ 04_receptor.json（受体 + 盒子）
输出: 05_docking.json — 每个化合物的结合能 + 最优构象

流程: SMILES → RDKit 加氢/3D → Meeko PDBQT → Vina 对接 → 结合能
"""
import argparse
import json
import sys
import os
import hashlib
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, str(Path(__file__).parent))
from common import write_json, read_json, get_outdir, log

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from vina import Vina
    DEPS = True
except Exception as e:
    DEPS = False
    _ERR = str(e)
    log(f"对接依赖缺失: {e}", level="warning")


def smiles_to_pdbqt(smiles, out_path):
    """SMILES → 3D → PDBQT"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "invalid_smiles"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        if AllChem.EmbedMolecule(mol, useRandomCoords=True) != 0:
            return False, "embed_failed"
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    preparator = MoleculePreparation()
    try:
        setups = preparator.prepare(mol)
    except Exception as e:
        return False, f"meeko_prepare:{e}"
    for setup in setups:
        s, ok, err = PDBQTWriterLegacy.write_string(setup)
        if ok:
            Path(out_path).write_text(s)
            return True, "ok"
    return False, "pdbqt_write_failed"


def _stable_id(name):
    """跨进程稳定的短 ID（避免 Python hash 随机化导致文件名冲突）"""
    return hashlib.md5(name.encode("utf-8")).hexdigest()[:8]


def dock_one(task):
    """单化合物对接（子进程）"""
    name, smiles, receptor_pdbqt, center, size, workdir = task
    try:
        sid = _stable_id(name)
        lig_pdbqt = os.path.join(workdir, f"lig_{sid}.pdbqt")
        ok, msg = smiles_to_pdbqt(smiles, lig_pdbqt)
        if not ok:
            return {"name": name, "affinity": None, "error": msg}

        v = Vina(sf_name="vina", seed=42, verbosity=0)
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(lig_pdbqt)
        v.compute_vina_maps(center=list(center), box_size=list(size))
        v.dock(exhaustiveness=8, n_poses=9)
        energies = v.energies(n_poses=9)
        best = float(energies[0][0])

        # 保存最优构象
        pose_path = os.path.join(workdir, f"pose_{sid}.pdbqt")
        v.write_poses(pose_path, n_poses=1, overwrite=True)

        return {"name": name, "affinity": best, "n_poses": len(energies),
                "energies": [float(e[0]) for e in energies], "pose_file": pose_path}
    except Exception as e:
        return {"name": name, "affinity": None, "error": str(e)}


def main():
    ap = argparse.ArgumentParser(description="批量分子对接")
    ap.add_argument("--outdir", default="./tcm_run")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--workers", type=int, default=0, help="并行进程数（0=自动）")
    args = ap.parse_args()
    out = get_outdir(args.outdir)

    if not DEPS:
        log("对接依赖不可用，跳过对接步骤", level="error")
        write_json(out / "05_docking.json", {"error": "deps_missing", "results": []})
        return 1

    comp = read_json(out / "03_clean_compounds.json")
    rec = read_json(out / "04_receptor.json")

    if not rec.get("receptor_pdbqt"):
        log("无受体文件，跳过对接", level="error")
        write_json(out / "05_docking.json", {"error": "no_receptor", "results": []})
        return 1

    dockable = [c for c in comp["compounds"] if c.get("dockable") and c.get("smiles")]
    if not dockable:
        log("无可对接化合物", level="warning")
        write_json(out / "05_docking.json", {"error": "no_ligands", "results": []})
        return 1

    center = rec["box_center"]
    size = rec["box_size"]
    receptor = rec["receptor_pdbqt"]
    workdir = str(out / "docking")
    Path(workdir).mkdir(parents=True, exist_ok=True)

    log(f"[05] 开始对接 {len(dockable)} 个化合物（exhaustiveness={args.exhaustiveness}）")

    n_workers = args.workers or min(8, os.cpu_count() or 4)
    tasks = [(c["name"], c["smiles"], receptor, center, size, workdir) for c in dockable]

    results = []
    try:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futs = {ex.submit(dock_one, t): t[0] for t in tasks}
            done = 0
            for f in as_completed(futs):
                r = f.result()
                results.append(r)
                done += 1
                tag = f"{r['affinity']:.2f} kcal/mol" if r.get("affinity") is not None else f"失败({r.get('error')})"
                log(f"[05] ({done}/{len(tasks)}) {r['name']}: {tag}")
    except Exception as e:
        # 进程池崩溃（多为内存不足 / 子进程被杀）：降级为串行重试
        log(f"[05] 并行对接异常（{type(e).__name__}: {e}），降级为串行执行", level="warning")
        results = []
        for i, t in enumerate(tasks, 1):
            r = dock_one(t)
            results.append(r)
            tag = f"{r['affinity']:.2f} kcal/mol" if r.get("affinity") is not None else f"失败({r.get('error')})"
            log(f"[05] (seq {i}/{len(tasks)}) {r['name']}: {tag}")

    # 按结合能排序（越低越好）
    ok_results = [r for r in results if r.get("affinity") is not None]
    ok_results.sort(key=lambda x: x["affinity"])
    failed = [r for r in results if r.get("affinity") is None]

    # 合并化合物元数据
    meta = {c["name"]: c for c in dockable}
    for i, r in enumerate(ok_results, 1):
        r["rank"] = i
        m = meta.get(r["name"], {})
        r["source_db"] = m.get("source_db")
        r["mw"] = m.get("mw")
        r["xlogp"] = m.get("xlogp")
        r["ob"] = m.get("ob")
        r["dl"] = m.get("dl")
        r["cross_sources"] = m.get("cross_sources", [])

    result = {
        "receptor_source": rec.get("source"),
        "box_center": center, "box_size": size, "box_mode": rec.get("box_mode"),
        "exhaustiveness": args.exhaustiveness,
        "n_docked": len(ok_results), "n_failed": len(failed),
        "best_compound": ok_results[0]["name"] if ok_results else None,
        "best_affinity": ok_results[0]["affinity"] if ok_results else None,
        "results": ok_results,
        "failed": failed,
    }
    write_json(out / "05_docking.json", result)

    if ok_results:
        log(f"[05] 完成：{len(ok_results)} 个成功，最优 {ok_results[0]['name']} "
            f"({ok_results[0]['affinity']:.2f} kcal/mol)")
    return 0 if ok_results else 1


if __name__ == "__main__":
    sys.exit(main())
