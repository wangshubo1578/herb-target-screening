#!/usr/bin/env python3
"""
07_clinical_trials.py — 查找涉及最优化合物的临床研究
输入: 05_docking.json（+ 03 提供化合物名）
输出: 07_clinical.json

数据源: ClinicalTrials.gov v2 API
  https://clinicaltrials.gov/api/v2/studies?query.intr=<化合物>&pageSize=20
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import CLINICALTRIALS_API, http_get, write_json, read_json, get_outdir, log


def base_names(name):
    """生成用于检索的名称变体：原名 + 去手性前缀/立体描述符的基础名。"""
    variants = [name]
    # 去立体化学前缀，如 (S)-, (R)-, (-)-, (+)-, (2R,3S)-
    stripped = re.sub(r"^\s*\([^)]*\)\s*", "", name)
    stripped = re.sub(r"^\s*[(\[]?[+\-][)\]]?\s*", "", stripped).strip()
    if stripped and stripped != name:
        variants.append(stripped)
    # 去长 IUPAC 命名中的空格后缀（保留首个词，如 "Stylopine hydrochloride" → "Stylopine"）
    first = re.split(r"[\s,]+", stripped or name)[0]
    if first and first not in variants and len(first) > 3:
        variants.append(first)
    return variants


def search_trials(compound, max_results=20):
    """按干预物检索临床研究。依次尝试名称变体，返回命中最多的一组。"""
    best = []
    for cand in base_names(compound):
        attempts = [
            {"query.intr": cand, "pageSize": max_results},
            {"query.term": cand, "pageSize": max_results},
        ]
        for params in attempts:
            try:
                r = http_get(CLINICALTRIALS_API, params=params, timeout=30)
                data = r.json()
            except Exception as e:
                log(f"[07] ClinicalTrials 查询失败 ({cand}): {e}", level="warning")
                continue
            trials = _parse(data)
            # 过滤：干预物中确实含该化合物名（避免全文误命中）
            key = cand.lower()
            trials = [t for t in trials
                      if key in " ".join(t.get("interventions", [])).lower()
                      or key in (t.get("title") or "").lower()]
            if len(trials) > len(best):
                best = trials
            if best:
                return best
    return best


def _parse(data):
    trials = []
    for s in data.get("studies", []):
        p = s.get("protocolSection", {})
        idm = p.get("identificationModule", {})
        st = p.get("statusModule", {})
        dm = p.get("designModule", {})
        cm = p.get("conditionsModule", {})
        im = p.get("armsInterventionsModule", {})
        sm = p.get("sponsorCollaboratorsModule", {})

        interventions = [iv.get("name") for iv in im.get("interventions", []) if iv.get("name")]

        trials.append({
            "nct_id": idm.get("nctId"),
            "title": idm.get("briefTitle"),
            "status": st.get("overallStatus"),
            "phases": dm.get("phases", []),
            "study_type": dm.get("studyType"),
            "start_date": st.get("startDateStruct", {}).get("date"),
            "completion_date": st.get("completionDateStruct", {}).get("date"),
            "conditions": cm.get("conditions", []),
            "interventions": interventions,
            "sponsor": sm.get("leadSponsor", {}).get("name"),
            "enrollment": dm.get("enrollmentInfo", {}).get("count"),
            "url": f"https://clinicaltrials.gov/study/{idm.get('nctId')}",
        })
    return trials


def main():
    ap = argparse.ArgumentParser(description="临床研究检索")
    ap.add_argument("--outdir", default="./tcm_run")
    ap.add_argument("--compound", help="指定化合物名（默认用对接最优者）")
    ap.add_argument("--max", type=int, default=20)
    args = ap.parse_args()
    out = get_outdir(args.outdir)

    dock = read_json(out / "05_docking.json")
    name = args.compound or dock.get("best_compound")
    if not name:
        log("[07] 无最优化合物，跳过临床检索", level="warning")
        write_json(out / "07_clinical.json", {"error": "no_best_compound", "trials": []})
        return 1

    log(f"[07] 检索化合物 '{name}' 的临床研究 ...")
    trials = search_trials(name, max_results=args.max)

    result = {"compound": name, "n_trials": len(trials), "trials": trials}
    write_json(out / "07_clinical.json", result)
    log(f"[07] 完成：找到 {len(trials)} 项临床研究")
    return 0

if __name__ == "__main__":
    sys.exit(main())
