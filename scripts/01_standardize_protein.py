#!/usr/bin/env python3
"""
01_standardize_protein.py — 蛋白名标准化
输入: --protein <用户输入的名称/基因名/别名>
输出: 01_protein.json — 规范名 + 各数据库所需名称格式 + 别名列表

用途: 用户可能输入基因名(MMP1)、蛋白全名(Interstitial collagenase)、别名(CLG)，
      统查 UniProt 得到标准信息，并生成：
        name_for_herb  = 基因名 (HERB 用)
        name_for_tcmsp = UniProt 推荐全名 (TCMSP 精确匹配用)
        name_aliases   = 所有别名 (用于失败重试)
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import UNIPROT_API, http_get, write_json, get_outdir, log


def query_uniprot(user_input, human_only=True):
    """多字段 OR 查询，提高命中率"""
    clauses = [
        f"gene:{user_input}",
        f'protein_name:"{user_input}"',
        f"id:{user_input}",
    ]
    q = "(" + " OR ".join(clauses) + ")"
    if human_only:
        q += " AND organism_id:9606 AND reviewed:true"
    params = {"query": q, "format": "json", "size": 10}
    r = http_get(UNIPROT_API, params=params, timeout=30)
    return r.json().get("results", [])


def extract(entry):
    genes = entry.get("genes", [])
    primary_gene = None
    syns = []
    if genes:
        primary_gene = genes[0].get("geneName", {}).get("value")
        for g in genes:
            primary_gene = primary_gene or g.get("geneName", {}).get("value")
            for s in g.get("synonyms", []):
                syns.append(s.get("value"))
    pd = entry.get("proteinDescription", {})
    rec = pd.get("recommendedName", {}).get("fullName", {}).get("value")
    alt = []
    for a in pd.get("alternativeNames", []):
        v = a.get("fullName", {}).get("value")
        if v:
            alt.append(v)
    return {
        "uniprot_ac": entry.get("primaryAccession"),
        "uniprot_id": entry.get("uniProtkbId"),
        "gene": primary_gene,
        "gene_synonyms": syns,
        "recommended_name": rec,
        "alternative_names": alt,
        "organism": entry.get("organism", {}).get("scientificName"),
        "sequence_length": entry.get("sequence", {}).get("length"),
        "annotation_score": entry.get("annotationScore", 0),
    }


def main():
    ap = argparse.ArgumentParser(description="蛋白名标准化")
    ap.add_argument("--protein", required=True, help="蛋白质名称/基因名/别名")
    ap.add_argument("--outdir", default="./tcm_run")
    ap.add_argument("--no-human", action="store_true", help="不限制人类来源")
    args = ap.parse_args()

    out = get_outdir(args.outdir)
    log(f"查询 UniProt: {args.protein}")

    results = query_uniprot(args.protein, human_only=not args.no_human)
    if not results:
        log("UniProt 无命中，降级：放宽查询条件", level="warning")
        results = query_uniprot(args.protein, human_only=False)

    if not results:
        info = {
            "user_input": args.protein,
            "uniprot_ac": None,
            "gene": args.protein,
            "recommended_name": args.protein,
            "name_for_herb": args.protein,
            "name_for_tcmsp": args.protein,
            "name_aliases": [args.protein],
            "ambiguous": False,
            "note": "UniProt 未命中，直接使用用户输入",
        }
    else:
        results.sort(key=lambda e: e.get("annotationScore", 0), reverse=True)
        info = extract(results[0])
        info["user_input"] = args.protein
        info["ambiguous"] = len(results) > 1
        info["candidates"] = [extract(e) for e in results]

        # 各库所需名称
        info["name_for_herb"] = info.get("gene") or args.protein
        info["name_for_tcmsp"] = info.get("recommended_name") or args.protein

        # 别名集合（用于 TCMSP 精确匹配重试）
        aliases = {args.protein, info.get("gene"), info.get("recommended_name"),
                   info.get("uniprot_id"), info.get("uniprot_ac")}
        aliases.update(info.get("gene_synonyms", []))
        aliases.update(info.get("alternative_names", []))
        info["name_aliases"] = [a for a in aliases if a]

    write_json(out / "01_protein.json", info)
    log(f"[01] AC={info.get('uniprot_ac')} gene={info.get('gene')} "
        f"tcmsp_name='{info.get('name_for_tcmsp')}' ambiguous={info.get('ambiguous')}")

    if info.get("ambiguous"):
        log(f"[01] 注意：命中 {len(info.get('candidates', []))} 个候选蛋白，"
            f"已默认选 annotationScore 最高者。如需其它请查看候选列表。", level="warning")
    return 0


if __name__ == "__main__":
    sys.exit(main())
