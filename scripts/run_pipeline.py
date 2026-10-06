#!/usr/bin/env python3
"""
run_pipeline.py — 一键编排完整流程
用法: python run_pipeline.py --protein MMP1 --outdir ./tcm_run

依次执行 01~08，任一步失败给出提示但不强行中断（除关键步骤）。
"""
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent

STEPS = [
    ("01_standardize_protein.py", "蛋白名标准化", True),
    ("02_fetch_compounds.py", "多库检索中药单体", True),
    ("03_clean_dedupe.py", "清洗去重", True),
    ("04_prepare_receptor.py", "受体准备", True),
    ("05_batch_docking.py", "批量分子对接", False),  # 失败可继续（报告仍可出）
    ("06_admet.py", "ADMET 查询", False),
    ("07_clinical_trials.py", "临床研究检索", False),
    ("08_report.py", "报告生成", False),
]


def run(step, outdir, extra):
    script = HERE / step
    cmd = [sys.executable, str(script), "--outdir", outdir] + extra
    print(f"\n{'='*60}\n▶ {step}\n{'='*60}", flush=True)
    r = subprocess.run(cmd)
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description="中药单体反向筛选全流程")
    ap.add_argument("--protein", required=True, help="蛋白质名称/基因名")
    ap.add_argument("--outdir", default="./tcm_run")
    ap.add_argument("--top", type=int, default=0, help="精简到 N 个（0=全部对接）")
    ap.add_argument("--pdb", help="手动指定 PDB ID")
    ap.add_argument("--workers", type=int, default=0, help="对接并行进程数")
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    extra_01 = ["--protein", args.protein]
    extra_03 = ["--top", str(args.top)] if args.top else []
    extra_04 = ["--pdb", args.pdb] if args.pdb else []
    extra_05 = ["--workers", str(args.workers)] if args.workers else []

    extras = {"01": extra_01, "03": extra_03, "04": extra_04, "05": extra_05}

    for step, desc, critical in STEPS:
        key = step[:2]
        rc = run(step, args.outdir, extras.get(key, []))
        if rc != 0:
            mark = "关键步骤" if critical else "非关键步骤"
            print(f"⚠ {desc} 返回码 {rc}（{mark}）", flush=True)
            if critical:
                print(f"✗ 关键步骤失败，流程终止于 {step}")
                return rc

    print(f"\n{'='*60}")
    print(f"✓ 全部完成！结果目录: {args.outdir}")
    print(f"  报告: {Path(args.outdir)/'report.html'}")
    print(f"        {Path(args.outdir)/'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
