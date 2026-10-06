#!/usr/bin/env python3
"""
08_report.py — 汇总生成报告
输入: 01~07 全部 JSON
输出: report.md + report.html（交互式可视化）

报告内容:
  1. 靶点信息
  2. 候选化合物与对接结果表（结合能排序）
  3. 最优化合物 ADMET
  4. 相关临床研究
"""
import argparse
import json
import sys
import html as html_lib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import read_json, get_outdir, log, write_json


def esc(x):
    return html_lib.escape(str(x if x is not None else "—"))


def fmt_coords(v, nd=1):
    """坐标列表格式化：[75.45385...] → [75.5, 81.6, -8.7]"""
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(f"{round(float(x), nd):g}" for x in v) + "]"
    return str(v if v is not None else "—")


def load_all(outdir):
    out = Path(outdir)
    d = {}
    for key, fn in [("protein", "01_protein.json"), ("raw", "02_raw_compounds.json"),
                    ("clean", "03_clean_compounds.json"), ("receptor", "04_receptor.json"),
                    ("docking", "05_docking.json"), ("admet", "06_admet.json"),
                    ("clinical", "07_clinical.json")]:
        p = out / fn
        d[key] = read_json(p) if p.exists() else {}
    return d


def src_label(r):
    """来源展示：主库 + 交叉印证库（去前缀 ID，仅留库名去重）"""
    libs = []
    if r.get("source_db"):
        libs.append(r["source_db"])
    for s in r.get("cross_sources", []):
        lib = s.split(":")[0]
        if lib and lib not in libs:
            libs.append(lib)
    return " + ".join(libs) if libs else "—"


# ---------------- Markdown ----------------
def build_markdown(d):
    p = d.get("protein", {})
    dock = d.get("docking", {})
    admet = d.get("admet", {})
    clin = d.get("clinical", {})
    clean = d.get("clean", {})

    lines = []
    A = lines.append
    A("# 中药单体反向筛选与分子对接报告\n")
    A(f"> 靶点：**{p.get('gene') or p.get('user_input')}** ｜ "
      f"UniProt {p.get('uniprot_ac') or '—'} ｜ 生成时间自动\n")

    # 1. 靶点
    A("## 1. 靶点信息\n")
    A("| 项目 | 值 |")
    A("|------|----|")
    A(f"| 用户输入 | {p.get('user_input')} |")
    A(f"| 基因名 | {p.get('gene')} |")
    A(f"| UniProt AC | {p.get('uniprot_ac')} |")
    A(f"| 推荐名 | {p.get('recommended_name')} |")
    A(f"| TCMSP 匹配名 | {p.get('name_for_tcmsp')} |")
    rec = d.get("receptor", {})
    A(f"| 结构来源 | {rec.get('source')} |")
    A(f"| 对接盒子 | 中心 {fmt_coords(rec.get('box_center'))}，尺寸 {fmt_coords(rec.get('box_size'))}（{rec.get('box_mode')}）|")
    A("")

    # 2. 检索与清洗
    A("## 2. 检索与清洗\n")
    src = clean.get("sources", {})
    A(f"- HERB 检索到成分：{src.get('HERB', 0)} 条")
    A(f"- TCMSP 检索到成分：{src.get('TCMSP', 0)} 条")
    A(f"- 原始合计：{clean.get('total_raw', 0)} 条")
    A(f"- 清洗后候选：**{clean.get('total_final', 0)} 个**（去重、剔除非中药成分/内源物/离子溶剂等）")
    A("")

    # 3. 对接结果
    A("## 3. 候选化合物对接结果\n")
    results = dock.get("results", [])
    if results:
        A(f"对接引擎：AutoDock Vina ｜ exhaustiveness={dock.get('exhaustiveness')} ｜ "
          f"成功 {dock.get('n_docked')} 个\n")
        A("| 排名 | 化合物 | 来源 | MW | 结合能 (kcal/mol) |")
        A("|------|--------|------|-----|-------------------|")
        for r in results:
            A(f"| {r['rank']} | **{r['name']}** | {src_label(r)} | "
              f"{r.get('mw') or '—'} | {r['affinity']:.2f} |")
        A("")
        A(f"> 结合能越低表示预测结合越强。**最优化合物：{dock.get('best_compound')} "
          f"（{dock.get('best_affinity'):.2f} kcal/mol）**\n")
    else:
        A("_无对接结果_\n")

    # 4. ADMET
    A(f"## 4. 最优化合物 ADMET：{admet.get('compound', '—')}\n")
    pc = admet.get("pubchem") or {}
    rd = admet.get("rdkit") or {}
    if pc or rd:
        A("| 参数 | 值 | 说明 |")
        A("|------|----|------|")
        A(f"| 分子式 | {pc.get('MolecularFormula', '—')} | — |")
        A(f"| 分子量 MW | {rd.get('MW', pc.get('MolecularWeight', '—'))} | "
          f"{'符合 Lipinski' if (rd.get('MW') or 0) <= 500 else '超出 500'} |")
        A(f"| LogP | {rd.get('LogP', pc.get('XLogP', '—'))} | "
          f"{'合理' if 0 < (rd.get('LogP') or 0) < 5 else '需关注'} |")
        A(f"| TPSA | {rd.get('TPSA', pc.get('TPSA', '—'))} | "
          f"{'口服吸收良好' if (rd.get('TPSA') or 999) < 140 else '偏高'} |")
        A(f"| 氢键供体 HBD | {rd.get('HBD', pc.get('HBondDonorCount', '—'))} | |")
        A(f"| 氢键受体 HBA | {rd.get('HBA', pc.get('HBondAcceptorCount', '—'))} | |")
        A(f"| 可旋转键 | {rd.get('RotBonds', '—')} | |")
        A(f"| QED（类药性） | {rd.get('QED', '—')} | 越接近 1 越好 |")
        A(f"| Lipinski 违规 | {rd.get('LipinskiViolations', '—')} | "
          f"{'通过 (≤1)' if rd.get('LipinskiPass') else '不通过'} |")
        A(f"| Veber 规则 | {'通过' if rd.get('VeberPass') else '不通过'} | 口服可及性 |")
        A(f"| PAINS 警报 | {'是' if rd.get('PAINS') else '否'} | 泛筛选干扰化合物 |")
        A(f"| BBB 透过性 | {rd.get('BBB_hint', '—')} | 经验判断 |")
        A("")
        A("> ADMET 数据来自 PubChem 基础理化性质 + RDKit 本地计算的经验规则预测，**非实验值**，仅供初步参考。\n")
    else:
        A("_无 ADMET 数据_\n")

    # 5. 临床研究
    A(f"## 5. 临床研究（{clin.get('compound', '—')}）\n")
    trials = clin.get("trials", [])
    if trials:
        A(f"来源：ClinicalTrials.gov ｜ 共 {len(trials)} 项\n")
        A("| NCT ID | 标题 | 状态 | 期别 | 起始 |")
        A("|--------|------|------|------|------|")
        for t in trials:
            phase = "/".join(t.get("phases") or []) or "—"
            A(f"| [{t.get('nct_id')}]({t.get('url')}) | {t.get('title')} | "
              f"{t.get('status')} | {phase} | {t.get('start_date') or '—'} |")
        A("")
    else:
        A("_未检索到相关临床研究_\n")

    A("---")
    A("*本报告由 herb-target-screening skill 生成。对接打分与 ADMET 均为计算预测，"
      "不构成实验结论，请结合实验验证。*")
    return "\n".join(lines)


# ---------------- HTML ----------------
def build_html(d):
    p = d.get("protein", {})
    dock = d.get("docking", {})
    admet = d.get("admet", {})
    clin = d.get("clinical", {})
    clean = d.get("clean", {})
    rec = d.get("receptor", {})

    results = dock.get("results", [])
    names = [r["name"] for r in results]
    affinities = [round(r["affinity"], 2) for r in results]

    # 对接结果表行
    rows = "".join(
        f"<tr><td>{r['rank']}</td><td class='cname'>{esc(r['name'])}</td>"
        f"<td>{esc(src_label(r))}</td><td>{esc(r.get('mw'))}</td>"
        f"<td class='{'best' if r['rank']==1 else ''}'>{r['affinity']:.2f}</td></tr>"
        for r in results)

    # 临床表行
    clin_rows = "".join(
        f"<tr><td><a href='{esc(t.get('url'))}' target='_blank'>{esc(t.get('nct_id'))}</a></td>"
        f"<td>{esc(t.get('title'))}</td><td>{esc(t.get('status'))}</td>"
        f"<td>{esc('/'.join(t.get('phases') or []) or '—')}</td>"
        f"<td>{esc(t.get('start_date') or '—')}</td></tr>"
        for t in clin.get("trials", []))

    pc = admet.get("pubchem") or {}
    rd = admet.get("rdkit") or {}
    admet_rows = [
        ("分子式", esc(pc.get("MolecularFormula", "—")), ""),
        ("分子量 MW", rd.get("MW", pc.get("MolecularWeight", "—")),
         "符合 Lipinski" if (rd.get("MW") or 0) <= 500 else "超出 500"),
        ("LogP", rd.get("LogP", pc.get("XLogP", "—")),
         "合理" if 0 < (rd.get("LogP") or 0) < 5 else "需关注"),
        ("TPSA", rd.get("TPSA", pc.get("TPSA", "—")),
         "口服吸收良好" if (rd.get("TPSA") or 999) < 140 else "偏高"),
        ("HBD / HBA", f"{rd.get('HBD','—')} / {rd.get('HBA','—')}", ""),
        ("可旋转键", rd.get("RotBonds", "—"), ""),
        ("QED 类药性", rd.get("QED", "—"), "越接近 1 越好"),
        ("Lipinski 违规", rd.get("LipinskiViolations", "—"),
         "通过" if rd.get("LipinskiPass") else "不通过"),
        ("Veber 规则", "通过" if rd.get("VeberPass") else "不通过", "口服可及性"),
        ("PAINS 警报", "是" if rd.get("PAINS") else "否", "泛筛选干扰"),
        ("BBB 透过性", esc(rd.get("BBB_hint", "—")), "经验判断"),
    ]
    admet_html = "".join(
        f"<tr><td>{k}</td><td class='v'>{v}</td><td class='hint'>{h}</td></tr>"
        for k, v, h in admet_rows)

    src = clean.get("sources", {})
    best = dock.get("best_compound", "—")
    best_aff = dock.get("best_affinity")

    return f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>中药单体反向筛选 · {esc(p.get('gene'))}</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
:root{{--bg:#0f1420;--panel:#182131;--line:#2a3a52;--txt:#e6edf7;--sub:#8fa3bf;
--green:#4fd18b;--red:#ff5d5d;--amber:#ffb454;--blue:#5aa9ff}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;
background:var(--bg);color:var(--txt);padding:22px}}
.wrap{{max-width:1200px;margin:0 auto}}
h1{{font-size:23px;margin-bottom:6px}}
.sub{{color:var(--sub);font-size:13px;margin-bottom:20px}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:18px}}
@media(max-width:900px){{.grid{{grid-template-columns:1fr}}}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden;margin-bottom:18px}}
.card h2{{font-size:13px;color:var(--sub);text-transform:uppercase;letter-spacing:.8px;
padding:12px 16px;border-bottom:1px solid var(--line)}}
.card .body{{padding:14px 16px}}
table{{width:100%;border-collapse:collapse;font-size:13px}}
th,td{{padding:8px 10px;text-align:left;border-bottom:1px solid var(--line)}}
th{{color:var(--sub);font-weight:600;font-size:12px;text-transform:uppercase}}
tr:hover td{{background:rgba(90,169,255,.06)}}
td.best{{color:var(--green);font-weight:700}}
.cname{{font-weight:600}}
.v{{font-weight:600;color:var(--blue)}}
.hint{{color:var(--sub);font-size:12px}}
.kpi{{display:flex;gap:26px;flex-wrap:wrap}}
.kpi div{{text-align:center}}
.kpi .n{{font-size:26px;font-weight:700;color:var(--green)}}
.kpi .l{{font-size:12px;color:var(--sub)}}
a{{color:var(--blue);text-decoration:none}}
.banner{{background:linear-gradient(90deg,rgba(79,209,139,.14),transparent);
border-left:3px solid var(--green);padding:11px 15px;border-radius:6px;
font-size:13px;line-height:1.65;margin-bottom:18px}}
#chart{{width:100%;height:520px}}
.note{{font-size:12px;color:var(--sub);line-height:1.7;padding:10px 16px}}
</style></head><body><div class="wrap">
<h1>中药单体反向筛选与分子对接报告</h1>
<div class="sub">靶点 <b style="color:var(--txt)">{esc(p.get('gene'))}</b> ·
UniProt {esc(p.get('uniprot_ac'))} · 结构来源 {esc(rec.get('source'))}</div>

<div class="banner">
最优化合物：<b style="color:var(--green)">{esc(best)}</b>
{f"（{best_aff:.2f} kcal/mol）" if best_aff is not None else ""} ·
候选化合物 {clean.get('total_final', 0)} 个 ·
成功对接 {dock.get('n_docked', 0)} 个
</div>

<div class="grid">
<div class="card"><h2>检索与清洗概览</h2><div class="body"><div class="kpi">
<div><div class="n">{src.get('HERB',0)+src.get('TCMSP',0)}</div><div class="l">原始命中</div></div>
<div><div class="n">{clean.get('total_final',0)}</div><div class="l">清洗后候选</div></div>
<div><div class="n">{dock.get('n_docked',0)}</div><div class="l">成功对接</div></div>
<div><div class="n">{len(clin.get('trials',[]))}</div><div class="l">临床研究</div></div>
</div></div></div>

<div class="card"><h2>对接参数</h2><div class="body" style="font-size:13px;line-height:1.9">
引擎：AutoDock Vina 1.2.7<br>
exhaustiveness：{dock.get('exhaustiveness')}<br>
盒子中心：{esc(fmt_coords(rec.get('box_center')))}<br>
盒子尺寸：{esc(fmt_coords(rec.get('box_size')))}（{esc(rec.get('box_mode'))}）<br>
HERB 靶点 ID：{esc(d.get('raw',{}).get('HERB',{}).get('target_id'))}
</div></div>
</div>

<div class="card"><h2>候选化合物结合能排名</h2><div id="chart"></div>
<div class="note">横轴结合能（kcal/mol，越负越强），纵轴化合物。红色为最优化合物。</div></div>

<div class="card"><h2>对接结果明细</h2>
<table><thead><tr><th>排名</th><th>化合物</th><th>来源库</th><th>MW</th><th>结合能</th></tr></thead>
<tbody>{rows}</tbody></table></div>

<div class="card"><h2>最优化合物 ADMET — {esc(admet.get('compound'))}</h2>
<table><thead><tr><th>参数</th><th>值</th><th>说明</th></tr></thead>
<tbody>{admet_html}</tbody></table>
<div class="note">数据来自 PubChem 理化性质 + RDKit 规则预测，非实验值，仅供初步参考。</div></div>

<div class="card"><h2>相关临床研究 — {esc(clin.get('compound'))}（{len(clin.get('trials',[]))} 项）</h2>
<table><thead><tr><th>NCT ID</th><th>标题</th><th>状态</th><th>期别</th><th>起始</th></tr></thead>
<tbody>{clin_rows if clin_rows else "<tr><td colspan=5>未检索到</td></tr>"}</tbody></table></div>

<div class="note">本报告由 herb-target-screening skill 生成。对接打分与 ADMET 均为计算预测，不构成实验结论。</div>
</div>

<script>
var names={json.dumps(names, ensure_ascii=False)};
var aff={json.dumps(affinities)};
var colors=names.map((n,i)=>i===0?'#ff5d5d':'#5aa9ff');
var sorted=names.map((n,i)=>({{n:n,a:aff[i],c:colors[i]}})).sort((x,y)=>x.a-y.a);
Plotly.newPlot('chart',[{{
  x:sorted.map(d=>d.a), y:sorted.map(d=>d.n), type:'bar', orientation:'h',
  marker:{{color:sorted.map(d=>d.c)}},
  text:sorted.map(d=>d.a.toFixed(2)), textposition:'outside',
  hovertemplate:'%{{y}}<br>结合能 %{{x}} kcal/mol<extra></extra>'
}}],{{
  paper_bgcolor:'#182131', plot_bgcolor:'#182131',
  font:{{color:'#e6edf7',size:11}},
  margin:{{l:160,r:60,t:20,b:45}},
  xaxis:{{title:'结合能 (kcal/mol)', gridcolor:'#2a3a52', zerolinecolor:'#2a3a52'}},
  yaxis:{{automargin:true, gridcolor:'#2a3a52'}},
  height:Math.max(420, sorted.length*22)
}},{{responsive:true, displayModeBar:false}});
</script>
</body></html>"""


def main():
    ap = argparse.ArgumentParser(description="生成报告")
    ap.add_argument("--outdir", default="./tcm_run")
    args = ap.parse_args()
    out = get_outdir(args.outdir)

    d = load_all(out)
    md = build_markdown(d)
    (out / "report.md").write_text(md, encoding="utf-8")
    log(f"[08] 写出 report.md（{len(md)} 字节）")

    html = build_html(d)
    (out / "report.html").write_text(html, encoding="utf-8")
    log(f"[08] 写出 report.html（{len(html)} 字节）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
