#!/usr/bin/env python3
"""Turn the JSON from bench_upstream_speed.py into a short HTML summary."""

from __future__ import annotations

import argparse
import html
import json
import statistics
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("out")
    ap.add_argument("--rdkit-version", default="")
    ap.add_argument("--rdoxide-version", default="")
    ap.add_argument("--solo", nargs="*", default=[], help="JSON files of files re-run alone (no parallel noise)")
    ap.add_argument("--substruct", nargs=2, metavar=("RDKIT_JSON", "RDOXIDE_JSON"), help="tools/bench_substruct.py outputs")
    ap.add_argument("--threads", nargs=2, metavar=("RDKIT_JSON", "RDOXIDE_JSON"), help="tools/bench_threads.py outputs")
    args = ap.parse_args()
    data = json.loads(Path(args.results).read_text())
    rows = []  # (test id, rdkit s, rdoxide s)
    file_rows = []  # (file, rdkit s, rdoxide s)
    only_rdkit = 0
    n_total = {"rdkit": 0, "rdoxide": 0}
    for f, r in data["files"].items():
        a, b = r["rdkit"]["tests"], r["rdoxide"]["tests"]
        for side, t in (("rdkit", a), ("rdoxide", b)):
            n_total[side] += sum(1 for v in t.values() if v[0] == "passed")
        common = [k for k in a if a[k][0] == "passed" and k in b and b[k][0] == "passed"]
        only_rdkit += sum(1 for k in a if a[k][0] == "passed" and not (k in b and b[k][0] == "passed"))
        for k in common:
            rows.append((f"{f}::{k}", a[k][1], b[k][1]))
        if common:
            file_rows.append(
                (f, sum(a[k][1] for k in common), sum(b[k][1] for k in common), r["rdkit"]["file_seconds"], r["rdoxide"]["file_seconds"])
            )
    t_rd = sum(r[1] for r in rows)
    t_ox = sum(r[2] for r in rows)
    # file wall time (imports included) only for files where both sides pass some test
    common_files = {f for f, *_ in file_rows}
    wall_rd = sum(data["files"][f]["rdkit"]["file_seconds"] for f in common_files)
    wall_ox = sum(data["files"][f]["rdoxide"]["file_seconds"] for f in common_files)
    # a slowdown only counts when it is above timer noise: >20 ms and >20 %
    slower = sorted(
        (r for r in rows if r[2] - r[1] > 0.02 and r[2] > 1.2 * r[1]), key=lambda r: r[2] - r[1], reverse=True
    )
    solo_rows = []
    for sp in args.solo:
        for f, r in json.loads(Path(sp).read_text())["files"].items():
            a, b = r["rdkit"]["tests"], r["rdoxide"]["tests"]
            c = [k for k in a if a[k][0] == "passed" and k in b and b[k][0] == "passed"]
            solo_rows.append((f.rsplit("/", 1)[-1], len(c), sum(a[k][1] for k in c), sum(b[k][1] for k in c)))
    slower_ratio = sorted((r for r in rows if r[2] > r[1] and r[2] > 0.05), key=lambda r: r[2] / max(r[1], 1e-6), reverse=True)
    faster = sorted(rows, key=lambda r: r[1] - r[2], reverse=True)
    ratios = [r[2] / r[1] for r in rows if r[1] > 0.001]
    median = statistics.median(ratios) if ratios else float("nan")

    def tr(r, extra=""):
        ratio = r[2] / r[1] if r[1] > 0 else float("inf")
        return f"<tr><td>{html.escape(r[0])}</td><td>{r[1]:.3f}</td><td>{r[2]:.3f}</td><td>{ratio:.1f}x</td></tr>"

    def table(title, items, n=10):
        body = "".join(tr(r) for r in items[:n]) or "<tr><td colspan=4>none</td></tr>"
        return (
            f"<h2>{title}</h2><table><tr><th>test</th><th>RDKit s</th><th>rdoxide s</th><th>rdoxide / RDKit</th></tr>{body}</table>"
        )

    solo_html = ""
    if solo_rows:
        body = "".join(
            f"<tr><td>{html.escape(n)} ({c} tests)</td><td>{a:.2f}</td><td>{b:.2f}</td><td>{b / a:.2f}x</td></tr>"
            for n, c, a, b in solo_rows
        )
        solo_html = f"<h2>The slowest files re-run alone (no parallel noise)</h2><table><tr><th>file</th><th>RDKit s</th><th>rdoxide s</th><th>rdoxide / RDKit</th></tr>{body}</table>"

    sub_html = thr_html = ""
    if args.substruct:
        a, b = (json.loads(Path(x).read_text()) for x in args.substruct)
        body = ""
        for q, ra in a.items():
            if q == "molecules":
                continue
            rb = b[q]
            body += (
                f"<tr><td>{html.escape(q)}</td><td>{ra['hits']}</td>"
                f"<td>{ra['has']:.2f}</td><td>{rb['has']:.2f}</td><td>{rb['has'] / ra['has']:.2f}x</td>"
                f"<td>{ra['matches']:.2f}</td><td>{rb['matches']:.2f}</td><td>{rb['matches'] / ra['matches']:.2f}x</td></tr>"
            )
        sub_html = (
            f"<h2>Substructure micro-benchmark ({a['molecules']} molecules, ms per 1000 checks)</h2>"
            "<table><tr><th>query</th><th>hits</th><th>Has RDKit</th><th>Has rdoxide</th><th>ratio</th>"
            f"<th>Matches RDKit</th><th>Matches rdoxide</th><th>ratio</th></tr>{body}</table>"
        )
    if args.threads:
        a, b = (json.loads(Path(x).read_text()) for x in args.threads)
        body = ""
        for op, ra in a.items():
            if op == "cores":
                continue
            rb = b[op]
            body += (
                f"<tr><td>{html.escape(op)}</td><td>{ra['1']:.0f}</td><td>{rb['1']:.0f}</td>"
                f"<td>{rb['all']:.0f}</td><td>{ra['1'] / rb['all']:.1f}x</td></tr>"
            )
        thr_html = (
            f"<h2>Multicore scaling ({b['cores']} cores, ms; RDKit runs these single-threaded)</h2>"
            "<table><tr><th>operation</th><th>RDKit</th><th>rdoxide, 1 thread</th><th>rdoxide, all cores</th>"
            f"<th>RDKit / rdoxide all</th></tr>{body}</table>"
        )

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>rdoxide vs RDKit speed</title>
<style>
:root{{--bg:#fff;--fg:#1c1c1e;--mut:#6b6b70;--line:#e3e3e6;--acc:#0b6bcb;--card:#f6f7f9}}
@media (prefers-color-scheme:dark){{:root{{--bg:#16171a;--fg:#ececef;--mut:#9b9ba3;--line:#2c2d33;--acc:#5aa9ff;--card:#1f2024}}}}
body{{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}}
main{{max-width:900px;margin:0 auto;padding:24px 16px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:16px;margin:28px 0 8px}}
p.sub{{color:var(--mut);margin:0 0 20px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px}}
.card b{{display:block;font-size:22px;color:var(--acc)}} .card span{{color:var(--mut);font-size:13px}}
table{{width:100%;border-collapse:collapse;font-size:13px}}
th,td{{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);word-break:break-word}}
td:nth-child(n+2),th:nth-child(n+2){{text-align:right;white-space:nowrap}}
.note{{color:var(--mut);font-size:13px}}
</style></head><body><main>
<h1>rdoxide vs original RDKit: upstream test-suite speed</h1>
<p class="sub">Same Python {html.escape(data['python'])}, same machine. RDKit {html.escape(args.rdkit_version)} (PyPI wheel) vs rdoxide {html.escape(args.rdoxide_version)} (release build). Only tests that pass on both are compared.</p>
<div class="cards">
<div class="card"><b>{t_rd:.1f} s → {t_ox:.1f} s</b><span>summed test time, {len(rows)} common passing tests</span></div>
<div class="card"><b>{(t_rd/t_ox if t_ox else float('inf')):.2f}x</b><span>rdoxide speed-up on that sum (RDKit time / rdoxide time)</span></div>
<div class="card"><b>{median:.2f}x</b><span>median per-test time ratio, rdoxide / RDKit (below 1 = rdoxide faster)</span></div>
<div class="card"><b>{wall_rd:.0f} s → {wall_ox:.0f} s</b><span>file wall time incl. imports, {len(common_files)} files with common passing tests</span></div>
<div class="card"><b>{len(slower)} / {len(rows)}</b><span>tests clearly slower on rdoxide (&gt;20 ms and &gt;20 %)</span></div>
<div class="card"><b>{len(rows)} / {n_total['rdoxide'] - len(rows)} / {only_rdkit}</b><span>tests passing on both / only on rdoxide / only on RDKit ({n_total['rdoxide']} pass on rdoxide, {n_total['rdkit']} on RDKit)</span></div>
</div>
{table('Tests that run slower on rdoxide (largest absolute loss)', slower)}
{solo_html}
{sub_html}
{thr_html}
{table('Tests where rdoxide gains most', faster)}
<p class="note">Method: tools/bench_upstream_speed.py runs each upstream test file in its own process under both interpreters (4 files in parallel), timing every test; tools/bench_report.py builds this page. Parallel runs add noise to short tests; the timings of the slowest ones are worth re-running alone.</p>
</main></body></html>"""
    Path(args.out).write_text(page)
    print(f"{len(rows)} common tests; RDKit {t_rd:.1f}s, rdoxide {t_ox:.1f}s; slower: {len(slower)}")


if __name__ == "__main__":
    main()
