#!/usr/bin/env python3
"""Regenerate docs/bench/upstream_speed.html with every section.

Runs all the benchmarks (upstream test timings, the slowest files again alone, substructure matching,
multicore scaling, fingerprints and bit vectors) on both libraries and feeds them to tools/bench_report.py,
so that a regeneration cannot silently drop a section.

    python tools/regen_speed_report.py --rdkit-python /path/to/rdkit-venv/bin/python

Run it with the interpreter that has rdoxide installed as a release build
(`maturin develop --release`); both interpreters must be the same Python version, and both need pandas
and pillow for the upstream tests that use them. Intermediate JSON goes to --workdir (default .work/speed).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "tools"
SOLO_FILES = 5


def run(cmd: list[str], out: Path | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    res = subprocess.run(cmd, stdout=subprocess.PIPE if out else None, text=True, cwd=ROOT)
    if res.returncode != 0:
        sys.exit(f"command failed ({res.returncode}): {' '.join(cmd)}")
    if out:
        out.write_text(res.stdout)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rdkit-python", required=True, help="interpreter with the original RDKit wheel")
    ap.add_argument("--rdoxide-python", default=sys.executable, help="interpreter with rdoxide (default: this one)")
    ap.add_argument("--workdir", default=str(ROOT / ".work" / "speed"))
    ap.add_argument("--out", default=str(ROOT / "docs" / "bench" / "upstream_speed.html"))
    ap.add_argument("--rdoxide-version", default="", help="label for the page (default: git short hash)")
    ap.add_argument("-j", type=int, default=4, help="files run in parallel for the main upstream timing")
    args = ap.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    rd, ox = args.rdkit_python, args.rdoxide_python

    rdkit_version = subprocess.run(
        [rd, "-c", "import rdkit; print(rdkit.__version__)"], capture_output=True, text=True, check=True
    ).stdout.strip()
    version = args.rdoxide_version or subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=ROOT
    ).stdout.strip()

    # 1. every upstream test file, both libraries (parallel files add noise to short tests)
    main_json = work / "speed.json"
    run([ox, str(TOOLS / "bench_upstream_speed.py"), "--rdkit-python", rd, "--out", str(main_json), "-j", str(args.j)])

    # 2. the files with the largest runtimes again, one at a time
    data = json.loads(main_json.read_text())
    ranked = []
    for f, r in data["files"].items():
        a = sum(v[1] for v in r["rdkit"]["tests"].values() if v[0] == "passed")
        b = sum(v[1] for v in r["rdoxide"]["tests"].values() if v[0] == "passed")
        if a > 0.05 and b > 0.05:
            ranked.append((a + b, f))
    solo = []
    for i, (_, f) in enumerate(sorted(ranked, reverse=True)[:SOLO_FILES], 1):
        path = work / f"solo{i}.json"
        run([ox, str(TOOLS / "bench_upstream_speed.py"), "--rdkit-python", rd, "--out", str(path), "-j", "1", "-k", f])
        solo.append(str(path))

    # 3. micro benchmarks, each on its own library
    for name in ("substruct", "threads", "micro", "bitvect"):
        run([rd, str(TOOLS / f"bench_{name}.py"), "rdkit"], work / f"{name}_rdkit.json")
        run([ox, str(TOOLS / f"bench_{name}.py"), "rdoxide"], work / f"{name}_rdoxide.json")

    pair = lambda n: [str(work / f"{n}_rdkit.json"), str(work / f"{n}_rdoxide.json")]  # noqa: E731
    cmd = [
        ox, str(TOOLS / "bench_report.py"), str(main_json), args.out,
        "--rdkit-version", rdkit_version, "--rdoxide-version", version,
        "--solo", *solo,
        "--substruct", *pair("substruct"),
        "--threads", *pair("threads"),
        "--fingerprints", *pair("micro"), *pair("bitvect"),
    ]  # fmt: skip
    run(cmd)
    page = Path(args.out).read_text()
    missing = [h for h in ("slowest files re-run alone", "Substructure micro-benchmark", "Multicore scaling",
                           "Fingerprint calculation", "Tests where rdoxide gains most") if h not in page]
    if missing:
        sys.exit(f"the regenerated page lacks sections: {missing}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
