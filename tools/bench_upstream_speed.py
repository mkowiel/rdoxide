#!/usr/bin/env python3
"""Time RDKit's own Python tests on original RDKit and on rdoxide.

Both interpreters must be the same Python version; the original RDKit lives in
its own virtual environment. Every upstream test file listed in
tests/upstream/files.txt runs in a fresh subprocess under each interpreter and
records the wall-clock time of each test. Only tests that pass on both sides
are compared.

Usage:
    python tools/bench_upstream_speed.py --rdkit-python PATH --out results.json [-j N] [-k SUBSTR]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from run_upstream_tests import UPSTREAM, test_files  # noqa: E402

RUNNER = ROOT / "tools" / "_timed_unittest_runner.py"


def run(python: str, path: Path, alias: bool, timeout: int) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "r.json"
        env = dict(os.environ, RDBASE=str(UPSTREAM), PYTHONDONTWRITEBYTECODE="1")
        env.pop("PYTHONPATH", None)
        if alias:
            env["PYTHONPATH"] = str(ROOT / "tools" / "_alias")
        cmd = [python, "-P", str(RUNNER), str(path), str(out)] + (["--alias"] if alias else [])
        try:
            subprocess.run(cmd, cwd=path.parent, env=env, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return {"tests": {}, "timeout": True, "file_seconds": float(timeout), "import_seconds": 0.0}
        return json.loads(out.read_text()) if out.exists() else {"tests": {}, "crash": True, "file_seconds": 0.0, "import_seconds": 0.0}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rdkit-python", required=True, help="python of the venv with the original RDKit")
    ap.add_argument("--out", required=True)
    ap.add_argument("-j", type=int, default=4)
    ap.add_argument("-k", default="")
    ap.add_argument("--timeout", type=int, default=600)
    args = ap.parse_args()

    files = [f for f in test_files() if args.k in str(f)]
    results: dict[str, dict] = {}

    def one(f: Path):
        rel = f.relative_to(UPSTREAM).as_posix()
        return rel, {
            "rdkit": run(args.rdkit_python, f, False, args.timeout),
            "rdoxide": run(sys.executable, f, True, args.timeout),
        }

    with cf.ThreadPoolExecutor(args.j) as ex:
        for rel, res in ex.map(one, files):
            results[rel] = res
    Path(args.out).write_text(json.dumps({"python": sys.version.split()[0], "files": results}))
    print(f"timed {len(files)} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
