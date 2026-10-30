#!/usr/bin/env python3
"""Run RDKit's own Python test-suite against rdoxide (a "ratchet" test).

Upstream test files are executed *in place* inside the upstream checkout
(`upstream/rdkit`, see `tools/sync_upstream.py fetch`) with `RDBASE` pointing
there, so the tests find their data files exactly as in RDKit's CTest setup.
The `rdkit` package that gets imported is the installed rdoxide one.

Results are compared with `tests/upstream/expected_pass.txt`:
  * a test listed there that now fails  -> regression, exit code 1
  * a test that passes but is not listed -> reported; `--update` adds it

Usage:
    python tools/run_upstream_tests.py [--update] [-j N] [-k PATTERN]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import glob
import os
import subprocess
import threading
import sys
import tempfile
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UPSTREAM = ROOT / "upstream" / "rdkit"
SUITE = ROOT / "tests" / "upstream"
FILES = SUITE / "files.txt"
EXPECTED = SUITE / "expected_pass.txt"
RUNNER = Path(__file__).resolve().parent / "_unittest_file_runner.py"


def test_files() -> list[Path]:
    out: list[Path] = []
    for line in FILES.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.extend(sorted(Path(p) for p in glob.glob(str(UPSTREAM / line), recursive=True)))
    return out


def run_file(path: Path, timeout: int) -> dict[str, tuple[str, str]]:
    """Run one upstream test file; return {test_id: (outcome, details)}."""
    rel = path.relative_to(UPSTREAM).as_posix()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "results.json"
        env = dict(
            os.environ,
            RDBASE=str(UPSTREAM),
            PYTHONDONTWRITEBYTECODE="1",
            # `import rdkit` resolves to rdoxide (also in subprocesses of tests)
            PYTHONPATH=str(ROOT / "tools" / "_alias"),
        )
        # -P: don't put the (upstream) cwd on sys.path, so `import rdoxide`
        # resolves to the installed rdoxide package
        cmd = [sys.executable, "-P", str(RUNNER), str(path), str(out)]
        try:
            proc = subprocess.run(cmd, cwd=path.parent, env=env, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return {f"{rel}::<timeout>": ("failed", "timeout")}
        if not out.exists():
            return {f"{rel}::<crash>": ("failed", proc.stderr.decode(errors="replace")[-2000:])}
        data = json.loads(out.read_text())
        return {f"{rel}::{k}": tuple(v) for k, v in data.items()} or {f"{rel}::<no tests>": ("failed", "")}


# these two files write the same output files (test_data/outNCI_few.sdf, ...), so they must not
# run at the same time
_SHARED_OUTPUT_LOCK = threading.Lock()
_SHARED_OUTPUT_FILES = ("Code/GraphMol/Wrap/rough_test.py", "Code/GraphMol/nbWrap/rough_test.py")


def run_file_locked(path: Path, timeout: int) -> dict[str, tuple[str, str]]:
    if str(path).endswith(_SHARED_OUTPUT_FILES):
        with _SHARED_OUTPUT_LOCK:
            return run_file(path, timeout)
    return run_file(path, timeout)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true", help="add newly passing tests to expected_pass.txt")
    ap.add_argument("-j", type=int, default=os.cpu_count() or 2)
    ap.add_argument("-k", default="", help="only run files whose path contains this string")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("-v", "--verbose", action="store_true", help="show failures and skips")
    args = ap.parse_args()

    if not UPSTREAM.exists():
        print("upstream/rdkit missing - run tools/sync_upstream.py fetch", file=sys.stderr)
        return 2
    files = [f for f in test_files() if args.k in str(f)]
    results: dict[str, tuple[str, str]] = {}
    with cf.ThreadPoolExecutor(args.j) as ex:
        for res in ex.map(lambda f: run_file_locked(f, args.timeout), files):
            results.update(res)

    expected = set()
    if EXPECTED.exists():
        expected = {ln.strip() for ln in EXPECTED.read_text().splitlines() if ln.strip() and not ln.startswith("#")}
    passed = {t for t, (o, _) in results.items() if o == "passed"}
    ran_files = {t.split("::", 1)[0] for t in results}
    regressions = sorted(t for t in expected if t.split("::", 1)[0] in ran_files and t not in passed)
    new = sorted(passed - expected)

    total = len(results)
    print(f"upstream tests: {len(passed)}/{total} passed in {len(files)} files")
    for t in regressions:
        print(f"REGRESSION: {t}\n    {results.get(t, ('', 'not run'))[1][-500:]}")
    if args.verbose:
        for t, (o, d) in sorted(results.items()):
            if o in ("failed", "skipped"):
                print(f"{o.upper()} {t}\n    {d.strip().splitlines()[-1] if d.strip() else ''}")
    if new:
        print(f"{len(new)} newly passing test(s)" + ("" if args.update else " (run with --update to record)"))
        for t in new:
            print(f"  NEW: {t}")
    if args.update:
        keep = {t for t in expected if t.split("::", 1)[0] not in ran_files}
        EXPECTED.write_text(
            "# Upstream RDKit Python tests known to pass on rdoxide (ratchet).\n"
            "# Maintained by tools/run_upstream_tests.py --update\n"
            + "".join(f"{t}\n" for t in sorted(keep | passed))
        )
    return 1 if regressions else 0


if __name__ == "__main__":
    sys.exit(main())
