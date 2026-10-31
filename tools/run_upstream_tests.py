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
EXPECTED_LONG = SUITE / "expected_pass_long.txt"
RUNNER = Path(__file__).resolve().parent / "_unittest_file_runner.py"


LONG = False  # set by --long: also run the tests upstream skips as "long"


def read_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")}


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
            # Contrib tests start `python` subprocesses
            PATH=str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", ""),
            **({"RDOXIDE_UPSTREAM_LONG": "1"} if LONG else {}),
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


NOT_RUN_REASONS = {
    "Code/Demos/": "C++/Boost.Python demos",
    "Contrib/mmpa/": "stale expectation: the test expects `c1ccc(cc1)[*:1]` but RDKit 2026.09.1 writes `c1ccc([*:1])cc1` (the official wheel fails it too)",
    "Contrib/ConformerParser/": "needs the separate C++ module `rdConformerParser` that Contrib builds",
    "Contrib/": "Contrib script",
    "Data/": "data-set scripts",
    "Docs/": "documentation scripts",
    "Projects/": "database command-line tool project",
    "Scripts/": "RDKit's own test driver",
    "External/ChemDraw/": "needs the proprietary ChemDraw SDK (binary CDX)",
    "Code/GraphMol/MolDraw2D/test_dir/": "drawing regression driver for the C++ tests",
    "rdkit/sping/": "legacy PostScript drawing tests (`sping`)",
    "rdkit/TestRunner.py": "RDKit's own test runner",
}


def write_report(path: Path, results, files) -> None:
    """Markdown list of failing, skipped and not-run upstream tests."""
    gaps = {}
    for ln in (SUITE / "known_gaps.txt").read_text().splitlines():
        if ln.strip() and not ln.startswith("#"):
            tid, _, why = ln.partition(" | ")
            gaps[tid.strip()] = why.strip()
    failed = sorted(t for t, (o, _) in results.items() if o == "failed")
    skipped = sorted((t, d) for t, (o, d) in results.items() if o == "skipped")
    ran = {str(f.relative_to(UPSTREAM)) for f in files}
    import re

    every = {
        p.relative_to(UPSTREAM).as_posix()
        for p in UPSTREAM.rglob("*.py")
        if re.match(r"(test|Test|UnitTest).*\.py$", p.name) or p.name.endswith(("_test.py", "Test.py", "_tests.py"))
    }
    not_run = sorted(every - ran)
    total = len(results)
    passed = sum(1 for o, _ in results.values() if o == "passed")
    out = [
        "# Upstream tests that do not pass",
        "",
        "Generated by `tools/run_upstream_tests.py --long --report docs/UPSTREAM_TESTS.md` (`--long` also runs",
        "upstream's long regression tests); reasons for the",
        "failures are kept in `tests/upstream/known_gaps.txt`.",
        "",
        f"{passed} of {total} collected tests pass in {len(files)} files; {len(failed)} fail, {len(skipped)} are skipped.",
        "",
        f"## Failing ({len(failed)})",
        "",
        "| Test | Why |",
        "|------|-----|",
    ]
    for t in failed:
        out.append(f"| `{t}` | {gaps.get(t, '**unclassified**')} |")
    out += ["", f"## Skipped by the tests themselves ({len(skipped)})", "", "| Test | Reason |", "|------|--------|"]
    for t, d in skipped:
        reason = (d.strip().splitlines() or [""])[-1].strip() or "skipped"
        out.append(f"| `{t}` | {reason} |")
    out += ["", f"## Test files that are not run ({len(not_run)})", "", "| File | Why |", "|------|-----|"]
    for f in not_run:
        why = next((r for pre, r in NOT_RUN_REASONS.items() if f.startswith(pre)), "")
        if not why:
            why = {
                "Code/GraphMol/Substruct/UnitTestSubstruct.py": "runs the C++ test executable `test2`",
                "Code/GraphMol/UnitTestQueryMol.py": "runs the C++ test executable `querytest.exe`",
                "rdkit/ML/Cluster/murtagh_test.py": "script of the Murtagh clustering demo (not a unittest)",
            }.get(f, "")
        out.append(f"| `{f}` | {why or 'not in `tests/upstream/files.txt`'} |")
    out.append("")
    path.write_text("\n".join(out))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true", help="add newly passing tests to expected_pass.txt")
    ap.add_argument("--long", action="store_true", help="enable upstream's long regression tests (`-l`); their ratchet is expected_pass_long.txt")
    ap.add_argument("-j", type=int, default=os.cpu_count() or 2)
    ap.add_argument("-k", default="", help="only run files whose path contains this string")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--report", metavar="FILE", help="write a markdown list of failing, skipped and unrun tests")
    ap.add_argument("-v", "--verbose", action="store_true", help="show failures and skips")
    args = ap.parse_args()
    global LONG
    LONG = args.long

    if not UPSTREAM.exists():
        print("upstream/rdkit missing - run tools/sync_upstream.py fetch", file=sys.stderr)
        return 2
    files = [f for f in test_files() if args.k in str(f)]
    results: dict[str, tuple[str, str]] = {}
    with cf.ThreadPoolExecutor(args.j) as ex:
        for res in ex.map(lambda f: run_file_locked(f, args.timeout), files):
            results.update(res)

    expected_normal = read_ids(EXPECTED)
    expected = expected_normal | (read_ids(EXPECTED_LONG) if LONG else set())
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
    if args.report:
        write_report(Path(args.report), results, files)
    if new:
        print(f"{len(new)} newly passing test(s)" + ("" if args.update else " (run with --update to record)"))
        for t in new:
            print(f"  NEW: {t}")
    if args.update:
        keep = {t for t in expected if t.split("::", 1)[0] not in ran_files}
        header = (
            "# Upstream RDKit Python tests known to pass on rdoxide (ratchet).\n"
            "# Maintained by tools/run_upstream_tests.py --update\n"
        )
        if LONG:
            long_header = (
                "# Long upstream tests (enabled with `run_upstream_tests.py --long`) known to pass.\n"
                "# Maintained by tools/run_upstream_tests.py --long --update\n"
            )
            EXPECTED_LONG.write_text(long_header + "".join(f"{t}\n" for t in sorted((keep | passed) - expected_normal)))
        else:
            EXPECTED.write_text(header + "".join(f"{t}\n" for t in sorted(keep | passed)))
    return 1 if regressions else 0


if __name__ == "__main__":
    sys.exit(main())
