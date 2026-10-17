"""Helper for bench_upstream_speed.py: run one upstream unittest file and write
per-test outcomes and wall-clock seconds as JSON.

Usage: _timed_unittest_runner.py <test file> <out.json> [--alias]
`--alias` makes `import rdkit` resolve to rdoxide; without it the real RDKit
is used.
"""

import importlib.util
import json
import multiprocessing
import sys
import time
import unittest
from pathlib import Path


def main(path: str, out: str, alias: bool) -> None:
    t_file = time.perf_counter()
    if alias:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import rdkit_alias

        rdkit_alias.install()
    try:
        multiprocessing.set_start_method("fork", force=True)
    except (ValueError, RuntimeError):
        pass
    results: dict[str, list] = {}
    import_seconds = 0.0
    try:
        spec = importlib.util.spec_from_file_location("__upstream_test__", path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["__upstream_test__"] = mod
        t0 = time.perf_counter()
        spec.loader.exec_module(mod)
        import_seconds = time.perf_counter() - t0
    except BaseException as e:  # noqa: BLE001
        results["<collection>"] = ["failed", 0.0]
    else:
        suite = unittest.defaultTestLoader.loadTestsFromModule(mod)

        class Timed(unittest.TestResult):
            def _id(self, test):
                return f"{type(test).__name__}::{getattr(test, '_testMethodName', str(test))}"

            def startTest(self, test):
                super().startTest(test)
                self._t = time.perf_counter()

            def _put(self, test, outcome):
                results[self._id(test)] = [outcome, time.perf_counter() - self._t]

            def addSuccess(self, test):
                self._put(test, "passed")

            def addFailure(self, test, err):
                self._put(test, "failed")

            def addError(self, test, err):
                self._put(test, "failed")

            def addSkip(self, test, reason):
                self._put(test, "skipped")

            def addExpectedFailure(self, test, err):
                self._put(test, "passed")

            def addUnexpectedSuccess(self, test):
                self._put(test, "failed")

        suite.run(Timed())
    with open(out, "w") as fh:
        json.dump(
            {
                "tests": results,
                "import_seconds": import_seconds,
                "file_seconds": time.perf_counter() - t_file,
            },
            fh,
        )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], "--alias" in sys.argv[3:])
