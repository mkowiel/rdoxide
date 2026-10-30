"""Helper for run_upstream_tests.py: run one unittest file as a standalone
module (not as part of upstream's `rdkit` package) and print JSON results."""

import importlib.util
import json
import multiprocessing
import os
import sys
import unittest
from pathlib import Path


def run_pytest_functions(path: str, results: dict) -> None:
    """Files with pytest-style `test_*` functions (no unittest cases)."""
    try:
        import pytest
    except ImportError:
        return

    class Plugin:
        def pytest_runtest_logreport(self, report):
            tail = report.nodeid.split("::", 1)[-1]
            key = f"pytest::{tail}"
            if report.when == "call":
                if report.passed:
                    results[key] = ["passed", ""]
                elif report.skipped:
                    results[key] = ["skipped", str(report.longrepr)[-200:]]
                else:
                    results[key] = ["failed", str(report.longrepr)[-2000:]]
            elif report.failed:
                results[key] = ["failed", str(report.longrepr)[-2000:]]
            elif report.skipped:
                results[key] = ["skipped", ""]

    pytest.main(["-q", "-p", "no:cacheprovider", "--rootdir", str(Path(path).parent), path],
                plugins=[Plugin()])


def main(path: str, out: str) -> None:
    # upstream's tests import `rdkit`; rdoxide's package is `rdoxide`
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import rdkit_alias

    rdkit_alias.install()
    results = {}
    if "/Contrib/" in path:
        # Contrib scripts import their neighbours (`import sascorer`); appended, so nothing is shadowed
        sys.path.append(str(Path(path).resolve().parent))
    # The test module is loaded under a private name, so child processes that
    # target functions in it can only work if they are forked, not spawned.
    try:
        multiprocessing.set_start_method("fork", force=True)
    except (ValueError, RuntimeError):
        pass
    try:
        spec = importlib.util.spec_from_file_location("__upstream_test__", path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["__upstream_test__"] = mod
        spec.loader.exec_module(mod)
        if os.environ.get("RDOXIDE_UPSTREAM_LONG") and hasattr(mod, "doLong"):
            mod.doLong = True  # upstream's `-l` switch for the long regression tests
    except BaseException as e:  # noqa: BLE001 - import failures are results too
        results["<collection>"] = ["failed", f"{type(e).__name__}: {e}"]
    else:
        suite = unittest.defaultTestLoader.loadTestsFromModule(mod)

        class Collector(unittest.TestResult):
            def _id(self, test):
                return f"{type(test).__name__}::{getattr(test, '_testMethodName', str(test))}"

            def addSuccess(self, test):
                results[self._id(test)] = ["passed", ""]

            def addFailure(self, test, err):
                results[self._id(test)] = ["failed", self._exc_info_to_string(err, test)[-2000:]]

            def addError(self, test, err):
                results[self._id(test)] = ["failed", self._exc_info_to_string(err, test)[-2000:]]

            def addSkip(self, test, reason):
                results[self._id(test)] = ["skipped", reason]

            def addExpectedFailure(self, test, err):
                results[self._id(test)] = ["passed", "expected failure"]

            def addUnexpectedSuccess(self, test):
                results[self._id(test)] = ["failed", "unexpected success"]

        if suite.countTestCases():
            suite.run(Collector())
        else:
            run_pytest_functions(path, results)
            if not results:
                # a script without test cases: importing it without error is the test
                results["<script>"] = ["passed", ""]
    with open(out, "w") as fh:
        json.dump(results, fh)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
