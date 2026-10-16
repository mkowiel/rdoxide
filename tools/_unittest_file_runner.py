"""Helper for run_upstream_tests.py: run one unittest file as a standalone
module (not as part of upstream's `rdkit` package) and print JSON results."""

import importlib.util
import json
import multiprocessing
import sys
import unittest


def main(path: str, out: str) -> None:
    results = {}
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

        suite.run(Collector())
    with open(out, "w") as fh:
        json.dump(results, fh)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
