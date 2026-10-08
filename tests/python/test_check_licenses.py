"""The SPDX allow-list logic of tools/check_licenses.py."""

import importlib.util
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "check_licenses", Path(__file__).resolve().parents[2] / "tools" / "check_licenses.py")
check_licenses = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_licenses)


@pytest.mark.parametrize(
    "expression, allowed",
    [
        ("MIT", True),
        ("MIT OR Apache-2.0", True),
        ("Apache-2.0 WITH LLVM-exception", True),
        ("(MIT OR Apache-2.0) AND Unicode-3.0", True),
        ("MIT/Apache-2.0", True),
        ("GPL-3.0-only", False),
        ("GPL-3.0-only OR MIT", True),  # one allowed alternative is enough
        ("GPL-3.0-only AND MIT", False),  # every required license must be allowed
        ("MPL-2.0", False),
        ("Apache-2.0 WITH Some-Other-Exception", False),
    ],
)
def test_expression_evaluation(expression, allowed):
    ok, _ = check_licenses.evaluate(check_licenses.tokenize(expression))
    assert ok is allowed
