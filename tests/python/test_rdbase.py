"""rdoxide.rdBase basics (upstream logging behaviour is covered by the verbatim
upstream rdkit/UnitTestLogging.py, run via tools/run_upstream_tests.py)."""

import pytest

import rdoxide
from rdoxide import RDLogger, rdBase


def test_versions():
    assert rdBase.rdkitVersion == "2026.09.1"
    assert rdoxide.__version__ == rdBase.rdkitVersion
    assert rdBase._version() == rdBase.rdkitVersion


def test_module_names():
    assert rdBase.__name__ == "rdoxide.rdBase"
    import rdoxide.rdBase as again

    assert again is rdBase


def test_block_logs():
    rdBase.EnableLog("rdApp.error")
    with rdBase.BlockLogs():
        assert "rdApp.error:disabled" in rdBase.LogStatus()
    assert "rdApp.error:enabled" in rdBase.LogStatus()


def test_rdlogger_levels():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    assert "rdApp.error:disabled" in rdBase.LogStatus()
    lg.setLevel(RDLogger.WARNING)
    assert "rdApp.warning:enabled" in rdBase.LogStatus()


def _nested(prefix, elem):
    # the class names are C++-mangled upstream and differ between standard libraries
    names = [n for n in dir(rdBase) if n.startswith(prefix) and elem in n]
    assert len(names) == 1, names
    return getattr(rdBase, names[0])


def test_nested_vector_types():
    for cls, tag in ((rdBase._vecti, "vectorIi"), (rdBase._vectj, "vectorIj"), (rdBase._vectd, "vectorId")):
        v = _nested("_vect", tag)()
        assert len(v) == 0
        inner = cls()
        inner.append(1)
        v.append(inner)
        v.append(cls())
        assert len(v) == 2
        assert list(v[0]) == [1]
        assert [len(x) for x in v] == [1, 0]
        assert len(v[0:1]) == 1
        del v[1]
        assert len(v) == 1


def test_list_vector_types():
    for tag in ("vectorIi", "vectorIj"):
        lst = _nested("_list", tag)()
        assert len(lst) == 0
        with pytest.raises((TypeError, AttributeError)):
            lst[0] = []
