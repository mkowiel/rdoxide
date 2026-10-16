"""rdoxide.rdBase basics (upstream logging behaviour is covered by the verbatim
upstream rdkit/UnitTestLogging.py, run via tools/run_upstream_tests.py)."""

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
