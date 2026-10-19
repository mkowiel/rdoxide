from rdoxide import Chem


def test_cdx_support_is_reported_absent():
    assert Chem.HasChemDrawCDXSupport() is False
