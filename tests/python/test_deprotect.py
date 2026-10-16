from rdoxide import Chem
from rdoxide.Chem import rdDeprotect


def test_defaults_annotate_the_molecule():
    m = rdDeprotect.Deprotect(Chem.MolFromSmiles("N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C"))
    assert Chem.MolToSmiles(m) == Chem.CanonSmiles("NCc1ccccc1N")
    assert list(m.GetPropsAsDict()["DEPROTECTIONS"]) == ["Boc", "Boc"]
    assert m.GetIntProp("DEPROTECTION_COUNT") == 2


def test_in_place_and_custom_list():
    data = rdDeprotect.DeprotectData(
        "amine", "[N;H0,H1:1]C(=O)C>>[N:1]", "Ac", "acetyl")
    assert data.isValid() and data.abbreviation == "Ac"
    m = Chem.MolFromSmiles("CC(=O)NC")
    assert rdDeprotect.DeprotectInPlace(m, [data])
    assert Chem.MolToSmiles(m) == "CN"
    assert not rdDeprotect.DeprotectInPlace(m, [data])
    assert len(rdDeprotect.GetDeprotections()) == 25
