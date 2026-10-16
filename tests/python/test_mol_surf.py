"""Labute ASA, TPSA and the VSA binnings (expectations from RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolDescriptors as D

# smiles: (Labute ASA, Labute ASA without implicit Hs, TPSA, TPSA with S and P)
CASES = {
    "CC(=O)Nc1ccc(O)cc1": (64.66690538413799, 63.94513703834635, 49.33, 49.33),
    "OC(=O)c1ccccc1N": (58.091991783146064, 57.2903955680267, 63.32, 63.32),
    "CS(=O)(=O)Nc1ccccc1": (66.1069833099167, 65.70431909632079, 46.17, 54.55),
    "C[N+](C)(C)C": (33.80388252781807, 32.7764874823146, 0.0, 0.0),
    "N#CC=O": (23.658294923738637, 22.46440175078022, 40.86, 40.86),
    "C1CN1": (19.64822417230497, 18.461173388146065, 21.94, 21.94),
    "CSC": (25.2273031675914, 24.417413147866565, 0.0, 25.3),
    "OP(=O)(O)O": (28.30662084449607, 27.11376439711333, 77.76, 87.57),
    "CCOC(=O)c1ccc(N)cc1": (71.14118748925014, 70.49924701278623, 52.32, 52.32),
    "c1ccncc1": (36.651051100443475, 35.70348392426772, 12.89, 12.89),
}


@pytest.mark.parametrize("smiles", CASES)
def test_labute_and_tpsa(smiles):
    labute, labute_no_hs, tpsa, tpsa_sp = CASES[smiles]
    mol = lambda: Chem.MolFromSmiles(smiles)  # noqa: E731 - results are cached on the molecule
    assert D.CalcLabuteASA(mol()) == pytest.approx(labute, rel=1e-12)
    assert D.CalcLabuteASA(mol(), includeHs=False) == pytest.approx(labute_no_hs, rel=1e-12)
    assert D.CalcTPSA(mol()) == pytest.approx(tpsa, abs=1e-9)
    assert D.CalcTPSA(mol(), includeSandP=True) == pytest.approx(tpsa_sp, abs=1e-9)


def test_contributions():
    m = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    contribs = D._CalcTPSAContribs(m)
    assert isinstance(contribs, tuple)
    assert contribs == (0.0, 0.0, 17.07, 12.03, 0.0, 0.0, 0.0, 0.0, 20.23, 0.0, 0.0)
    atoms, h = D._CalcLabuteASAContribs(Chem.MolFromSmiles("CCO"))
    assert len(atoms) == 3
    assert h == pytest.approx(1.2612803353779227, rel=1e-12)


def test_results_are_cached_on_the_molecule():
    m = Chem.MolFromSmiles("CCO")
    first = D.CalcLabuteASA(m, includeHs=True)
    # like upstream, the cached value is returned whatever includeHs says
    assert D.CalcLabuteASA(m, includeHs=False) == first
    assert D.CalcLabuteASA(m, includeHs=False, force=True) == pytest.approx(18.682762484374372, rel=1e-12)
    assert m.HasProp("_labuteASA")


def test_vsa_binnings():
    m = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    slogp = D.SlogP_VSA_(m)
    assert len(slogp) == 12
    assert slogp[:5] == pytest.approx(
        [5.316788604006331, 11.013707124192212, 4.794537184071822, 0.0, 6.923737199690624], rel=1e-12)
    assert len(D.SMR_VSA_(m)) == 10
    assert D.SlogP_VSA_(m, [0.0, 0.5]) == pytest.approx(
        [21.125032912270363, 36.876591748220626, 5.749511833283905], rel=1e-12)
    # an empty sequence selects the default bins
    assert D.SlogP_VSA_(m, ()) == slogp


def test_custom_prop_vsa():
    m = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1S")
    for atom in m.GetAtoms():
        atom.SetDoubleProp("p", atom.GetIdx() / 10.0)
    assert D.CustomProp_VSA_(m, "p", [0.1, 0.5]) == pytest.approx(
        [6.923737199690624, 21.70589179211322, 46.579413057040206], rel=1e-12)
    # atoms without the property count as 1
    assert D.CustomProp_VSA_(m, "missing", [0.5, 1.5]) == pytest.approx(
        [0.0, 75.20904204884405, 0.0], rel=1e-12)


def test_custom_prop_vsa_property_types():
    ints = Chem.MolFromSmiles("CCO")
    for atom in ints.GetAtoms():
        atom.SetIntProp("q", 1)
    with pytest.raises(RuntimeError, match="bad any cast"):
        D.CustomProp_VSA_(ints, "q", [0.5, 1.5])
    strings = Chem.MolFromSmiles("CCO")
    for atom in strings.GetAtoms():
        atom.SetProp("q", "2.5")
    assert D.CustomProp_VSA_(strings, "q", [0.5, 1.5])[2] > 0
    with pytest.raises(TypeError):
        D.CustomProp_VSA_(strings, "q", 3)


def test_versions():
    assert D._CalcLabuteASA_version == "1.0.2"
    assert D._CalcTPSA_version == "2.0.0"
