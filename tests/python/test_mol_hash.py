"""rdMolHash bindings (the full set of hashes is checked by the upstream testMolHash.py)."""
import pytest

from rdoxide import Chem
from rdoxide.Chem import rdMolHash

HF = rdMolHash.HashFunction


def test_enum_and_defaults():
    assert int(HF.AnonymousGraph) == 1 and int(HF.HetAtomProtomerv2) == 19
    assert rdMolHash.excludeFromTautomerismProp == "_MolHashExcludeFromTautomerism"
    m = Chem.MolFromSmiles("OC(=O)C")
    assert rdMolHash.MolHash(m, HF.MolFormula) == "C2H4O2"
    # the molecule is not modified
    assert Chem.MolToSmiles(m) == "CC(=O)O"


def test_excluded_atoms_stay_out_of_the_tautomer_system():
    m = Chem.MolFromSmiles("CC(=O)CC=O")
    plain = rdMolHash.MolHash(m, HF.HetAtomTautomerv2)
    for a in m.GetAtoms():
        a.SetBoolProp(rdMolHash.excludeFromTautomerismProp, True)
    assert rdMolHash.MolHash(m, HF.HetAtomTautomerv2) != plain


def test_cxsmiles():
    m = Chem.MolFromSmiles("C[C@@H](O)[C@@H](C)C |o1:1,3|")
    assert rdMolHash.MolHash(m, HF.CanonicalSmiles, True).endswith("|o1:1,3|") or "|" in rdMolHash.MolHash(m, HF.CanonicalSmiles, True)
    assert "|" not in rdMolHash.MolHash(m, HF.CanonicalSmiles, False)


def test_bad_function():
    with pytest.raises(ValueError):
        rdMolHash.MolHash(Chem.MolFromSmiles("CC"), 99)
