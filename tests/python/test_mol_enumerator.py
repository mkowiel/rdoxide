import pytest
from rdoxide import Chem
from rdoxide.Chem import rdMolEnumerator


def test_sru_enumeration_with_defaults():
    mol = Chem.MolFromSmiles("FCN(CC)-* |Sg:n:2,5:2-3:ht|")
    bundle = rdMolEnumerator.Enumerate(mol)
    assert len(bundle) == 2 and bundle.Size() == 2
    assert [Chem.MolToSmiles(m) for m in bundle] == ["*N(CC)N(*)CF", "*N(CC)N(*)N(*)CF"]
    assert bundle.GetMol(1).GetNumAtoms() == bundle[1].GetNumAtoms()
    with pytest.raises(IndexError):
        bundle[2]
    assert rdMolEnumerator.Enumerate(mol, 1).Size() == 1


def test_params_and_operations():
    mol = Chem.MolFromSmiles("*-CO-* |$star_e;;;star_e$,Sg:n:2,1::ht|")
    ps = rdMolEnumerator.MolEnumeratorParams(rdMolEnumerator.EnumeratorType.RepeatUnit)
    assert ps.maxToEnumerate == 1000
    ps.maxToEnumerate = 2
    with pytest.raises(AttributeError):
        ps.notAnAttribute = 1
    assert rdMolEnumerator.Enumerate(mol, ps).Size() == 2
    ps = rdMolEnumerator.MolEnumeratorParams()
    ps.SetEnumerationOperator(rdMolEnumerator.EnumeratorType.LinkNode)
    assert rdMolEnumerator.Enumerate(mol, ps).Size() == 0
    with pytest.raises(ValueError):
        rdMolEnumerator.Enumerate(mol, rdMolEnumerator.MolEnumeratorParams())


def test_mol_bundle():
    bundle = Chem.MolBundle()
    assert bundle.AddMol(Chem.MolFromSmiles("CC")) == 1
    assert bundle.AddMol(Chem.MolFromSmiles("CCC")) == 2
    assert [m.GetNumAtoms() for m in bundle] == [2, 3]
