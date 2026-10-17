"""ExplicitBitVect atom properties and the matching query."""

import pytest
from rdoxide import Chem, DataStructs
from rdoxide.Chem import rdqueries


def test_atom_bit_vect_round_trip():
    m = Chem.MolFromSmiles("CC")
    for atom in m.GetAtoms():
        bv = DataStructs.ExplicitBitVect(100)
        bv.SetBit(atom.GetIdx())
        atom.SetExplicitBitVectProp("prop", bv)
    for atom in m.GetAtoms():
        assert atom.GetExplicitBitVectProp("prop").GetBit(atom.GetIdx())
    with pytest.raises(KeyError):
        m.GetAtomWithIdx(0).GetExplicitBitVectProp("missing")
    m.GetAtomWithIdx(0).SetIntProp("n", 1)
    with pytest.raises(ValueError):
        m.GetAtomWithIdx(0).GetExplicitBitVectProp("n")


def test_bit_vect_query_uses_tanimoto_tolerance():
    m = Chem.MolFromSmiles("CC")
    for atom in m.GetAtoms():
        bv = DataStructs.ExplicitBitVect(4)
        bv.SetBit(atom.GetIdx())
        atom.SetExplicitBitVectProp("prop", bv)

    def matches(bit, tolerance):
        query = Chem.MolFromSmarts("C")
        probe = DataStructs.ExplicitBitVect(4)
        probe.SetBit(bit)
        for atom in query.GetAtoms():
            atom.ExpandQuery(rdqueries.HasBitVectPropWithValueQueryAtom("prop", probe, tolerance=tolerance))
        return m.GetSubstructMatches(query)

    assert matches(1, 0.0) == ((1,),)
    assert matches(0, 0.0) == ((0,),)
    assert matches(0, 1.0) == ((0,), (1,))


def test_atoms_matching_query():
    bv = DataStructs.ExplicitBitVect(4)
    bv.SetBit(0)
    bv.SetBit(2)
    qa = rdqueries.HasBitVectPropWithValueQueryAtom("prop", bv, tolerance=0.0)
    m = Chem.MolFromSmiles("CC")
    m.GetAtomWithIdx(0).SetExplicitBitVectProp("prop", bv)
    assert tuple(a.GetIdx() for a in m.GetAtomsMatchingQuery(qa)) == (0,)
