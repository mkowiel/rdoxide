"""SubstanceGroups in binary pickles; expectations validated against RDKit 2026.09.1."""

import pickle

from rdkit import Chem
from rdkit.Geometry import Point3D


def sup_mol():
    mol = Chem.MolFromSmiles("CCOC")
    group = Chem.CreateMolSubstanceGroup(mol, "SUP")
    for a in (1, 2):
        group.AddAtomWithIdx(a)
    for b in (0, 2):
        group.AddBondWithIdx(b)
    group.AddAttachPoint(1, 0, "1")
    group.AddAttachPoint(2, -1, "")
    group.AddBracket([Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(2, 2, 2)])
    group.AddCState(0, Point3D(1.5, 2.5, 3.5))
    group.SetProp("LABEL", "x")
    group.SetUnsignedProp("COMPNO", 3)
    return mol


def test_sup_group_round_trips_through_binary():
    mol = sup_mol()
    for copy in (Chem.Mol(mol.ToBinary()), pickle.loads(pickle.dumps(mol))):
        (group,) = Chem.GetMolSubstanceGroups(copy)
        assert group.GetProp("TYPE") == "SUP"
        assert group.GetProp("LABEL") == "x"
        assert group.GetUnsignedProp("COMPNO") == 3
        assert list(group.GetAtoms()) == [1, 2]
        assert list(group.GetBonds()) == [0, 2]
        assert [(a.aIdx, a.lvIdx, a.id) for a in group.GetAttachPoints()] == [(1, 0, "1"), (2, -1, "")]
        (state,) = group.GetCStates()
        assert state.bondIdx == 0
        assert tuple(state.vector) == (1.5, 2.5, 3.5)
        (bracket,) = group.GetBrackets()
        assert [tuple(p) for p in bracket] == [(0, 0, 0), (1, 1, 1), (2, 2, 2)]
        assert copy.ToBinary() == mol.ToBinary()


def test_non_sup_cstates_carry_no_vector():
    mol = Chem.MolFromSmiles("CCOC")
    group = Chem.CreateMolSubstanceGroup(mol, "SRU")
    group.AddAtomWithIdx(1)
    group.AddAtomWithIdx(2)
    group.AddBondWithIdx(0)
    group.AddBondWithIdx(2)
    group.AddCState(2, Point3D(1, 2, 3))
    copy = Chem.Mol(mol.ToBinary())
    (state,) = Chem.GetMolSubstanceGroups(copy)[0].GetCStates()
    assert state.bondIdx == 2
    assert tuple(state.vector) == (0, 0, 0)


def test_wide_indices_for_large_molecules():
    mol = Chem.MolFromSmiles("C" * 300)
    group = Chem.CreateMolSubstanceGroup(mol, "GEN")
    for i in range(0, 300, 3):
        group.AddAtomWithIdx(i)
    (copy_group,) = Chem.GetMolSubstanceGroups(Chem.Mol(mol.ToBinary()))
    assert list(copy_group.GetAtoms()) == list(range(0, 300, 3))
