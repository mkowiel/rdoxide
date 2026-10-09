"""GetMolFrags, CombineMols, RenumberAtoms and friends (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Geometry import Point3D

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

SALT = "CC(=O)[O-].[NH3+]C"


def test_get_mol_frags():
    m = Chem.MolFromSmiles(SALT)
    assert rdmolops.GetMolFrags(m) == ((0, 1, 2, 3), (4, 5))
    frags, mapping = [], []
    mols = rdmolops.GetMolFrags(m, asMols=True, frags=frags, fragsMolAtomMapping=mapping)
    assert [Chem.MolToSmiles(x) for x in mols] == ["CC(=O)[O-]", "C[NH3+]"]
    assert frags == [0, 0, 0, 0, 1, 1]
    assert mapping == [(0, 1, 2, 3), (4, 5)]


def test_formal_charge_and_shortest_path():
    assert rdmolops.GetFormalCharge(Chem.MolFromSmiles("C[NH3+]")) == 1
    assert rdmolops.GetFormalCharge(Chem.MolFromSmiles(SALT)) == 0
    assert rdmolops.GetShortestPath(Chem.MolFromSmiles("C1CC1CCO"), 0, 5) == (0, 2, 3, 4, 5)
    with pytest.raises(ValueError, match="bad atom index"):
        rdmolops.GetShortestPath(Chem.MolFromSmiles(SALT), 0, 99)
    with pytest.raises(RuntimeError, match="aid1 != aid2"):
        rdmolops.GetShortestPath(Chem.MolFromSmiles(SALT), 0, 0)


def test_renumber_atoms_keeps_stereo():
    m = Chem.MolFromSmiles("O[C@H](C)C[C@@H](N)C")
    r = rdmolops.RenumberAtoms(m, [6, 5, 4, 3, 2, 1, 0])
    assert Chem.MolToSmiles(r) == Chem.MolToSmiles(m)
    assert [a.GetSymbol() for a in r.GetAtoms()] == ["C", "N", "C", "C", "C", "C", "O"]
    with pytest.raises(ValueError, match="atomCounts shorter than the number of atoms"):
        rdmolops.RenumberAtoms(m, [0, 1])


def test_renumber_atoms_ring_stereo_props():
    m = Chem.MolFromSmiles("[C@@H]1(F)CC[C@H](Cl)CC1")
    r = rdmolops.RenumberAtoms(m, [2, 0, 1, 3, 4, 5, 6, 7])
    assert Chem.MolToSmiles(r) == Chem.MolToSmiles(m)
    assert r.GetRingInfo().NumRings() == 1


def test_renumber_atoms_moves_conformers():
    m = Chem.MolFromMolBlock(BLOCKS["caffeine"], removeHs=False)
    order = list(range(m.GetNumAtoms()))[::-1]
    r = rdmolops.RenumberAtoms(m, order)
    n = m.GetNumAtoms()
    for new, old in enumerate(order):
        assert tuple(r.GetConformer().GetAtomPosition(new)) == tuple(m.GetConformer().GetAtomPosition(old))
    assert r.GetNumAtoms() == n and r.GetConformer().Is3D()


def test_remove_stereochemistry():
    m = Chem.MolFromSmiles("F/C=C/Cl")
    rdmolops.RemoveStereochemistry(m)
    assert Chem.MolToSmiles(m) == "FC=CCl"
    m = Chem.MolFromSmiles("C[C@H](F)Cl")
    rdmolops.RemoveStereochemistry(m)
    assert Chem.MolToSmiles(m) == "CC(F)Cl"


def test_combine_mols_offsets_conformers():
    e = Chem.MolFromMolBlock(BLOCKS["ethanol"], removeHs=False)
    c = Chem.MolFromMolBlock(BLOCKS["caffeine"], removeHs=False)
    x = rdmolops.CombineMols(e, c, Point3D(1.5, -2, 3))
    assert x.GetNumAtoms() == e.GetNumAtoms() + c.GetNumAtoms()
    assert x.GetNumConformers() == 1
    first = c.GetConformer().GetAtomPosition(0)
    got = x.GetConformer().GetAtomPosition(e.GetNumAtoms())
    assert (got.x, got.y, got.z) == pytest.approx((first.x + 1.5, first.y - 2, first.z + 3))
    plain = rdmolops.CombineMols(Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("O"))
    assert Chem.MolToSmiles(plain) == "CC.O"


def test_combine_mols_bad_offset():
    with pytest.raises(TypeError, match="did not match C\\+\\+ signature"):
        rdmolops.CombineMols(Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("O"), (1, 2, 3))
