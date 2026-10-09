"""Bond wedging helpers (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdmolops

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)


def wedged(mol):
    """(begin, end, direction) of the bonds that carry a wedge."""
    return [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondDir()))
        for b in mol.GetBonds()
        if b.GetBondDir() != Chem.BondDir.NONE
    ]


def depicted(smiles):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    return mol


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H](F)Cl", [(1, 0, "BEGINWEDGE")]),
        ("N[C@@H](C)C(=O)O", [(1, 2, "BEGINDASH")]),
        ("C[C@H]1CC[C@@H](C)CC1", [(1, 0, "BEGINDASH"), (4, 5, "BEGINDASH")]),
        ("C[C@@]12CC[C@H](C1)C2", [(1, 0, "BEGINWEDGE"), (4, 3, "BEGINWEDGE")]),
        ("C1CC[C@@H]2CCCC[C@H]2C1", [(3, 2, "BEGINWEDGE"), (8, 7, "BEGINDASH")]),
    ],
)
def test_wedge_mol_bonds_picks_one_bond_per_center(smiles, expected):
    mol = depicted(smiles)
    rdmolops.WedgeMolBonds(mol, mol.GetConformer())
    assert wedged(mol) == expected


def test_wedge_two_bonds_if_possible():
    smiles = "F[C@](Cl)(Br)I"
    default = depicted(smiles)
    rdmolops.WedgeMolBonds(default, default.GetConformer())
    assert wedged(default) == [(1, 0, "BEGINDASH")]

    params = rdmolops.BondWedgingParameters()
    assert params.wedgeTwoBondsIfPossible is False
    params.wedgeTwoBondsIfPossible = True
    two = depicted(smiles)
    rdmolops.WedgeMolBonds(two, two.GetConformer(), params)
    assert wedged(two) == [(1, 0, "BEGINDASH"), (1, 4, "BEGINWEDGE")]

    # degree-three centres keep a single wedge
    three = depicted("C[C@H](F)Cl")
    rdmolops.WedgeMolBonds(three, three.GetConformer(), params)
    assert wedged(three) == [(1, 0, "BEGINWEDGE")]


def test_wedge_mol_bonds_checks_the_conformer():
    mol = depicted("C[C@H](F)Cl")
    other = depicted("C[C@@H](F)Cl")
    with pytest.raises(RuntimeError, match="conformer does not belong to molecule"):
        rdmolops.WedgeMolBonds(mol, other.GetConformer())
    with pytest.raises(RuntimeError, match="no owner"):
        rdmolops.WedgeMolBonds(mol, Chem.Conformer(4))


def test_wedge_bond():
    mol = depicted("C[C@H](F)Cl")
    rdmolops.WedgeBond(mol.GetBondWithIdx(0), 1, mol.GetConformer())
    assert wedged(mol) == [(0, 1, "BEGINWEDGE")]
    rdmolops.WedgeBond(mol.GetBondWithIdx(2), 1, mol.GetConformer())
    assert wedged(mol) == [(0, 1, "BEGINWEDGE"), (1, 3, "BEGINWEDGE")]

    other = depicted("C[C@H](F)Cl")
    with pytest.raises(RuntimeError, match="do not belong to same molecule"):
        rdmolops.WedgeBond(mol.GetBondWithIdx(0), 1, other.GetConformer())
    with pytest.raises(RuntimeError, match="Test Assert"):
        rdmolops.WedgeBond(mol.GetBondWithIdx(0), 0, mol.GetConformer())


def test_wedge_bond_ignores_non_single_bonds():
    mol = depicted("C/C=C/C")
    before = wedged(mol)
    rdmolops.WedgeBond(mol.GetBondWithIdx(1), 1, mol.GetConformer())
    assert wedged(mol) == before


def wedged_molblock():
    """A mol block whose first two bonds around the centre carry wedges."""
    mol = depicted("C[C@H](F)Cl")
    rdmolops.WedgeBond(mol.GetBondWithIdx(0), 1, mol.GetConformer())
    rdmolops.WedgeBond(mol.GetBondWithIdx(2), 1, mol.GetConformer())
    return Chem.MolToMolBlock(mol)


MOLBLOCK = wedged_molblock()


@pytest.mark.parametrize(
    "all_types, verify, expected",
    [
        (True, True, [(1, 3, "BEGINWEDGE")]),
        (True, False, [(0, 1, "BEGINWEDGE"), (1, 3, "BEGINWEDGE")]),
        (False, True, [(1, 3, "BEGINWEDGE")]),
        (False, False, [(0, 1, "BEGINWEDGE"), (1, 3, "BEGINWEDGE")]),
    ],
)
def test_reapply_mol_block_wedging(all_types, verify, expected):
    mol = Chem.MolFromMolBlock(MOLBLOCK)
    rdmolops.ReapplyMolBlockWedging(mol, all_types, verify)
    assert wedged(mol) == expected
    default = Chem.MolFromMolBlock(MOLBLOCK)
    rdmolops.ReapplyMolBlockWedging(default)
    assert wedged(default) == [(0, 1, "BEGINWEDGE"), (1, 3, "BEGINWEDGE")]


def test_reapply_mol_block_wedging_handles_either_double_bonds():
    block = """
  Mrv

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.8000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  3
  2  3  1  0
  1  4  1  1
M  END
"""
    mol = Chem.MolFromMolBlock(block)
    assert wedged(mol) == []
    rdmolops.ReapplyMolBlockWedging(mol)
    assert wedged(mol) == [(0, 1, "EITHERDOUBLE"), (0, 3, "BEGINWEDGE")]
    rdmolops.ReapplyMolBlockWedging(mol, True, True)
    assert wedged(mol) == [(0, 1, "UNKNOWN")]
    assert mol.GetBondWithIdx(0).GetStereo() == Chem.BondStereo.STEREOANY


def test_remove_non_explicit_3d_chirality():
    mol = Chem.MolFromMolBlock(BLOCKS["chiral"], removeHs=False)
    tags = [str(a.GetChiralTag()) for a in mol.GetAtoms()][:3]
    assert tags == ["CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CCW", "CHI_UNSPECIFIED"]
    flagged = Chem.MolFromMolBlock(BLOCKS["chiral_chfcl"], removeHs=False)
    for a in flagged.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rdmolops.AssignAtomChiralTagsFromStructure(flagged)
    assert [a.HasProp("_NonExplicit3DChirality") for a in flagged.GetAtoms()][:3] == [True, True, False]
    rdmolops.RemoveNonExplicit3DChirality(flagged)
    assert [str(a.GetChiralTag()) for a in flagged.GetAtoms()][:3] == ["CHI_UNSPECIFIED"] * 3
    assert not any(a.HasProp("_NonExplicit3DChirality") for a in flagged.GetAtoms())
