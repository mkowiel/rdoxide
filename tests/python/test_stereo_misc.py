"""Stereo helpers, atom helpers and enum exports (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdchem, rdmolops

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

WEDGED = """
  Mrv

  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8000    0.5000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8000    0.5000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.9000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.9000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
"""


def tags(mol):
    return [str(a.GetChiralTag()) for a in mol.GetAtoms()]


def test_needs_hs_and_atom_predicates():
    m = Chem.MolFromSmiles("c1ccccc1CC(=O)N")
    assert rdmolops.NeedsHs(m)
    assert not rdmolops.NeedsHs(Chem.AddHs(m))
    assert [rdmolops.CountAtomElec(a) for a in m.GetAtoms()] == [1, 1, 1, 1, 1, 1, -1, 1, 5, 2]
    assert [rdmolops.AtomHasConjugatedBond(a) for a in m.GetAtoms()] == [True] * 6 + [False] + [True] * 3


def test_free_atoms_raise_like_upstream():
    with pytest.raises(RuntimeError, match="calcImplicitValence"):
        rdmolops.CountAtomElec(Chem.Atom(6))
    with pytest.raises(RuntimeError, match="no owner"):
        rdmolops.AtomHasConjugatedBond(Chem.Atom(6))


def test_assign_chiral_types_from_bond_dirs():
    m = Chem.MolFromMolBlock(WEDGED)
    assert tags(m)[0] == "CHI_TETRAHEDRAL_CCW"
    for a in m.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rdmolops.AssignChiralTypesFromBondDirs(m)
    # a two-dimensional centre with a missing wedge direction stays unspecified
    assert set(tags(m)) == {"CHI_UNSPECIFIED"}
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        rdmolops.AssignChiralTypesFromBondDirs(m, 5)


def test_assign_chiral_tags_from_mol_parity():
    m = Chem.MolFromMolBlock(WEDGED)
    atom = m.GetAtomWithIdx(0)
    for parity, tag in ((1, "CHI_TETRAHEDRAL_CW"), (2, "CHI_TETRAHEDRAL_CCW"), (0, "CHI_UNSPECIFIED")):
        atom.SetIntProp("molParity", parity)
        rdmolops.AssignAtomChiralTagsFromMolParity(m)
        assert str(atom.GetChiralTag()) == tag
    atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    rdmolops.AssignAtomChiralTagsFromMolParity(m, replaceExistingTags=False)
    assert str(atom.GetChiralTag()) == "CHI_TETRAHEDRAL_CW"


def test_chirality_from_3d_structure():
    m = Chem.MolFromMolBlock(BLOCKS["chiral_chfcl"], removeHs=False)
    for a in m.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rdmolops.AssignAtomChiralTagsFromStructure(m)
    assert tags(m)[:3] == ["CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW", "CHI_UNSPECIFIED"]
    for a in m.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rdmolops.AssignStereochemistryFrom3D(m)
    assert Chem.MolToSmiles(m) == "[H]C([H])([H])[C@]([H])(F)Cl"
    # no conformer: nothing to do
    rdmolops.AssignStereochemistryFrom3D(Chem.MolFromSmiles("CC"))


def test_detect_bond_stereochemistry():
    m = Chem.MolFromMolBlock(BLOCKS["trans_butene"], removeHs=False)
    for b in m.GetBonds():
        b.SetStereo(Chem.BondStereo.STEREONONE)
    rdmolops.DetectBondStereochemistry(m)
    dirs = [str(b.GetBondDir()) for b in m.GetBonds()][:3]
    assert dirs == ["ENDUPRIGHT", "NONE", "ENDUPRIGHT"]
    other = Chem.Mol(m)
    rdmolops.DetectBondStereoChemistry(other, other.GetConformer())
    rdmolops.SetBondStereoFromDirections(other)
    assert [str(b.GetStereo()) for b in other.GetBonds()][:3] == ["STEREONONE", "STEREOTRANS", "STEREONONE"]
    with pytest.raises(RuntimeError, match="does not belong to molecule"):
        rdmolops.DetectBondStereoChemistry(other, m.GetConformer())
    with pytest.raises(RuntimeError, match="no owner"):
        rdmolops.DetectBondStereoChemistry(other, Chem.Conformer(3))


def test_set_terminal_atom_coords_errors():
    m = Chem.MolFromMolBlock(BLOCKS["chiral_chfcl"], removeHs=False)
    rdmolops.SetTerminalAtomCoords(m, 3, 1)
    with pytest.raises(RuntimeError, match="degenerate atoms"):
        rdmolops.SetTerminalAtomCoords(m, 1, 1)
    with pytest.raises(RuntimeError, match="Range Error"):
        rdmolops.SetTerminalAtomCoords(m, 100, 1)
    with pytest.raises(RuntimeError, match="bad atom degree"):
        rdmolops.SetTerminalAtomCoords(m, 0, 1)
    with pytest.raises(RuntimeError, match="no bond between atoms"):
        rdmolops.SetTerminalAtomCoords(m, 3, 5)
    with pytest.raises(RuntimeError, match="idx1"):
        rdmolops.SetTerminalAtomCoords(m, 3, 50)


def test_cleanup_functions():
    m = Chem.MolFromSmiles("C[C@H](C)C", sanitize=False)
    m.GetAtomWithIdx(1).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    rdmolops.CleanupChirality(m)
    assert str(m.GetAtomWithIdx(1).GetChiralTag()) == "CHI_UNSPECIFIED"
    rdmolops.CleanupAtropisomers(Chem.MolFromSmiles("CC"))
    organo = Chem.MolFromSmiles("C[Fe]", sanitize=False)
    rdmolops.CleanupOrganometallics(organo)
    assert Chem.MolToSmiles(organo) == "[CH3][Fe]"


def test_atom_labels_round_trip():
    m = Chem.MolFromSmiles("C=CC#N")
    atom = m.GetAtomWithIdx(0)
    assert rdchem.GetAtomRLabel(atom) == 0
    rdchem.SetAtomRLabel(atom, 5)
    assert rdchem.GetAtomRLabel(atom) == 5
    rdchem.SetAtomRLabel(atom, 0)
    assert not atom.HasProp("_MolFileRLabel")
    for bad in (100, -1):
        with pytest.raises(RuntimeError, match="rlabel out of range"):
            rdchem.SetAtomRLabel(atom, bad)
    for get, set_, prop in (
        (rdchem.GetAtomAlias, rdchem.SetAtomAlias, "molFileAlias"),
        (rdchem.GetAtomValue, rdchem.SetAtomValue, "molFileValue"),
    ):
        assert get(atom) == ""
        set_(atom, "foo")
        assert get(atom) == "foo" == atom.GetProp(prop)
        set_(atom, "")
        assert not atom.HasProp(prop)
    rdchem.SetSupplementalSmilesLabel(atom, "<xxx>")
    assert rdchem.GetSupplementalSmilesLabel(atom) == "<xxx>"
    assert Chem.MolToSmiles(m) == "C<xxx>=CC#N"
    rdchem.SetSupplementalSmilesLabel(atom, "")
    assert rdchem.GetSupplementalSmilesLabel(atom) == ""


def test_num_pi_electrons():
    assert [rdchem.GetNumPiElectrons(a) for a in Chem.MolFromSmiles("C=CC#N").GetAtoms()] == [1, 1, 2, 2]
    assert rdchem.GetNumPiElectrons(Chem.MolFromSmiles("c1ccccc1").GetAtomWithIdx(0)) == 1
    with pytest.raises(RuntimeError, match="no owner"):
        rdchem.GetNumPiElectrons(Chem.Atom(6))


def test_enum_members_are_exported():
    assert rdchem.CHI_TETRAHEDRAL_CW is Chem.ChiralType.CHI_TETRAHEDRAL_CW
    assert int(rdchem.CHI_ALLENE) == 5
    assert (int(rdchem.IMPLICIT), int(rdchem.EXPLICIT)) == (0, 1)
    assert rdchem.STEREO_AND is Chem.StereoGroupType.STEREO_AND
    assert [int(x) for x in (rdchem.COMPOSITE_AND, rdchem.COMPOSITE_OR, rdchem.COMPOSITE_XOR)] == [0, 1, 2]
    assert [int(rdchem.ALLOW_INCOMPLETE_OCTETS), int(rdchem.UNCONSTRAINED_ANIONS)] == [1, 16]
    assert rdchem.AllProps is Chem.PropertyPickleOptions.AllProps
    assert repr(rdchem.KEKULE_ALL) == "rdoxide.Chem.rdchem.ResonanceFlags.KEKULE_ALL"
