"""Enhanced stereo helpers (expectations from RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdmolops

ABS, OR, AND = (
    Chem.StereoGroupType.STEREO_ABSOLUTE,
    Chem.StereoGroupType.STEREO_OR,
    Chem.StereoGroupType.STEREO_AND,
)


def with_groups(smiles, *groups, assign=True):
    """A molecule with stereo groups given as (type, atom indices, read id)."""
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    if assign:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mol.SetStereoGroups(
        [Chem.CreateStereoGroup(kind, mol, atomIds=list(atoms), readId=read_id) for kind, atoms, read_id in groups]
    )
    return mol


def atom_notes(mol):
    return [a.GetProp("atomNote") if a.HasProp("atomNote") else None for a in mol.GetAtoms()]


def bond_notes(mol):
    return [b.GetProp("bondNote") if b.HasProp("bondNote") else None for b in mol.GetBonds()]


def test_add_stereo_annotations_uses_groups_and_cip_codes():
    mol = with_groups("C[C@H](F)[C@@H](Cl)Br", (ABS, [1], 0), (OR, [3], 1))
    rdmolops.AddStereoAnnotations(mol)
    assert atom_notes(mol) == [None, "abs (S)", None, "or1", None, None]

    mol = with_groups("C[C@H](F)[C@@H](Cl)Br", (OR, [1], 1), (AND, [3], 2))
    rdmolops.AddStereoAnnotations(mol)
    assert atom_notes(mol) == [None, "or1", None, "and1", None, None]


def test_add_stereo_annotations_labels_ungrouped_centers_and_double_bonds():
    mol = with_groups("F/C=C/[C@H](Cl)C", (ABS, [3], 0))
    rdmolops.AddStereoAnnotations(mol)
    assert atom_notes(mol)[3] == "abs (R)"
    assert bond_notes(mol) == [None, "(E)", None, None, None]

    plain = with_groups("C[C@H](O)[C@@H](C)O")
    rdmolops.AddStereoAnnotations(plain)
    assert atom_notes(plain) == [None, "(S)", None, "(R)", None, None]


def test_add_stereo_annotations_custom_labels():
    mol = with_groups("C[C@H](F)[C@@H](Cl)Br", (ABS, [1], 0), (OR, [3], 1))
    rdmolops.AddStereoAnnotations(mol, "A{id}{cip}", "O{id}", "", "", "B{cip}")
    assert atom_notes(mol) == [None, "A0S", None, "O1", None, None]
    # an empty label switches that kind of annotation off
    mol = with_groups("F/C=C/[C@H](Cl)C", (ABS, [3], 0))
    rdmolops.AddStereoAnnotations(mol, "A{id}{cip}", "O{id}", "", "", "B{cip}")
    assert (atom_notes(mol)[3], bond_notes(mol)[1]) == ("A0R", "BE")


def test_simplify_enhanced_stereo():
    mol = with_groups("C[C@H](F)[C@@H](Cl)Br", (AND, [1, 3], 1))
    rdmolops.SimplifyEnhancedStereo(mol)
    assert mol.GetProp("molNote") == "AND enantiomer"
    assert len(mol.GetStereoGroups()) == 0
    assert not any(a.HasProp("_CIPCode") for a in mol.GetAtoms())

    mol = with_groups("C[C@H](F)Cl", (OR, [1], 1))
    rdmolops.SimplifyEnhancedStereo(mol, removeAffectedStereoGroups=False)
    assert mol.GetProp("molNote") == "OR enantiomer"
    assert len(mol.GetStereoGroups()) == 1


@pytest.mark.parametrize(
    "groups",
    [
        [(ABS, [1], 0), (OR, [3], 1)],  # several groups
        [(ABS, [1], 0)],  # an ABS group is not an enantiomer label
        [],
    ],
)
def test_simplify_enhanced_stereo_leaves_other_molecules_alone(groups):
    smiles = "C[C@H](F)[C@@H](Cl)Br" if len(groups) != 1 else "C[C@H](F)Cl"
    mol = with_groups(smiles, *groups)
    count = len(mol.GetStereoGroups())
    rdmolops.SimplifyEnhancedStereo(mol)
    assert not mol.HasProp("molNote")
    assert len(mol.GetStereoGroups()) == count


def test_find_meso_centers():
    mol = Chem.MolFromSmiles("C[C@H](O)[C@@H](C)O")
    assert rdmolops.FindMesoCenters(mol) == ((1, 3),)
    assert [a.GetIntProp("_mesoOtherAtom") for a in mol.GetAtoms() if a.HasProp("_mesoOtherAtom")] == [3, 1]
    assert rdmolops.FindMesoCenters(mol, includeIsotopes=False, includeAtomMaps=True) == ((1, 3),)
    assert rdmolops.FindMesoCenters(Chem.MolFromSmiles("C[C@H](O)[C@H](C)O")) == ()
    assert rdmolops.FindMesoCenters(Chem.MolFromSmiles("CC")) == ()
    # centers on the same ring are meso even though their ranks differ
    ring = Chem.MolFromSmiles("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O")
    assert rdmolops.FindMesoCenters(ring)[0] == (1, 2)
    assert len(rdmolops.FindMesoCenters(ring)) == 10


def test_translate_chiral_flag_to_stereo_groups():
    def groups(mol):
        return [(str(g.GetGroupType()), [a.GetIdx() for a in g.GetAtoms()]) for g in mol.GetStereoGroups()]

    def flagged(flag, *existing):
        mol = with_groups("C[C@H](F)[C@@H](Cl)Br", *existing, assign=False)
        if flag is not None:
            mol.SetIntProp("_MolFileChiralFlag", flag)
        return mol

    mol = flagged(1)
    rdmolops.TranslateChiralFlagToStereoGroups(mol)
    assert groups(mol) == [("STEREO_ABSOLUTE", [1, 3])]
    assert not mol.HasProp("_MolFileChiralFlag")

    mol = flagged(0)
    rdmolops.TranslateChiralFlagToStereoGroups(mol)
    assert groups(mol) == [("STEREO_AND", [1, 3])]
    mol = flagged(0)
    rdmolops.TranslateChiralFlagToStereoGroups(mol, Chem.StereoGroupType.STEREO_OR)
    assert groups(mol) == [("STEREO_OR", [1, 3])]

    # centers that are already grouped stay where they are
    mol = flagged(1, (OR, [1], 1))
    rdmolops.TranslateChiralFlagToStereoGroups(mol)
    assert groups(mol) == [("STEREO_OR", [1]), ("STEREO_ABSOLUTE", [3])]
    # an existing ABS group takes the new centers
    mol = flagged(0, (ABS, [1], 0))
    rdmolops.TranslateChiralFlagToStereoGroups(mol)
    assert groups(mol) == [("STEREO_ABSOLUTE", [1]), ("STEREO_AND", [3])]

    mol = flagged(None)
    rdmolops.TranslateChiralFlagToStereoGroups(mol)
    assert groups(mol) == []

    with pytest.raises(TypeError, match="did not match C\\+\\+ signature"):
        rdmolops.TranslateChiralFlagToStereoGroups(Chem.MolFromSmiles("CC"), 5)


def test_cleanup_stereo_groups_drops_unspecified_centers():
    mol = with_groups("C[C@H](F)[C@@H](Cl)Br", (ABS, [1, 3], 0))
    mol.GetAtomWithIdx(1).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rdmolops.CleanupStereoGroups(mol)
    assert [[a.GetIdx() for a in g.GetAtoms()] for g in mol.GetStereoGroups()] == [[3]]
    mol.GetAtomWithIdx(3).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rdmolops.CleanupStereoGroups(mol)
    assert len(mol.GetStereoGroups()) == 0
