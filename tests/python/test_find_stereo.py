"""FindPotentialStereo / StereoInfo and the EnumerateStereoisomers module.

Expected values were produced with the reference RDKit.
"""

from rdoxide import Chem
from rdoxide.Chem import rdchem
from rdoxide.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)


def _summary(smiles, **kwargs):
    infos = Chem.FindPotentialStereo(Chem.MolFromSmiles(smiles), **kwargs)
    return [
        (str(i.type), str(i.specified), i.centeredOn, str(i.descriptor), list(i.controllingAtoms))
        for i in infos
    ]


def test_atom_and_bond():
    noatom = rdchem.StereoInfo.NOATOM
    assert _summary("C[C@H](F)C=CC") == [
        ("Atom_Tetrahedral", "Specified", 1, "Tet_CCW", [0, 2, 3]),
        ("Bond_Double", "Unspecified", 3, "NoValue", [1, noatom, 5, noatom]),
    ]


def test_ring_stereo():
    assert _summary("C[C@H]1CC[C@@H](C)CC1") == [
        ("Atom_Tetrahedral", "Specified", 1, "Tet_CW", [0, 2, 7]),
        ("Atom_Tetrahedral", "Specified", 4, "Tet_CW", [3, 5, 6]),
    ]


def test_cumulene():
    assert [s[2] for s in _summary("CC=C=CC")] == [1, 2]


def test_enums_and_setters():
    info = Chem.FindPotentialStereo(Chem.MolFromSmiles("CC(F)Cl"))[0]
    assert info.type == rdchem.StereoType.Atom_Tetrahedral
    assert info.specified == rdchem.StereoSpecified.Unspecified
    assert info.descriptor == rdchem.StereoDescriptor.NoValue
    info.specified = rdchem.StereoSpecified.Specified
    assert info.specified == rdchem.StereoSpecified.Specified


def test_clean_it_removes_spurious_stereo():
    mol = Chem.MolFromSmiles("C[C@H](C)F", sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.FindPotentialStereo(mol, cleanIt=True)
    assert mol.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED


def test_enumerate_stereoisomers():
    mol = Chem.MolFromSmiles("CC(F)C=CC")
    assert sorted(Chem.MolToSmiles(m) for m in EnumerateStereoisomers(mol)) == [
        "C/C=C/[C@@H](C)F",
        "C/C=C/[C@H](C)F",
        "C/C=C\\[C@@H](C)F",
        "C/C=C\\[C@H](C)F",
    ]
    opts = StereoEnumerationOptions(maxIsomers=2)
    assert len(list(EnumerateStereoisomers(mol, opts))) == 2
