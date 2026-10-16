import pytest

from rdkit import Chem


@pytest.fixture
def new_stereo():
    was = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(False)
    yield
    Chem.SetUseLegacyStereoPerception(was)


def stereo_bonds(m):
    return [str(b.GetStereo()) for b in m.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE]


def test_double_bond_stereo(new_stereo):
    assert Chem.MolToSmiles(Chem.MolFromSmiles("F/C=C/Cl")) == "F/C=C/Cl"
    assert stereo_bonds(Chem.MolFromSmiles("F/C=C/Cl")) == ["STEREOTRANS"]
    assert stereo_bonds(Chem.MolFromSmiles("F/C=C\\Cl")) == ["STEREOCIS"]
    assert stereo_bonds(Chem.MolFromSmiles("FC=CCl")) == []
    assert stereo_bonds(Chem.MolFromSmiles("C/C=C/C=C/C")) == ["STEREOTRANS"] * 2


def test_tetrahedral_centers(new_stereo):
    m = Chem.MolFromSmiles("C[C@H](F)Cl")
    assert [a.GetChiralTag() for a in m.GetAtoms()][1] == Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    # not a stereocenter: the tag is removed
    m = Chem.MolFromSmiles("C[C@H](C)C")
    assert all(a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for a in m.GetAtoms())
    assert Chem.MolToSmiles(Chem.MolFromSmiles("N[C@@H](C)C(=O)O")) == "C[C@H](N)C(=O)O"


def test_find_potential_stereo_matches_assignment(new_stereo):
    m = Chem.MolFromSmiles("C[C@H](O)[C@@H](O)C")
    info = Chem.FindPotentialStereo(m)
    assert [i.centeredOn for i in info] == [1, 3]
    assert all(i.specified == Chem.StereoSpecified.Specified for i in info)


def test_assign_stereochemistry_force(new_stereo):
    m = Chem.MolFromSmiles("FC=CCl", sanitize=False)
    m.UpdatePropertyCache(strict=False)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    assert m.HasProp("_StereochemDone")
