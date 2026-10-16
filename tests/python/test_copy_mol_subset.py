"""CopyMolSubset, PathToSubmol and their option classes (expectations from RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdDepictor, rdmolops
from rdoxide.Chem.rdmolops import SubsetInfo, SubsetMethod, SubsetOptions

S = Chem.MolToSmiles


def pairs(mapping):
    return list(mapping)


def test_subset_options_defaults_and_setters():
    options = SubsetOptions()
    assert (options.sanitize, options.clearComputedProps, options.copyAsQuery) == (False, False, False)
    assert options.copyCoordinates is True
    assert options.conformerIdx == 4294967295
    assert options.method == SubsetMethod.BONDS_BETWEEN_ATOMS
    options.method = SubsetMethod.BONDS
    assert options.method == SubsetMethod.BONDS
    with pytest.raises(TypeError, match="did not match C\\+\\+ signature"):
        options.method = 3
    with pytest.raises(OverflowError):
        options.conformerIdx = -3
    assert set(SubsetMethod.names) == {"BONDS_BETWEEN_ATOMS", "BONDS"}


def test_atoms_and_bonds_form():
    mol = Chem.MolFromSmiles("CCOC")
    assert S(rdmolops.CopyMolSubset(mol, [0, 1, 2], [0, 1])) == "CCO"
    info = SubsetInfo()
    sub = rdmolops.CopyMolSubset(mol, [0, 1, 2], [0, 1], info)
    assert sub.GetNumAtoms() == 3
    assert pairs(info.atomMapping) == [(0, 0), (1, 1), (2, 2)]
    assert pairs(info.bondMapping) == [(0, 0), (1, 1)]
    # keywords and options work in every position
    info = SubsetInfo()
    sub = rdmolops.CopyMolSubset(mol=mol, atomIndices=[0, 1], bondIndices=[0], subsetInfo=info)
    assert (sub.GetNumAtoms(), pairs(info.bondMapping)) == (2, [(0, 0)])
    assert rdmolops.CopyMolSubset(mol, [0, 1, 2], [0, 1], SubsetOptions()).GetNumAtoms() == 3


def test_selecting_everything_is_a_plain_copy():
    mol = Chem.MolFromSmiles("CCOC")
    info = SubsetInfo()
    copy = rdmolops.CopyMolSubset(mol, [0, 1, 2, 3], [0, 1, 2], info)
    assert S(copy) == "CCOC"
    assert pairs(info.atomMapping) == [] and pairs(info.bondMapping) == []


def test_bad_indices_and_conflicting_bonds():
    mol = Chem.MolFromSmiles("CCOC")
    with pytest.raises(ValueError, match="subset bonds contain atoms not contained"):
        rdmolops.CopyMolSubset(mol, [0, 2], [0])
    with pytest.raises(IndexError, match="10"):
        rdmolops.CopyMolSubset(mol, [0, 10], [0])
    with pytest.raises(IndexError, match="10"):
        rdmolops.CopyMolSubset(mol, [0], [10])
    with pytest.raises(OverflowError, match="negative"):
        rdmolops.CopyMolSubset(mol, [-1], [0])
    with pytest.raises(TypeError, match="unsigned int"):
        rdmolops.CopyMolSubset(mol, "a", [0])
    with pytest.raises(AttributeError, match="__iter__"):
        rdmolops.CopyMolSubset(mol, 5)
    with pytest.raises(TypeError, match="did not match C\\+\\+ signature"):
        rdmolops.CopyMolSubset(mol)


def test_path_form():
    mol = Chem.MolFromSmiles("CCOC")
    assert rdmolops.CopyMolSubset(mol, [1, 2]).GetNumAtoms() == 2
    assert rdmolops.CopyMolSubset(mol, (1, 2, 3)).GetNumAtoms() == 3
    assert rdmolops.CopyMolSubset(mol, None).GetNumAtoms() == 0
    # indices outside the molecule are ignored
    assert rdmolops.CopyMolSubset(mol, [1, 2, 99]).GetNumAtoms() == 2

    options = SubsetOptions()
    options.method = SubsetMethod.BONDS
    info = SubsetInfo()
    sub = rdmolops.CopyMolSubset(mol, [1, 2], info, options)
    assert S(sub) == "COC"
    assert pairs(info.atomMapping) == [(1, 0), (2, 1), (3, 2)]
    assert pairs(info.bondMapping) == [(1, 0), (2, 1)]
    info = SubsetInfo()
    sub = rdmolops.CopyMolSubset(mol, path=[0, 1], subsetInfo=info, options=SubsetOptions())
    assert sub.GetNumAtoms() == 2


def test_uint_map_behaves_like_the_boost_map():
    mol = Chem.MolFromSmiles("CCOC")
    info = SubsetInfo()
    rdmolops.CopyMolSubset(mol, [0, 1, 2], info)
    mapping = info.atomMapping
    mapping[9] = 3
    mapping[5] = 6
    del mapping[5]
    del mapping[100]  # missing keys are ignored
    assert len(info.atomMapping) == 4  # writes go through to the info
    assert 1 in mapping and 7 not in mapping
    assert mapping[1] == 1
    assert list(mapping) == [(0, 0), (1, 1), (2, 2), (9, 3)]
    with pytest.raises(KeyError, match="Invalid key"):
        mapping[100]


def test_options_control_what_is_copied():
    mol = Chem.MolFromSmiles("F/C=C/C(Cl)[C@H](O)N")
    rdDepictor.Compute2DCoords(mol)
    mol.SetProp("foo", "bar")
    mol.SetProp("cp", "x", computed=True)

    sub = rdmolops.CopyMolSubset(mol, [0, 1, 2])
    assert sub.GetNumConformers() == 1
    assert sub.GetPropsAsDict(True, True)["cp"] == "x"
    assert "foo" not in sub.GetPropsAsDict(True, True)

    options = SubsetOptions()
    options.copyCoordinates = False
    assert rdmolops.CopyMolSubset(mol, [0, 1, 2], options).GetNumConformers() == 0
    options = SubsetOptions()
    options.clearComputedProps = True
    assert rdmolops.CopyMolSubset(mol, [0, 1, 2], options).GetPropsAsDict(True, True) == {}

    options = SubsetOptions()
    options.copyAsQuery = True
    query = rdmolops.CopyMolSubset(Chem.MolFromSmiles("CCOC"), [0, 1, 2], options)
    assert query.GetAtomWithIdx(0).HasQuery() and query.GetBondWithIdx(0).HasQuery()
    assert Chem.MolToSmarts(query) == "[#6]-[#6]-[#8]"

    options = SubsetOptions()
    options.sanitize = True
    with pytest.raises(Chem.AtomKekulizeException):
        rdmolops.CopyMolSubset(Chem.MolFromSmiles("c1ccccc1O", sanitize=False), [0, 1, 2], options)


def test_double_bond_stereo_follows_the_subset():
    mol = Chem.MolFromSmiles("F/C=C/C(Cl)[C@H](O)N")
    sub = rdmolops.CopyMolSubset(mol, [0, 1, 2, 3])
    assert S(sub) == "C/C=C/F"
    stereo = [(str(b.GetStereo()), list(b.GetStereoAtoms())) for b in sub.GetBonds()]
    assert stereo[1] == ("STEREOTRANS", [0, 3])
    # losing a stereo atom without a replacement drops the stereo
    assert S(rdmolops.CopyMolSubset(mol, [0, 1, 2])) == "C=CF"


def test_path_to_submol():
    mol = Chem.MolFromSmiles("OC1CC1C[C@H](F)Cl")
    rdDepictor.Compute2DCoords(mol)
    atom_map = {}
    sub = rdmolops.PathToSubmol(mol, [0, 1], False, atom_map)
    assert S(sub) == "CCO"
    assert atom_map == {0: 0, 1: 1, 2: 2}
    assert sub.GetNumConformers() == 1
    assert S(rdmolops.PathToSubmol(mol, [0, 1, 2, 3])) == "CCCCO"
    assert rdmolops.PathToSubmol(mol, []).GetNumAtoms() == 0
    atom_map = {}
    query = rdmolops.PathToSubmol(mol, [3, 4, 5], True, atom_map)
    assert Chem.MolToSmarts(query) == "[#6]-[#6]-[#6@]-[#9]"
    assert atom_map == {3: 0, 4: 1, 5: 2, 6: 3}
    prefilled = {7: 8}
    rdmolops.PathToSubmol(mol, [0], False, prefilled)
    assert prefilled == {0: 0, 1: 1}
    with pytest.raises(TypeError, match="Expecting an object of type dict"):
        rdmolops.PathToSubmol(mol, [0, 1], False, 5)
