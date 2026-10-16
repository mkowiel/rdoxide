"""rdCIPLabeler and what uses it under the new stereo perception."""

import pytest
from rdkit import Chem, RDConfig
from rdkit.Chem import rdCIPLabeler, rdFingerprintGenerator


@pytest.fixture
def new_stereo():
    Chem.SetUseLegacyStereoPerception(False)
    try:
        yield
    finally:
        Chem.SetUseLegacyStereoPerception(True)


def labels(mol):
    atoms = {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}
    bonds = {b.GetIdx(): b.GetProp("_CIPCode") for b in mol.GetBonds() if b.HasProp("_CIPCode")}
    return atoms, bonds


def test_assign_all():
    mol = Chem.MolFromSmiles(r"C[C@H](F)C(Cl)=C/C=C\C")
    legacy = labels(mol)  # the legacy perception labels the atoms on parsing
    assert legacy[0] == {1: "S"}
    rdCIPLabeler.AssignCIPLabels(mol)
    atoms, bonds = labels(mol)
    assert atoms == legacy[0]
    assert bonds == {6: "Z"}
    assert mol.HasProp("_CIPComputed")


def test_pseudoasymmetric_and_auxiliary_labels():
    mol = Chem.MolFromSmiles("CC[C@H](C)CCCCC[C@H]1CC[C@@H](C)CC1")
    rdCIPLabeler.AssignCIPLabels(mol)
    assert labels(mol)[0] == {2: "S", 9: "s", 12: "S"}


def test_selected_atoms_and_bonds():
    mol = Chem.MolFromSmiles("C[C@H](Cl)CC[C@H](Cl)C")
    for atom in mol.GetAtoms():
        atom.ClearProp("_CIPCode")
    rdCIPLabeler.AssignCIPLabels(mol, [1], None)
    assert labels(mol)[0] == {1: "S"}
    assert not mol.HasProp("_CIPComputed")
    with pytest.raises(ValueError):
        rdCIPLabeler.AssignCIPLabels(mol, [100], None)


def test_max_recursive_iterations():
    mol = Chem.MolFromSmiles("C[C@H](F)Cl")
    with pytest.raises(RuntimeError, match="Max Iterations Exceeded"):
        rdCIPLabeler.AssignCIPLabels(mol, maxRecursiveIterations=1)
    rdCIPLabeler.AssignCIPLabels(mol, maxRecursiveIterations=1000)
    assert labels(mol)[0] == {1: "R"}


def test_find_mol_chiral_centers_new_stereo(new_stereo):
    mol = Chem.MolFromSmiles("C[C@H](F)Cl")
    assert Chem.FindMolChiralCenters(mol, useLegacyImplementation=False) == [(1, "R")]
    mol = Chem.MolFromSmiles("OC(=O)[C@H]1CC[C@@H](CC1)O[C@@H](F)Cl")
    assert Chem.FindMolChiralCenters(mol, useLegacyImplementation=False) == [
        (3, "S"),
        (6, "r"),
        (10, "S"),
    ]


def test_morgan_chirality_new_stereo(new_stereo):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, includeChirality=True)
    r = gen.GetSparseCountFingerprint(Chem.MolFromSmiles("C[C@H](F)Cl"))
    s = gen.GetSparseCountFingerprint(Chem.MolFromSmiles("C[C@@H](F)Cl"))
    assert r.GetNonzeroElements() != s.GetNonzeroElements()
    # a molecule seen several times gets the same answer
    assert r.GetNonzeroElements() == gen.GetSparseCountFingerprint(
        Chem.MolFromSmiles("C[C@H](F)Cl")
    ).GetNonzeroElements()


def test_atom_pair_chirality_new_stereo(new_stereo):
    gen = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality=True)
    r = gen.GetSparseCountFingerprint(Chem.MolFromSmiles("C[C@H](F)Cl"))
    s = gen.GetSparseCountFingerprint(Chem.MolFromSmiles("C[C@@H](F)Cl"))
    assert r.GetNonzeroElements() != s.GetNonzeroElements()


def test_morgan_double_bond_stereo_new_stereo(new_stereo):
    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, includeChirality=True, atomInvariantsGenerator=None
    )
    e = gen.GetSparseCountFingerprint(Chem.MolFromSmiles(r"C/C=C/C"))
    z = gen.GetSparseCountFingerprint(Chem.MolFromSmiles(r"C/C=C\C"))
    assert e.GetNonzeroElements() != z.GetNonzeroElements()
