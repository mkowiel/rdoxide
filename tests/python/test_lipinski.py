"""Lipinski-style counts (expectations from RDKit 2026.09.1)."""

import pytest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D

E = D.NumRotatableBondsOptions


def test_simple_counts():
    m = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    assert D.CalcNumLipinskiHBA(m) == 3
    assert D.CalcNumLipinskiHBD(m) == 2
    assert D.CalcNumHBD(m) == 2
    assert D.CalcNumHBA(m) == 2
    assert D.CalcNumHeteroatoms(m) == 3
    assert D.CalcNumAmideBonds(m) == 1
    assert D.CalcNumRings(m) == 1
    assert D.CalcFractionCSP3(m) == pytest.approx(0.125)
    assert D.CalcNumAromaticRings(m) == 1
    assert D.CalcNumAromaticCarbocycles(m) == 1
    assert D.CalcNumHeterocycles(m) == 0


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # [Default, NonStrict, Strict, StrictLinkages]
        ("CC(=O)Nc1ccc(O)cc1", [1, 2, 1, 2]),
        ("CC(C)(C)C(=O)OC[C@H](N)C(=O)O", [3, 5, 3, 5]),
        ("O=C(NC)c1ccc(-c2ccccc2)cc1", [2, 3, 2, 1]),
        ("c1ccccc1-c1ccccc1C", [1, 1, 1, 0]),
    ],
)
def test_rotatable_bond_definitions(smiles, expected):
    m = Chem.MolFromSmiles(smiles)
    got = [
        D.CalcNumRotatableBonds(m),
        D.CalcNumRotatableBonds(m, E.NonStrict),
        D.CalcNumRotatableBonds(m, E.Strict),
        D.CalcNumRotatableBonds(m, E.StrictLinkages),
    ]
    assert got == expected


def test_rotatable_bond_plain_arguments():
    # a plain value goes through the bool overload first: truthiness selects Strict or NonStrict
    m = Chem.MolFromSmiles("CC(C)(C)C(=O)OC[C@H](N)C(=O)O")
    assert D.CalcNumRotatableBonds(m, True) == 3
    assert D.CalcNumRotatableBonds(m, False) == 5
    assert D.CalcNumRotatableBonds(m, 2) == 3
    assert D.CalcNumRotatableBonds(m, 0) == 5
    assert D.CalcNumRotatableBonds(m, -1) == 3
    # out-of-range enum values take the StrictLinkages branch, as upstream does
    assert D.CalcNumRotatableBonds(m, E(7)) == 5


def test_enum_values():
    assert int(E.Default) == -1
    assert int(E.NonStrict) == 0
    assert int(E.Strict) == 1
    assert int(E.StrictLinkages) == 2
    assert str(E.Strict) == "Strict"


def test_ring_families():
    m = Chem.MolFromSmiles("C12CC3CC(CC(C3)C1)C2")  # adamantane
    assert D.CalcNumRings(m) == 4
    assert D.CalcNumAliphaticRings(m) == 4
    assert D.CalcNumSaturatedCarbocycles(m) == 4
    assert D.CalcNumBridgeheadAtoms(m) == 4
    m = Chem.MolFromSmiles("c1ccc2ccccc2c1")  # naphthalene
    assert D.CalcNumRings(m) == 2
    assert D.CalcNumAromaticRings(m) == 2
    assert D.CalcNumAliphaticRings(m) == 0
    m = Chem.MolFromSmiles("C1=CC2CC1C2")  # norbornene
    assert D.CalcNumSaturatedCarbocycles(m) == 1
    assert D.CalcNumAliphaticCarbocycles(m) == 3


def test_spiro_and_bridgehead_atoms():
    m = Chem.MolFromSmiles("C1CCC2(CC1)CCCC2")
    atoms = []
    assert D.CalcNumSpiroAtoms(m, atoms) == 1
    assert atoms == [3]
    assert D.CalcNumBridgeheadAtoms(m) == 0
    m = Chem.MolFromSmiles("C1CC2CCC1C2")
    atoms = []
    assert D.CalcNumBridgeheadAtoms(m, atoms) == 2
    assert atoms == [2, 5]
    assert D.CalcNumSpiroAtoms(m) == 0


def test_heterocycle_counts():
    m = Chem.MolFromSmiles("c1ccc2c(c1)OC(=O)c1ccccc1-2")
    assert D.CalcNumHeterocycles(m) == 1
    assert D.CalcNumAromaticHeterocycles(m) == 1
    assert D.CalcNumAromaticCarbocycles(m) == 2
    assert D.CalcNumRings(m) == 3


@pytest.mark.parametrize(
    "smiles, assigned, unspecified",
    [
        ("C[C@H](N)C(=O)O", 1, 0),
        ("CC(N)C(=O)O", 1, 1),
        ("F/C=C/Cl", 0, 0),
        ("C[C@H](N)C(=O)O.O=C(O)[C@@H](N)C", 2, 0),
    ],
)
def test_stereo_centers(smiles, assigned, unspecified):
    m = Chem.MolFromSmiles(smiles)
    assert D.CalcNumAtomStereoCenters(m) == assigned
    assert D.CalcNumUnspecifiedAtomStereoCenters(m) == unspecified
    # without stereo perception done yet, the counts are the same
    mm = Chem.Mol(m)
    mm.ClearProp("_StereochemDone")
    assert D.CalcNumAtomStereoCenters(mm) == assigned
    assert D.CalcNumUnspecifiedAtomStereoCenters(mm) == unspecified


def test_versions():
    assert D._CalcNumLipinskiHBD_version == "2.0.0"
    assert D._CalcNumRotatableBonds_version == "3.2.0"
    assert D._CalcNumSpiroAtoms_version == "1.0.0"
