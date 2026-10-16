"""MolToSmiles / MolFragmentToSmiles / MolToRandomSmilesVect /
SmilesWriteParams. Expected values checked against RDKit 2026.09.1 (this file
also passes with the official wheel); cases adapted from upstream
Code/GraphMol/Wrap/rough_test.py and rdkit/Chem/UnitTestSmiles.py
(BSD-3-Clause)."""

import pytest

from rdoxide import Chem


@pytest.mark.parametrize(
    "smi,expected",
    [
        ("OCC", "CCO"),
        ("c1ccccc1O", "Oc1ccccc1"),
        ("C1=CC=CC=C1", "c1ccccc1"),
        ("OC(=O)[C@@H](N)Cc1ccccc1", "N[C@@H](Cc1ccccc1)C(=O)O"),
        ("F/C=C/Cl", "F/C=C/Cl"),
        ("Cl\\C=C\\F", "F/C=C/Cl"),
        ("C[C@H]1CC[C@@H](C)CC1", "C[C@H]1CC[C@@H](C)CC1"),
        ("[Na+].[Cl-].C1CC1", "C1CC1.[Cl-].[Na+]"),
        ("[2H]C", "[2H]C"),
        ("C[N+](C)(C)C", "C[N+](C)(C)C"),
        ("[CH3:1]O", "O[CH3:1]"),
        ("c1cc[nH]c1", "c1cc[nH]c1"),
        ("C1CC2CCC1CC2", "C1CC2CCC1CC2"),
    ],
)
def test_canonical(smi, expected):
    assert Chem.MolToSmiles(Chem.MolFromSmiles(smi)) == expected


def test_options():
    m = Chem.MolFromSmiles("OC(=O)[C@@H](N)Cc1ccccc1")
    assert Chem.MolToSmiles(m, isomericSmiles=False) == "NC(Cc1ccccc1)C(=O)O"
    assert Chem.MolToSmiles(m, canonical=False) == "OC(=O)[C@@H](N)Cc1ccccc1"
    assert Chem.MolToSmiles(Chem.Mol(m), kekuleSmiles=True) == "N[C@@H](CC1=CC=CC=C1)C(=O)O"
    assert Chem.MolToSmiles(m, rootedAtAtom=6) == "c1(C[C@H](N)C(=O)O)ccccc1"
    assert Chem.MolToSmiles(m, allBondsExplicit=True) == "N-[C@@H](-C-c1:c:c:c:c:c:1)-C(=O)-O"
    assert Chem.MolToSmiles(m, allHsExplicit=True) == "[NH2][C@@H]([CH2][c]1[cH][cH][cH][cH][cH]1)[C](=[O])[OH]"


def test_params_object():
    m = Chem.MolFromSmiles("OC(=O)[C@@H](N)Cc1ccccc1")
    ps = Chem.SmilesWriteParams()
    assert ps.doIsomericSmiles and ps.canonical and ps.cleanStereo and ps.includeDativeBonds
    assert ps.rootedAtAtom == -1 and not ps.doKekule
    ps.doIsomericSmiles = False
    assert Chem.MolToSmiles(m, ps) == "NC(Cc1ccccc1)C(=O)O"
    assert Chem.MolToSmiles(m, params=ps) == "NC(Cc1ccccc1)C(=O)O"
    with pytest.raises(AttributeError):
        ps.notAField = 1


def test_output_order():
    m = Chem.MolFromSmiles("OCC")
    Chem.MolToSmiles(m)
    d = m.GetPropsAsDict(True, True)
    assert list(d["_smilesAtomOutputOrder"]) == [2, 1, 0]
    assert list(d["_smilesBondOutputOrder"]) == [1, 0]


def test_ring_closures_and_dots():
    m = Chem.MolFromSmiles("C12C3C4C1C5C2C3C45")  # cubane
    assert Chem.MolToSmiles(m) == "C12C3C4C1C1C2C3C41"
    assert Chem.MolToSmiles(Chem.MolFromSmiles("C.C.O")) == "C.C.O"


def test_dative_bonds():
    m = Chem.MolFromSmiles("[NH3]->[Pt]")
    assert Chem.MolToSmiles(m) == "[NH3]->[Pt]"
    ps = Chem.SmilesWriteParams()
    ps.includeDativeBonds = False
    assert Chem.MolToSmiles(m, ps) == "[NH3][Pt]"


def test_fragment_smiles():
    m = Chem.MolFromSmiles("OC(=O)[C@@H](N)Cc1ccccc1")
    assert Chem.MolFragmentToSmiles(m, atomsToUse=[0, 1, 2, 3]) == "CC(=O)O"
    assert Chem.MolFragmentToSmiles(m, [0, 1, 2], bondsToUse=[0]) == "O.CO"
    assert Chem.MolFragmentToSmiles(m, [0, 1, 2], atomSymbols=["[A%d]" % i for i in range(m.GetNumAtoms())]) == "[A0][A1]=[A2]"
    ps = Chem.SmilesWriteParams()
    # note: upstream rejects the positional (mol, params, atoms) form; rdoxide accepts it
    assert Chem.MolFragmentToSmiles(m, params=ps, atomsToUse=[5, 6, 7]) == "ccC"
    with pytest.raises(ValueError):
        Chem.MolFragmentToSmiles(m, atomsToUse=[])


def test_random_smiles():
    m = Chem.MolFromSmiles("OC(=O)[C@@H](N)Cc1ccccc1")
    a = Chem.MolToRandomSmilesVect(m, 5, randomSeed=0xF00D)
    b = Chem.MolToRandomSmilesVect(m, 5, randomSeed=0xF00D)
    assert a == b and len(a) == 5
    assert all(Chem.MolToSmiles(Chem.MolFromSmiles(s)) == Chem.MolToSmiles(m) for s in a)
