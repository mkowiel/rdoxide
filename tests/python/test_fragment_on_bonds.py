"""FragmentOnBonds, FragmentOnSomeBonds and FragmentOnBRICSBonds (expectations from RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdmolops

M = Chem.MolFromSmiles
S = Chem.MolToSmiles


def test_fragment_on_bonds_labels_dummies_with_atom_indices():
    res = rdmolops.FragmentOnBonds(M("CCOC(=O)c1ccccc1"), [1, 3])
    assert S(res) == "[1*]OC(=[4*])c1ccccc1.[2*]CC.[3*]=O"


def test_fragment_on_bonds_options():
    m = M("CCOC(=O)c1ccccc1")
    assert S(rdmolops.FragmentOnBonds(m, [1], dummyLabels=[(7, 8)])) == "[7*]OC(=O)c1ccccc1.[8*]CC"
    assert S(rdmolops.FragmentOnBonds(m, [1], addDummies=False)) == "CC.O=C(O)c1ccccc1"
    res = rdmolops.FragmentOnBonds(m, [1], bondTypes=[Chem.BondType.DOUBLE])
    assert sorted(str(b.GetBondType()) for b in res.GetBonds() if b.GetBeginAtom().GetAtomicNum() == 0 or b.GetEndAtom().GetAtomicNum() == 0) == ["DOUBLE", "DOUBLE"]


def test_fragment_on_bonds_cuts_per_atom():
    cuts = [0] * 5
    rdmolops.FragmentOnBonds(M("CCOCC"), [1, 2], cutsPerAtom=cuts)
    assert cuts == [0, 1, 2, 1, 0]


def test_fragment_on_bonds_keeps_stereo():
    assert S(rdmolops.FragmentOnBonds(M("C[C@H](F)Cl"), [0])) == "*[C@H](F)Cl.[1*]C"
    assert S(rdmolops.FragmentOnBonds(M("F/C=C/C(Cl)/C=C/F"), [3])) == "[3*]Cl.[4*]C(/C=C/F)/C=C/F"
    # aromatic NH keeps its hydrogen when no dummies are added
    assert S(rdmolops.FragmentOnBonds(M("c1cc[nH]c1CC"), [5], addDummies=False)) == "C.Cc1ccc[nH]1"


def test_fragment_on_bonds_errors():
    m = M("CCOC")
    with pytest.raises(ValueError, match="bondIndices contains duplicates"):
        rdmolops.FragmentOnBonds(m, [0, 0])
    with pytest.raises(ValueError, match="list element larger than allowed value"):
        rdmolops.FragmentOnBonds(m, [10])
    with pytest.raises(ValueError, match="empty bond indices"):
        rdmolops.FragmentOnBonds(m, [])
    with pytest.raises(ValueError, match="bondTypes shorter than bondIndices"):
        rdmolops.FragmentOnBonds(m, [0, 1], bondTypes=[Chem.BondType.SINGLE])
    with pytest.raises(ValueError, match="cutsPerAtom shorter than the number of atoms"):
        rdmolops.FragmentOnBonds(m, [0], cutsPerAtom=[0, 0])
    with pytest.raises(RuntimeError, match="bad dummyLabel vector"):
        rdmolops.FragmentOnBonds(m, [0, 1], dummyLabels=[(1, 2)])


def test_fragment_on_some_bonds():
    frags = rdmolops.FragmentOnSomeBonds(M("CCCCO"), [0, 1, 2, 3], 2)
    assert [S(f) for f in frags] == [
        "*C[2*].[1*]C.[1*]CCO",
        "*CC[3*].[1*]C.[2*]CO",
        "[1*]C[3*].[2*]CC.[2*]CO",
        "*CCC[4*].[1*]C.[3*]O",
        "[1*]CC[4*].[2*]CC.[3*]O",
        "[2*]C[4*].[3*]CCC.[3*]O",
    ]
    frags, cuts = rdmolops.FragmentOnSomeBonds(M("CCO"), [0, 1], 1, returnCutsPerAtom=True)
    assert [S(f) for f in frags] == ["*CO.[1*]C", "[1*]O.[2*]CC"]
    assert cuts == ((1, 1, 0), (0, 1, 1))


def test_fragment_on_some_bonds_limits():
    assert rdmolops.FragmentOnSomeBonds(M("CCC"), [0, 1], 0) == ()
    assert rdmolops.FragmentOnSomeBonds(M("CCC"), [0, 1], 3) == ()
    with pytest.raises(ValueError, match="up to 63 bonds"):
        rdmolops.FragmentOnSomeBonds(M("C" * 70), list(range(64)))


@pytest.mark.parametrize(
    "smi, expected",
    [
        ("CCOC(=O)c1ccccc1", "[1*]C([6*])=O.[16*]c1ccccc1.[3*]O[3*].[4*]CC"),
        ("C1CCCCC1NC(=O)C", "[1*]C(C)=O.[15*]C1CCCCC1.[5*]N[5*]"),
        ("c1ccccc1", "c1ccccc1"),
    ],
)
def test_fragment_on_brics_bonds(smi, expected):
    assert S(rdmolops.FragmentOnBRICSBonds(M(smi))) == expected
