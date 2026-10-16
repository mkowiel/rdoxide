"""MolToSmarts / MolFragmentToSmarts / Atom.GetSmarts / Bond.GetSmarts.
Expected values checked against RDKit 2026.09.1 (this file also passes with
the official wheel); cases adapted from upstream rough_test.py and
smaTest.cpp (BSD-3-Clause)."""

import pytest

from rdoxide import Chem


@pytest.mark.parametrize(
    "sma,expected",
    [
        ("[C;H2]-,=[N,O]", "[C&H2]-,=[N,O]"),
        ("[!$(C=O)]", "[!$(C=O)]"),
        ("[$(C=O),$(C#N)]", "[$(C=O),$(C#N)]"),
        ("[D{2-3}]", "[D{2-3}]"),
        ("[k{5-6}]", "[k{5-6}]"),
        ("[#6]~[#7]", "[#6]~[#7]"),
        ("C-!@C", "C-&!@C"),
        ("[c,n;H1]", "[c,n;H1]"),
        ("[2H+:3]", "[#1&2*&+:3]"),
        ("[CH2]", "[C&H2]"),
        ("[R2]", "[R2]"),
        ("[^2]", "[^2]"),
        ("[+{0-1}]", "[+{0-1}]"),
        ("C/C=C/C", "C/C=C/C"),
        ("[NH3+]->[Pt]", "[N&H3&+]->[#78]"),
    ],
)
def test_mol_to_smarts(sma, expected):
    assert Chem.MolToSmarts(Chem.MolFromSmarts(sma)) == expected


def test_smarts_from_smiles():
    m = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
    assert Chem.MolToSmarts(m) == "[#6]-[#6@H](-[#7])-[#6](=[#8])-[#8]"
    assert Chem.MolToSmarts(m, isomericSmiles=False) == "[#6]-[#6](-[#7])-[#6](=[#8])-[#8]"
    assert Chem.MolToSmarts(Chem.MolFromSmiles("c1ccccc1[Fe+2]")) == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[Fe+2]"


def test_rooted_and_params():
    m = Chem.MolFromSmarts("CCO")
    assert Chem.MolToSmarts(m, rootedAtAtom=2) == "OCC"
    ps = Chem.SmilesWriteParams()
    ps.rootedAtAtom = 1
    assert Chem.MolToSmarts(m, ps) == "C(C)O"


def test_fragment():
    m = Chem.MolFromSmarts("[C;H2]-,=[N,O]C")
    assert Chem.MolFragmentToSmarts(m, atomsToUse=[0, 1]) == "[C&H2]-,=[N,O]"
    assert Chem.MolFragmentToSmarts(m, [1, 2]) == "[N,O]C"
    assert Chem.MolFragmentToSmarts(m, [0, 2]) == "[C&H2].C"
    with pytest.raises(ValueError):
        Chem.MolFragmentToSmarts(m, [])


def test_atom_and_bond_get_smarts():
    m = Chem.MolFromSmarts("[C;H2]-,=[N,O]")
    assert [a.GetSmarts() for a in m.GetAtoms()] == ["[C&H2]", "[N,O]"]
    assert m.GetBondWithIdx(0).GetSmarts() == "-,="
    m = Chem.MolFromSmiles("c1ccccc1O")
    assert m.GetAtomWithIdx(0).GetSmarts() == "c"
    assert m.GetBondWithIdx(0).GetSmarts() == ""
    assert m.GetBondWithIdx(0).GetSmarts(allBondsExplicit=True) == ":"
    assert Chem.AtomFromSmarts("[C,N]").GetSmarts() == "[C,N]"
    assert Chem.BondFromSmarts("!@").GetSmarts() == "!@"
    with pytest.raises(RuntimeError):
        Chem.Atom(6).GetSmarts()
