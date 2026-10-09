"""CXSMILES output (MolToCXSmiles); expectations from the reference RDKit."""
import pytest

from rdkit import Chem
from rdkit.Chem import rdmolfiles

F = rdmolfiles.CXSmilesFields


def test_enhanced_stereo_and_radicals():
    m = Chem.MolFromSmiles("C[C@@H](O)[C@H](C)F |o1:1,3|")
    assert Chem.MolToCXSmiles(m) == Chem.MolToCXSmiles(m, Chem.SmilesWriteParams())
    assert "|o1:" in Chem.MolToCXSmiles(m)
    assert "|" not in Chem.MolToCXSmiles(m, Chem.SmilesWriteParams(), int(F.CX_NONE))
    r = Chem.MolFromSmiles("[CH2]C")
    assert Chem.MolToCXSmiles(r).endswith("|^1:0|")


def test_labels_and_atom_props():
    m = Chem.MolFromSmiles("*C* |$star_e;;star_e$|")
    assert Chem.MolToCXSmiles(m).startswith("*C* |$star_e;;star_e$")
    m = Chem.MolFromSmiles("CCO")
    m.GetAtomWithIdx(0).SetProp("foo", "bar")
    assert Chem.MolToCXSmiles(m).endswith("|atomProp:0.foo.bar|")


def test_coordinates():
    m = Chem.MolFromMolBlock(Chem.MolToMolBlock(Chem.MolFromSmiles("CC")) if False else """
  test

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
""")
    assert Chem.MolToCXSmiles(m) == "CC |(0,0,;1.5,0,)|"
    assert Chem.MolToCXSmiles(m, Chem.SmilesWriteParams(), int(F.CX_ALL_BUT_COORDS)) == "CC"


def test_bool_overload_and_empty():
    m = Chem.MolFromSmiles("C[C@H](F)Cl")
    assert Chem.MolToCXSmiles(m, isomericSmiles=False) == "CC(F)Cl"
    assert Chem.MolToCXSmiles(Chem.Mol()) == ""


@pytest.mark.parametrize("cx", [
    "*CCO* |$star_e;;;;star_e$,Sg:n:1,2,3::eu:::|",
    "*C* |$star_e;;star_e$,Sg:n:1::ht:::|",
    "*CC |$star_e;;$,SgD:0:querydata:val::::|",
    "*CC(*)C(*)N* |$star_e;;;star_e;;star_e;;star_e$,Sg:n:6,1,2,4::hh:6,0:2,4:|",
])
def test_substance_groups_round_trip(cx):
    m = Chem.MolFromSmiles(cx)
    assert m is not None
    assert Chem.MolToCXSmiles(m) == cx


def test_link_nodes_round_trip():
    # the outer atoms are written when the atom has more than two neighbors
    cx = "OC1CCC(F)C1 |LN:4:1.3.3.6|"
    m = Chem.MolFromSmiles(cx)
    assert m.GetProp("_molLinkNodes") == "1 3 2 5 4 5 7"
    assert Chem.MolToCXSmiles(m) == cx
    m = Chem.MolFromSmiles("OC1CCCC1 |LN:4:1.3|")
    assert m.GetProp("_molLinkNodes") == "1 3 2 5 4 5 6"
