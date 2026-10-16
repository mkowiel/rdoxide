"""SDMolSupplier. Expected values were checked against RDKit 2026.09.1."""
import os

import pytest

from rdoxide import Chem

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "molfiles")

PROPLISTS = """
  Mrv

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <atom.iprop.foo>
7

> <atom.dprop.bar>
1.5

> <atom.prop.baz>
[x] x

$$$$
"""


def test_iterate_and_random_access():
    suppl = Chem.SDMolSupplier(os.path.join(DATA, "first_200.props.sdf"))
    assert len(suppl) == 200
    assert Chem.MolToSmiles(suppl[0]) == "CC1=CC(=O)C=CC1=O"
    assert Chem.MolToSmiles(suppl[-1]) == "CC(=O)Nc1cc(C)cc(C)c1"
    assert suppl[3].GetProp("NUM_RINGS") == "1"
    assert suppl[0].GetProp("AMW") == "122.12344"
    assert suppl.GetItemText(1).startswith("\n     RDKit          2D\n\n 20 23")
    mols = [m for m in suppl]
    assert len(mols) == 200 and all(m is not None for m in mols)
    with pytest.raises(IndexError):
        suppl[200]


def test_next_and_reset():
    suppl = Chem.SDMolSupplier(os.path.join(DATA, "first_200.props.sdf"))
    m1 = next(suppl)
    next(suppl)
    suppl.reset()
    assert Chem.MolToSmiles(next(suppl)) == Chem.MolToSmiles(m1)


def test_multiline_properties():
    m = Chem.SDMolSupplier(os.path.join(DATA, "BlankPropLines.sdf"))[0]
    assert m.GetProp("MultiLineProperty1") == "foo\nbar\n \nbaz"
    assert m.GetProp("MultiLineProperty2") == "foo\n>  <fooprop>\nbaz\n "


def test_truncated_property_name():
    suppl = Chem.SDMolSupplier(os.path.join(DATA, "esters.prop_name_trunc.sdf"))
    assert [Chem.MolToSmiles(m) for m in suppl] == ["CC(=O)[O][Y]", "CCC(=O)[O][Y]"]
    assert len(suppl) == 2


def test_set_data_and_property_lists():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(PROPLISTS)
    m = next(suppl)
    assert m.GetAtomWithIdx(0).GetPropsAsDict() == {"foo": 7, "bar": 1.5}
    assert suppl.atEnd()
    with pytest.raises(StopIteration):
        next(suppl)
    with pytest.raises(IndexError):
        suppl[5]
    suppl.SetProcessPropertyLists(False)
    assert not suppl.GetProcessPropertyLists()


def test_bad_files():
    with pytest.raises(OSError):
        Chem.SDMolSupplier("/nonexistent/file.sdf")
    with Chem.SDMolSupplier(os.path.join(DATA, "first_200.props.sdf")) as suppl:
        assert next(suppl) is not None
