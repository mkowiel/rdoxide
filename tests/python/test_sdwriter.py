"""SDWriter. Expected values were checked against RDKit 2026.09.1."""
import io
import os

import pytest

from rdoxide import Chem

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "molfiles")


def mols():
    return [m for m in Chem.SDMolSupplier(os.path.join(DATA, "first_200.props.sdf"))][:5]


def test_roundtrip_through_file_object():
    sio = io.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols():
        w.write(m)
    assert w.NumMols() == 5
    w.flush()
    text = sio.getvalue()
    assert text.count("$$$$\n") == 5
    assert ">  <AMW>  (1) \n122.12344\n\n" in text
    back = Chem.SDMolSupplier()
    back.SetData(text)
    assert [Chem.MolToSmiles(m) for m in back] == [Chem.MolToSmiles(m) for m in mols()]


def test_props_kekulize_v3000(tmp_path):
    fn = str(tmp_path / "out.sdf")
    with Chem.SDWriter(fn) as w:
        w.SetProps(["AMW", "missing"])
        assert w.GetKekulize() and not w.GetForceV3000()
        w.SetForceV3000(True)
        w.write(mols()[0])
    text = open(fn).read()
    assert "V3000" in text
    assert "NUM_RINGS" not in text and ">  <AMW>  (1) \n" in text
    with pytest.raises(RuntimeError):
        w.write(mols()[0])


def test_get_text():
    m = Chem.MolFromSmiles("c1ccccc1O")
    m.SetProp("foo", "bar")
    m.SetIntProp("n", 3)
    text = Chem.SDWriter.GetText(m)
    assert text.endswith("M  END\n>  <foo>  \nbar\n\n>  <n>  \n3\n\n$$$$\n")
    assert Chem.SDWriter.GetText(m, molid=4).count("(5) ") == 2
    m.SetProp("bad", "a\n\nb")
    assert ">  <bad>" not in Chem.SDWriter.GetText(m)
