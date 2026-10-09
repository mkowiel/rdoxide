"""SmilesWriter. Expected values were checked against RDKit 2026.09.1."""
import io

import pytest

from rdkit import Chem


def mols():
    ms = [Chem.MolFromSmiles(s) for s in ("OC[C@H](N)C", "c1ccccc1", "CC(=O)O.[Na+]")]
    ms[0].SetProp("_Name", "alanol")
    ms[0].SetIntProp("n", 3)
    ms[0].SetDoubleProp("d", 1.5)
    ms[1].SetProp("s", "x")
    return ms


def written(*args, props=None, **kwargs):
    sio = io.StringIO()
    with Chem.SmilesWriter(sio, *args, **kwargs) as w:
        if props is not None:
            w.SetProps(props)
        for m in mols():
            w.write(m)
        assert w.NumMols() == 3
    return sio.getvalue()


def test_default_layout():
    assert written() == "SMILES Name \nC[C@@H](N)CO alanol\nc1ccccc1 1\nCC(=O)O.[Na+] 2\n"


def test_property_columns():
    text = written(props=["n", "s", "d", "missing"])
    assert text.splitlines()[:2] == ["SMILES Name n s d missing", "C[C@@H](N)CO alanol 3  1.5 "]
    assert text.splitlines()[2] == "c1ccccc1 1  x  "


def test_options():
    assert written("\t", "N", False).startswith("C[C@@H](N)CO\talanol\n")
    assert written(" ", "").splitlines()[0] == "SMILES "
    assert written(",", "Id", True, False).splitlines()[1] == "CC(N)CO,alanol"
    assert written(" ", "Name", True, True, True).splitlines()[2] == "C1=CC=CC=C1 1"


def test_file_and_closed_writer(tmp_path):
    fn = str(tmp_path / "out.smi")
    w = Chem.SmilesWriter(fn, nameHeader="")
    w.SetProps(["s"])
    w.write(mols()[1])
    w.close()
    assert open(fn).read() == "SMILES s\nc1ccccc1 x\n"
    for call in (lambda: w.write(mols()[0]), w.flush):
        with pytest.raises(RuntimeError, match="no output stream"):
            call()
    w.close()
    with pytest.raises(OSError, match="Bad output file"):
        Chem.SmilesWriter(str(tmp_path / "missing" / "x.smi"))
