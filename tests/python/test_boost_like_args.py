"""Arguments that boost::python accepts: 0/1 for bools, path-like file names,
bytes for text data."""
import pathlib

from rdoxide import Chem

SDF = """
     RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""


def test_int_flags_are_accepted():
    m = Chem.MolFromSmiles("C1=CN=CC=C1")
    m.SetProp("cprop", "foo", 1)
    assert m.GetPropsAsDict(True, True, True)["cprop"] == "foo"
    m.SetBoolProp("flag", 1, 0)
    assert m.GetBoolProp("flag") is True
    Chem.Kekulize(m, clearAromaticFlags=1)
    assert len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 6
    Chem.AssignStereochemistry(m, 1, 1, 1)


def test_sd_supplier_takes_pathlib_and_bytes(tmp_path):
    path = tmp_path / "one.sdf"
    path.write_text(SDF)
    assert len(Chem.SDMolSupplier(path)) == 1
    assert len(Chem.SDMolSupplier(str(path), 0)) == 1
    sup = Chem.SDMolSupplier()
    sup.SetData(SDF.encode())
    assert len(sup) == 1
    sup = Chem.SDMolSupplier()
    sup.SetData(SDF)
    assert len(sup) == 1


def test_smiles_supplier_takes_pathlib_and_bytes(tmp_path):
    path = tmp_path / "mols.smi"
    path.write_text("CC ethane\nCCC propane\n")
    assert len(Chem.SmilesMolSupplier(pathlib.Path(path), " ", 0, 1, False)) == 2
    sup = Chem.SmilesMolSupplier()
    sup.SetData(b"CC ethane\n", titleLine=False)
    assert len(sup) == 1


def test_sd_writer_takes_pathlib(tmp_path):
    path = tmp_path / "out.sdf"
    writer = Chem.SDWriter(path)
    writer.write(Chem.MolFromSmiles("CC"))
    writer.close()
    assert path.read_text().count("$$$$") == 1
