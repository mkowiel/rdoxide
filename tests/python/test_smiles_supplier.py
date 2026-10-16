"""SmilesMolSupplier (behaviour checked against RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

TABLE = ("smiles id logP\nCCO e1 0.5\nc1ccccc1 b1 2.1\n# comment\n\nCC(=O)O a1 -0.3\n"
         "badsmiles x 1\nC[N+](=O)[O-] n 0\n")


def smiles(mol):
  return None if mol is None else Chem.MolToSmiles(mol)


def test_iteration_names_and_properties():
  suppl = Chem.SmilesMolSupplierFromText(TABLE)
  mols = [m for m in suppl]
  assert [smiles(m) for m in mols] == ["CCO", "c1ccccc1", "CC(=O)O", None, "C[N+](=O)[O-]"]
  assert mols[0].GetProp("_Name") == "e1"
  assert mols[0].GetProp("logP") == "0.5"
  assert mols[1].GetProp("_Name") == "b1"


def test_random_access_and_length():
  suppl = Chem.SmilesMolSupplierFromText(TABLE)
  assert len(suppl) == 5
  assert smiles(suppl[2]) == "CC(=O)O"
  assert smiles(suppl[-1]) == "C[N+](=O)[O-]"
  assert suppl[3] is None
  with pytest.raises(IndexError):
    suppl[5]
  assert suppl.GetItemText(1) == "c1ccccc1 b1 2.1"


def test_reset_and_next():
  suppl = Chem.SmilesMolSupplierFromText(TABLE)
  first = smiles(next(suppl))
  next(suppl)
  suppl.reset()
  assert smiles(next(suppl)) == first


def test_no_title_line_uses_column_names():
  suppl = Chem.SmilesMolSupplierFromText("CCO e1 0.5\nCCN e2 0.7\n", titleLine=False)
  mol = suppl[1]
  assert mol.GetProp("_Name") == "e2"
  assert mol.GetProp("Column_2") == "0.7"


def test_name_column_minus_one_names_molecules_after_their_line():
  suppl = Chem.SmilesMolSupplierFromText("smi\nCCO\nCCC\n", nameColumn=-1)
  assert [m.GetProp("_Name") for m in suppl] == ["1", "2"]


def test_custom_columns_and_delimiter():
  suppl = Chem.SmilesMolSupplierFromText("smiles\tname\tv\nCCO\te1\t3\n", delimiter="\t")
  assert suppl[0].GetProp("v") == "3"
  suppl = Chem.SmilesMolSupplierFromText("id smiles\nx1 CCO\n", smilesColumn=1, nameColumn=0)
  assert smiles(suppl[0]) == "CCO" and suppl[0].GetProp("_Name") == "x1"


def test_sanitize_flag():
  text = "s n\nCN(C)(C)(C)C b\n"
  assert Chem.SmilesMolSupplierFromText(text)[0] is None
  assert Chem.SmilesMolSupplierFromText(text, sanitize=False)[0] is not None


def test_line_endings_and_missing_final_newline():
  assert len(Chem.SmilesMolSupplierFromText("s n\r\nCCO a\r\nCCC b\r\n")) == 2
  assert len(Chem.SmilesMolSupplierFromText("s n\nCCO a\nCCC b")) == 2


@pytest.mark.parametrize("text", ["", "s n\n", "# a\n# b\n"])
def test_inputs_without_molecules(text):
  suppl = Chem.SmilesMolSupplierFromText(text)
  assert len(suppl) == 0
  # upstream quirk: iterating an empty supplier after len() raises
  with pytest.raises(RuntimeError, match="ran out of lines"):
    [m for m in suppl]


def test_set_data_and_files(tmp_path):
  suppl = Chem.SmilesMolSupplier()
  suppl.SetData("x y\nCCO a\n")
  assert len(suppl) == 1
  path = tmp_path / "a.smi"
  path.write_text(TABLE)
  with Chem.SmilesMolSupplier(str(path), delimiter=" ") as from_file:
    assert len(from_file) == 5
    assert smiles(from_file[1]) == "c1ccccc1"


def test_bad_files(tmp_path):
  with pytest.raises(OSError, match="Bad input file"):
    Chem.SmilesMolSupplier(str(tmp_path / "missing.smi"))
  (tmp_path / "empty.smi").write_text("")
  with pytest.raises(OSError, match="Invalid input file"):
    Chem.SmilesMolSupplier(str(tmp_path / "empty.smi"))
