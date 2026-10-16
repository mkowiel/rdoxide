"""MultithreadedSDMolSupplier / MultithreadedSmilesMolSupplier (expectations from the reference RDKit)."""

from pathlib import Path

import pytest
from rdoxide import Chem

DATA = Path(__file__).parent.parent / "data" / "molfiles"
NAMES = ["48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213", "220", "229", "256"]


@pytest.mark.parametrize("threads", [1, 3, 0])
def test_sd_molecules_in_file_order(threads):
    sup = Chem.MultithreadedSDMolSupplier(str(DATA / "NCI_aids_few.sdf"), numWriterThreads=threads)
    assert [m.GetProp("_Name") for m in sup] == NAMES
    assert sup.GetLastRecordId() == 16
    assert sup.GetLastItemText().rstrip().endswith("$$$$")


def test_context_manager_closes_the_supplier():
    with Chem.MultithreadedSDMolSupplier(str(DATA / "NCI_aids_few.sdf")) as sup:
        assert next(sup).GetProp("_Name") == "48"
    with pytest.raises(StopIteration):
        next(sup)


def test_smiles_supplier_properties_and_unparsable_record():
    sup = Chem.MultithreadedSmilesMolSupplier(
        str(DATA / "fewSmi.csv"), delimiter=",", smilesColumn=1, nameColumn=0, titleLine=False, numWriterThreads=2
    )
    mols = list(sup)
    assert len(mols) == 11 and mols[-1] is None  # the empty record after the final newline
    assert [m.GetProp("_Name") for m in mols[:10]] == [str(i) for i in range(1, 11)]
    assert mols[0].GetProp("Column_2") == "34.14"


def test_smiles_without_name_column():
    sup = Chem.MultithreadedSmilesMolSupplier(str(DATA / "first_200.tpsa.csv"), ",", 0, -1)
    mol = next(sup)
    assert mol.GetProp("_Name") == "4294967295"  # upstream never counts the lines
    assert mol.GetProp("TPSA") == "34.14"


def test_property_list_switch_and_unparsable_sd_record():
    sup = Chem.MultithreadedSDMolSupplier(str(DATA / "NCI_aids_few.sdf"))
    assert sup.GetProcessPropertyLists() is True
    sup.SetProcessPropertyLists(False)
    assert sup.GetProcessPropertyLists() is False
    assert list(Chem.MultithreadedSDMolSupplier(str(DATA / "sdErrors1.sdf"))) == [None]


def test_missing_file():
    with pytest.raises(OSError):
        Chem.MultithreadedSDMolSupplier("/does/not/exist.sdf")
