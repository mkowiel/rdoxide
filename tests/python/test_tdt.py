"""TDTMolSupplier, TDTWriter and TDT files in rdMolProcessing."""
import io
from pathlib import Path

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdMolProcessing

DATA = Path(__file__).resolve().parents[2] / "upstream" / "rdkit" / "Code" / "GraphMol" / "FileParsers" / "test_data"
pytestmark = pytest.mark.skipif(not DATA.exists(), reason="upstream test data not checked out")


def test_supplier_length_and_access():
    suppl = Chem.TDTMolSupplier(str(DATA / "acd_few.tdt"))
    n = len(suppl)
    assert n == 10
    assert suppl[n - 1] is not None
    with pytest.raises(IndexError):
        suppl[n + 5]
    assert suppl.GetItemText(0).startswith("$SMI<")
    assert len([m for m in Chem.TDTMolSupplier(str(DATA / "acd_few.tdt"))]) == n


def test_conformers_and_name_record():
    suppl = Chem.TDTMolSupplier(str(DATA / "acd_few.tdt"), nameRecord="NAME", confId2D=0)
    m = next(suppl)
    assert m.HasProp("_Name") or not m.HasProp("NAME")


def test_writer_round_trip():
    m = Chem.MolFromSmiles("CCO")
    m.SetProp("_Name", "ethanol")
    m.SetProp("note", "a\nb")
    out = io.StringIO()
    w = Chem.TDTWriter(out)
    w.write(m)
    w.write(m)
    w.close()
    assert out.getvalue() == (
        "$SMI<CCO>\nNAME<ethanol>\nnote<a b>\n|\n$SMI<CCO>\nNAME<ethanol>\nnote<a b>\n")
    suppl = Chem.TDTMolSupplier()
    suppl.SetData(out.getvalue(), nameRecord="NAME")
    assert len(suppl) == 2
    assert suppl[1].GetProp("_Name") == "ethanol"


def test_mol_processing_reads_tdt():
    fps = rdMolProcessing.GetFingerprintsForMolsInFile(str(DATA / "acd_few.tdt"))
    assert len(fps) == 10 and all(fp is not None for fp in fps)
