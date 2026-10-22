"""rdinchi / Chem.inchi: InChIKeys against tests/data/inchi_corpus.jsonl (made
with RDKit 2026.09.1 / InChI 1.07.3 by tools/differential/make_inchi_corpus.py)
and the stubs of the not-yet-ported conversions."""
import json
from pathlib import Path

import pytest

from rdoxide import Chem
from rdoxide.Chem import inchi, rdinchi

CORPUS = Path(__file__).resolve().parents[1] / "data" / "inchi_corpus.jsonl"


def corpus():
    return [json.loads(line) for line in CORPUS.read_text().splitlines()]


def test_version():
    assert rdinchi.GetInchiVersion() == "1.07.3"
    assert Chem.GetInchiVersion() == "1.07.3"


def test_inchi_to_inchi_key_matches_rdkit():
    for rec in corpus():
        assert rdinchi.InchiToInchiKey(rec["inchi"]) == rec["key"], rec["inchi"]
        # the Python wrapper maps errors (empty key) to None
        assert Chem.InchiToInchiKey(rec["inchi"]) == (rec["key"] or None)


def test_inchi_to_inchi_key_logs_errors(capfd):
    from rdoxide import RDLogger
    RDLogger.EnableLog("rdApp.error")
    assert rdinchi.InchiToInchiKey("InChI=1S/CH4/h1H4/f/h") == ""
    assert "Invalid standard InChI string in generating InChI Key" in capfd.readouterr().err


def test_module_surface():
    assert Chem.INCHI_AVAILABLE is True and inchi.INCHI_AVAILABLE is True
    for name in ("MolToInchi", "MolFromInchi", "MolFromInchiAndAuxInfo", "MolToInchiKey"):
        assert callable(getattr(inchi, name))
    assert issubclass(Chem.InchiReadWriteError, Exception)
    for name in ("InchiToMol", "MolToInchi", "MolBlockToInchi", "InchiToInchiKey", "MolToInchiKey",
                 "GetInchiVersion"):
        assert callable(getattr(rdinchi, name))


def test_mol_to_inchi():
    m = Chem.MolFromSmiles("CC(=O)[O-]")
    inchi, code, msg, log, aux = rdinchi.MolToInchi(m)
    assert inchi == "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/p-1"
    assert code in (0, 1)
    assert aux.startswith("AuxInfo=1/1/")
    inchi, code, msg, log, aux = rdinchi.MolToInchi(m, options="/FixedH")
    assert inchi.endswith("/fC2H3O2/q-1")
    assert rdinchi.MolToInchiKey(m) == "QTBSBXVTEAMEQO-UHFFFAOYSA-M"


def test_mol_block_to_inchi():
    m = Chem.MolFromSmiles("CCO")
    inchi, code, msg, log, aux = rdinchi.MolBlockToInchi(Chem.MolToMolBlock(m))
    assert inchi == "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    assert code == 0


def test_inchi_to_mol():
    m, code, msg, log = rdinchi.InchiToMol("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", True, True)
    assert code == 0 and Chem.MolToSmiles(m) == "CCO"
    m = Chem.MolFromInchi("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1")
    assert Chem.MolToSmiles(m) == "C[C@H](N)C(=O)O"
    assert Chem.MolFromInchi("not an inchi") is None
