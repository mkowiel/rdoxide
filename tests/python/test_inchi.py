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
    # like RDKit built without InChI support until MolToInchi/MolFromInchi are ported
    assert Chem.INCHI_AVAILABLE is False and inchi.INCHI_AVAILABLE is False
    assert not hasattr(inchi, "MolToInchi")
    assert issubclass(Chem.InchiReadWriteError, Exception)
    for name in ("InchiToMol", "MolToInchi", "MolBlockToInchi", "InchiToInchiKey", "MolToInchiKey",
                 "GetInchiVersion"):
        assert callable(getattr(rdinchi, name))


def test_mol_to_inchi():
    m = Chem.MolFromSmiles("CC(=O)[O-]")
    inchi, code, msg, log, aux = rdinchi.MolToInchi(m)
    assert inchi == "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/p-1"
    assert code == 0
    assert aux.startswith("AuxInfo=1/1/")
    inchi, code, msg, log, aux = rdinchi.MolToInchi(m, options="/FixedH")
    assert inchi.startswith("InChI=1/C2H3O2/")
    assert rdinchi.MolToInchiKey(m) == "QTBSBXVTEAMEQO-UHFFFAOYSA-M"


def test_mol_block_to_inchi():
    m = Chem.MolFromSmiles("CCO")
    inchi, code, msg, log, aux = rdinchi.MolBlockToInchi(Chem.MolToMolBlock(m))
    assert inchi == "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    assert code == 0


def test_inchi_to_mol_not_ported():
    with pytest.raises(NotImplementedError):
        rdinchi.InchiToMol("InChI=1S/CH4/h1H4", sanitize=True, removeHs=True)
