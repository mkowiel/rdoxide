"""MolToMolBlock & friends. tests/data/molfiles/writer_expected.json holds
the exact output of RDKit 2026.09.1 for each input."""
import json
import os

import pytest

from rdkit import Chem

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "molfiles")
CASES = json.load(open(os.path.join(DATA, "writer_expected.json")))


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_matches_rdkit(case):
    m = Chem.MolFromMolBlock(case["input"])
    assert getattr(Chem, case["function"])(m, **case["kwargs"]) == case["expected"]


def test_writer_params():
    m = Chem.MolFromMolBlock(CASES[0]["input"])
    p = Chem.MolWriterParams()
    assert (p.includeStereo, p.kekulize, p.forceV3000, p.precision) == (True, True, False, 6)
    p.forceV3000 = True
    p.precision = 2
    block = Chem.MolToMolBlock(m, p)
    assert "V3000" in block and " 1.30 0.75 0.00 0" in block
    assert Chem.MolToMolBlock(m, forceV3000=True) == Chem.MolToV3KMolBlock(m)
    with pytest.raises(AttributeError):
        p.notAnAttribute = 1


def test_mol_file_roundtrip(tmp_path):
    m = Chem.MolFromMolBlock(CASES[0]["input"])
    fn = str(tmp_path / "x.mol")
    Chem.MolToMolFile(m, fn)
    assert Chem.MolToSmiles(Chem.MolFromMolFile(fn)) == Chem.MolToSmiles(m)
    Chem.MolToV3KMolFile(m, fn)
    assert "V3000" in open(fn).read()
    m2 = Chem.MolFromSmiles("c1ccccc1O")
    assert Chem.MolToMolBlock(m2, includeStereo=False).startswith("\n     RDKit          \n\n  7  7")
