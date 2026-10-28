import json
import pathlib

import pytest

from rdoxide import Chem
from rdoxide.Avalon import pyAvalonTools as p

DATA = pathlib.Path(__file__).parent.parent / "data" / "avalon_pyfp.jsonl"
CASES = [json.loads(line) for line in DATA.read_text().splitlines()]


def _mol(rec):
    if rec["molblock"] is None:
        return Chem.MolFromSmiles(rec["smiles"])
    return Chem.MolFromMolBlock(rec["molblock"], removeHs=False)


def test_all_cases_match_reference():
    bad = []
    for rec in CASES:
        m = _mol(rec)
        kw = dict(nBits=rec["nBits"], isQuery=rec["isQuery"], bitFlags=rec["bitFlags"])
        got = {
            "onBits": list(p.GetAvalonFP(m, **kw).GetOnBits()),
            "counts": {i: c for i, c in enumerate(p.GetAvalonCountFP(m, **kw)) if c},
            "words": list(p.GetAvalonFPAsWords(m, **kw)),
        }
        want = {"onBits": rec["onBits"], "counts": {int(k): v for k, v in rec["counts"].items()},
                "words": rec["words"]}
        if got != want:
            bad.append((rec["smiles"], rec["kind"], kw))
    assert not bad, f"{len(bad)} of {len(CASES)} differ, e.g. {bad[:3]}"


def test_defaults_and_molblock_input():
    m = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    fp = p.GetAvalonFP(m)
    assert fp.GetNumBits() == 512
    assert fp.GetNumOnBits() > 0
    mb = Chem.MolToMolBlock(m)
    assert list(p.GetAvalonFP(mb, False).GetOnBits()) == list(fp.GetOnBits())
    assert p.GetAvalonFP("CCO", True).GetNumOnBits() > 0
    with pytest.raises(TypeError):
        p.GetAvalonFP("CCO")
