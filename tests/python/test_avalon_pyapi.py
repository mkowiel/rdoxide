import json
import pathlib

import pytest

from rdoxide import Chem
from rdoxide.Avalon import pyAvalonTools as p

DATA = pathlib.Path(__file__).parent.parent / "data" / "avalon_pyapi.jsonl"
CASES = [json.loads(line) for line in DATA.read_text().splitlines()]


def _run(rec):
    kind, inp = rec["kind"], rec["input"]
    if kind == "canon_mol":
        return p.GetCanonSmiles(Chem.MolFromSmiles(inp), rec["flags"])
    if kind == "canon_mol2d":
        return p.GetCanonSmiles(Chem.MolFromMolBlock(inp), rec["flags"])
    if kind == "canon_smiles":
        return p.GetCanonSmiles(inp, True, rec["flags"])
    if kind == "canon_molblock":
        return p.GetCanonSmiles(inp, False, rec["flags"])
    if kind == "coords_smiles":
        return p.Generate2DCoords(inp, True)
    if kind == "coords_molblock":
        return p.Generate2DCoords(inp, False)
    if kind in ("coords_mol", "coords_molh"):
        m = Chem.MolFromSmiles(inp)
        if kind == "coords_molh":
            m = Chem.AddHs(m)
        # the coordinates depend on the platform's rand() in the reference; only the outcome is compared
        cid = p.Generate2DCoords(m)
        pos = m.GetConformer(cid).GetPositions()
        assert len(pos) == rec["natoms"] and abs(pos).max() < 1e3
        return {"cid": cid, "natoms": m.GetNumAtoms()}
    nb, q, fl = rec["args"]
    if kind == "fp_smiles":
        return list(p.GetAvalonFP(inp, True, nb, q, False, fl).GetOnBits())
    if kind == "fp_molblock":
        return list(p.GetAvalonFP(inp, False, nb, q, False, fl).GetOnBits())
    if kind == "words_smiles":
        return list(p.GetAvalonFPAsWords(inp, True, nb, q, False, fl))
    if kind == "count_smiles":
        cf = p.GetAvalonCountFP(inp, True, nb, q, fl)
        return {str(i): c for i, c in enumerate(cf) if c}
    raise AssertionError(kind)


def test_all_cases_match_reference():
    bad = []
    for rec in CASES:
        if rec.get("error"):
            with pytest.raises(Exception):
                _run(rec)
            continue
        got = _run(rec)
        if rec["kind"] in ("coords_mol", "coords_molh"):
            if got != {"cid": rec["cid"], "natoms": rec["natoms"]}:
                bad.append((rec["kind"], rec["input"][:60]))
        elif got != rec["out"]:
            bad.append((rec["kind"], rec["input"][:60]))
    assert not bad, f"{len(bad)} of {len(CASES)} differ, e.g. {bad[:5]}"


def test_constants_and_defaults():
    assert p.avalonSSSBits == 0x7FFF
    assert p.avalonSimilarityBits == 0xF07FFF
    assert p.GetCanonSmiles(Chem.MolFromSmiles("OCC")) == "CCO"
    assert p.GetCanonSmiles("OCC", True) == "CCO"
