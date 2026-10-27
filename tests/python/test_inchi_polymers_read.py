"""Reading polymer InChI (`InChI=1B/.../z...`) against `tests/data/inchi_polymers_read.jsonl`
(tools/make_inchi_polymer_read_corpus.py, official RDKit 2026.09.1)."""

import json
from pathlib import Path

from rdoxide import Chem
from rdoxide.Chem import rdinchi

DATA = Path(__file__).resolve().parents[1] / "data" / "inchi_polymers_read.jsonl"
CASES = [json.loads(line) for line in DATA.read_text().splitlines()]


def _run(case):
    m, ret, msg, log = rdinchi.InchiToMol(case["inchi"], False, False)
    return {"molblock": Chem.MolToMolBlock(m) if m is not None else None, "ret": ret,
            "message": msg, "log": log}


def agreement():
    return sum(
        _run(c) == {k: c[k] for k in ("molblock", "ret", "message", "log")} for c in CASES
    )


# number of agreeing cases; raise it as the polymer reading lands
FLOOR = 219


def test_polymer_inchi_read_progress_floor():
    assert agreement() >= FLOOR
