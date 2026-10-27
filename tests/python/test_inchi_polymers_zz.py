"""Polymer units capped by two Zz atoms (CRU folding / frame shift territory) against
`tests/data/inchi_polymers_zz.jsonl` (tools/make_inchi_polymer_zz_corpus.py, official RDKit 2026.09.1)."""

import json
from pathlib import Path

import pytest

from rdoxide.Chem import rdinchi

DATA = Path(__file__).resolve().parents[1] / "data" / "inchi_polymers_zz.jsonl"
CASES = [json.loads(line) for line in DATA.read_text().splitlines()]


def _run(case, opt):
    r = rdinchi.MolBlockToInchi(case["molblock"], opt)
    return [r[0], r[1], r[2], r[3], r[4]]


def agreement(opt):
    return sum(_run(c, opt) == c["results"][opt] for c in CASES)


# floors of agreeing cases per option set; raise them as the polymer stages land
FLOORS = {"": 115, "-Polymers": 115, "-Polymers -NPZz": 115, "-Polymers -FoldCRU": 115, "-Polymers105": 115}


@pytest.mark.parametrize("opt", sorted(FLOORS))
def test_zz_progress_floor(opt):
    assert agreement(opt) >= FLOORS[opt]
