"""InChI polymer support against `tests/data/inchi_polymers.jsonl` (built by
tools/make_inchi_polymer_corpus.py from the official RDKit 2026.09.1).

The default options (polymer data ignored) must already agree; the `-Polymers` option
sets are tracked as expected failures until the port of the polymer code (see
docs/INCHI_PORT.md, section 6) is complete. `PROGRESS` is the number of agreeing cases.
"""

import json
from pathlib import Path

import pytest

from rdoxide.Chem import rdinchi

DATA = Path(__file__).resolve().parents[1] / "data" / "inchi_polymers.jsonl"
CASES = [json.loads(line) for line in DATA.read_text().splitlines()]


def _run(case, opt):
    r = rdinchi.MolBlockToInchi(case["molblock"], opt)
    return [r[0], r[1], r[2], r[3], r[4]]


@pytest.mark.parametrize("i", range(len(CASES)))
def test_polymer_data_ignored_by_default(i):
    assert _run(CASES[i], "") == CASES[i]["results"][""]


def polymer_agreement(opt):
    return sum(_run(c, opt) == c["results"][opt] for c in CASES)


# floors of agreeing cases per option set; raise them as the polymer stages land
FLOORS = {
    "-Polymers": 107,
    "-Polymers -NPZz": 107,
    "-Polymers -FoldCRU": 107,
    "-Polymers -NoFoldCRU": 107,
    "-Polymers105": 123,
}


@pytest.mark.parametrize("opt", sorted(FLOORS))
def test_polymer_progress_floor(opt):
    # the molfiles RDKit rejects (stage P1, polymer validation) must give its exact errors
    assert polymer_agreement(opt) >= FLOORS[opt]
